#include "Genome.hpp"
#include "utils.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <iostream>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <string>

#define SAHAP_GENOME_DEBUG 0

enum _objectives        { OBJ_NONE,   MEC,   Poisson };
const char *objName[] = {"OBJ_NONE", "MEC", "Poisson"};
#ifndef OBJECTIVE
#define OBJECTIVE MEC // choices for now are MEC and Poisson
#endif
#if (OBJECTIVE != MEC && OBJECTIVE != Poisson)
#error "invalid objective"
#endif

// How is the temperature schedule adjusted to be dynamic?
enum _schedule            { SCHED_NONE,   RETREAT,   Betz };
const char *schedName[] = {"SCHED_NONE", "RETREAT", "Betz"};
#ifndef SCHEDULE
#define SCHEDULE RETREAT // choices for now are RETREAT and Betz
#endif
#if (SCHEDULE != RETREAT && SCHEDULE != Betz)
#error "invalid schedule"
#endif


using namespace std;

namespace SAHap {

Genome::Genome(InputFile file)
{
	long int seed = GetFancySeed(true);
	cout << "Genome seed " << seed << endl;
	this->randomEngine = mt19937(seed);
	this->file = file;
	this->haplotypes = vector<Haplotype>(this->file.ploidy, Haplotype(this->file.index.size()));
	this->shuffle();
}

Genome::~Genome() {
}

dnacnt_t Genome::mec() {
	dnaweight_t out = 0;
	for (size_t i = 0; i < haplotypes.size(); i++) {
		out = out + haplotypes[i].mec();
	}

	return out;
}

double Genome::mecScore() {
	double maxMec = this->haplotypes.size() * this->haplotypes[0].size() * 60;
	return this->mec() / maxMec;
}

double Genome::totalCoverage() {
    double coverage = 0;
    for (size_t i = 0; i < haplotypes.size(); i++)
	coverage += this->haplotypes[i].meanCoverage();
    return coverage;
}

double Genome::siteCostScore() {
	double out = 0;

	for (size_t i = 0; i < haplotypes.size(); i++) {
		// -logT_p(lambda_i, k_i)
#if OBJECTIVE == MEC
		out = out + haplotypes[i].mec();
#elif OBJECTIVE == Poisson
		out = out + haplotypes[i].siteCost();
#else
#error "No objective chosen"
#endif

	}

	// cout << "siteCost: " << out << endl;

	double maxCost = this->haplotypes.size() * this->haplotypes[0].size();
	return out / maxCost;
}

double Genome::score(dnaweight_t mec) {
	double maxMec = this->haplotypes.size() * this->haplotypes[0].size();
	return mec / maxMec;
}

double Genome::score() {
	return this->score(this->mec());
}

void Genome::shuffle() {
	uniform_int_distribution<size_t> distribution(0, this->haplotypes.size() - 1);

	this->maxIterations = this->haplotypes[0].size() * 100;

	if (this->initialized) {
		auto ploidy = this->haplotypes.size();
		auto length = this->haplotypes[0].size();
		this->haplotypes.clear();
		this->haplotypes = vector<Haplotype>(ploidy, Haplotype(length));
	}

	for (auto& r : this->file.reads) {
		this->haplotypes[distribution(this->randomEngine)].add(&r);
	}

	this->initialized = true;

	// for (auto& ch : this->haplotypes) {
	// 	cout << ch << endl;
	// }
}

double Genome::acceptance(double newScore, double curScore) {
	if (newScore < curScore) return 1;
	if (this->t == 0) return 0;
	double energyDiff = curScore - newScore;
	// cout << "Acceptance(" << energyDiff << ") = " << exp(energyDiff / this->t) << endl;

	// cout << "newScore: " << newScore << ", energyDiff: " << energyDiff << ", acceptance: " << exp(energyDiff / this->t) << endl;
	return exp(energyDiff / this->t);
}

bool Genome::done() {
	// Are we tired yet?
	// return this->pBad.getAverage() == 0 || this->t <= 1e-5;
	return this->curIteration >= this->maxIterations;
}

void Genome::setParameters(double tInitial, double tEnd, iteration_t maxIterations) {
	if (this->curIteration) {
		// Probably warn the user in some way that the simulation is in process
	}

	this->tInitial = tInitial;
	// this->tEnd = tEnd;
	this->tDecay = -log(tEnd / this->tInitial);
	this->maxIterations = maxIterations;

	cout << "decay is " << this->tDecay << endl;
}

void Genome::setTemperature(double t) {
	this->t = t;
}

void Genome::move() {
	// Perform a random move, saving enough information so we can revert later

	auto ploidy = this->haplotypes.size();
	size_t moveFrom = rand() % ploidy;
	size_t moveTo;

	while (!this->haplotypes[moveFrom].readSize()) {
		if (ploidy == 2) {
			moveFrom = !moveFrom;
		} else {
			moveFrom = rand() % ploidy;
		}
	}

	if (ploidy == 2) {
		moveTo = !moveFrom;
	} else {
		size_t moveOffset = rand() % (ploidy - 1);
		moveTo = moveFrom + moveOffset + 1;
	}

	Read * r = this->haplotypes[moveFrom].pick(this->randomEngine);
#if SAHAP_GENOME_DEBUG
	if (!r) {
		cerr << "DEBUG: Haplotype " << moveFrom << " has no reads remaining" << endl;
		std::raise(SIGINT);
		r = this->haplotypes[moveFrom].pick(this->randomEngine);
	}
#endif
	this->haplotypes[moveTo].add(r);
	this->haplotypes[moveFrom].remove(r);

	this->lastMove.from = moveFrom;
	this->lastMove.to = moveTo;
	this->lastMove.read = r;
}

void Genome::revertMove() {
	const auto& move = this->lastMove;
	this->haplotypes[move.to].remove(move.read);
	this->haplotypes[move.from].add(move.read);
}

void Genome::iteration() {
	// Run an iteration
	auto oldScore = this->siteCostScore();
	this->move();
	auto newScore = this->siteCostScore();

	uniform_real_distribution<double> distribution(0, 1);
	double chanceToKeep = this->acceptance(newScore, oldScore);
	double randomIndex = distribution(this->randomEngine);

	bool isGood = oldScore > newScore;
	bool accept = randomIndex <= chanceToKeep;
	assert(!isGood || accept);

	if (isGood || accept) {
		// always accept
	} else {
		// reject
		this->revertMove();
	}

	this->fAccept.record(isGood);

	if (!(isGood || oldScore == newScore)) {
		this->totalBad++;
		if (accept) {
			this->totalBadAccepted++;
		}
		this->pBad.record(chanceToKeep);
	}

	/*
	uniform_real_distribution<double> d(0, 3);
	if (d(this->randomEngine) <= 1) {
		this->pBad.record(false);
	} else {
		this->pBad.record(true);
	}
	*/
}

double Genome::getTemperature(iteration_t iteration) {
	double s = iteration / (double)this->maxIterations;
	double temp = (double)this->tInitial * exp(-this->tDecay * s);

	// cout << iteration << ": " << temp << endl;
	return temp;
}

void Genome::optimize(bool debug) {
	int TARGET_MEC = this->haplotypes[0].size() * this->totalCoverage() * READ_ERROR_RATE;
	// Reset state
	this->t = this->tInitial;
	this->pBad.total = this->fAccept.total = 0;
	this->pBad.pos = this->fAccept.pos = 0;
	this->pBad.sum = this->fAccept.sum = 0;

	this->totalBad = 0;
	this->totalBadAccepted = this->totalGood = 0;

	auto start_time = duration_cast<seconds>(system_clock::now().time_since_epoch());
	assert(this->haplotypes.size() == 2); // otherwise need to change a few things below that assume only 0 and 1 exist.
	assert(this->haplotypes[0].size() == this->haplotypes[1].size());
	printf("Performing %ld meta-iterations of %d each using schedule %s,\n",
	    (long)(this->maxIterations/META_ITER), META_ITER, schedName[SCHEDULE]);
	printf("optimizing objective %s across %lu sites with total coverage %g, target MEC %d\n",
	    objName[OBJECTIVE], this->haplotypes[0].size(), this->totalCoverage(), TARGET_MEC);
	while (!this->done()) {
		this->t = this->getTemperature(this->curIteration);
		double fracTime = this->curIteration*1.0/this->maxIterations;
		double pBad = this->pBad.getAverage();
		int MEC = (int)this->mec();
#if SCHEDULE==RETREAT
		static int when;
		++when;
		double retreat = 0.0; // percent
		if(when % (REPORT_INTERVAL/10) == 0) {
		    if(((fracTime>0.3||pBad<.2) && MEC > 8*TARGET_MEC) || ((fracTime>0.5||pBad<.1) && MEC > 4*TARGET_MEC))
			retreat = 0.01;
		    if((fracTime>0.95) && MEC > 1.3*TARGET_MEC) retreat = fracTime; // 100% retreat
		    if(retreat) {
			cout << "Retreat " << 100*retreat << "%, from " << this->curIteration * 100.0 / this->maxIterations;
			this->curIteration -= retreat * this->maxIterations;
			cout << "% to "  << this->curIteration * 100.0 / this->maxIterations << "%" << endl;
		    }
		}
#elif SCHEDULE==Betz
		static double computedTdecay, LOWER=atof(getenv("LOWER")),
		    ACCEPT_TARGET=atof(getenv("ACCEPT_TARGET")), PBAD_TARGET=atof(getenv("PBAD_TARGET"));
		if(!computedTdecay){
		    computedTdecay = this->tDecay;
		    printf("Betz values: LOWER %g ACCEPT_TARGET %g PBAD_TARGET %g\n",LOWER,ACCEPT_TARGET,PBAD_TARGET);
		}
		if(this->curIteration > this->fAccept.LENGTH) {
		    this->tDecay = computedTdecay * (LOWER +
			min(fabs(ACCEPT_TARGET-this->fAccept.getAverage()), fabs(PBAD_TARGET-this->pBad.getAverage())));
		}
#endif
		this->iteration();
		this->curIteration++;

		if (debug && curIteration % REPORT_INTERVAL == 0) {
			auto now_time = duration_cast<seconds>(system_clock::now().time_since_epoch());
			printf("%dk (%.4f%%,%ds)  T %.3g  fA %.3g  pBad %.3g  MEC %d", (int)this->curIteration/1000,
			    this->curIteration*100.0/this->maxIterations, (int)(now_time - start_time).count(),
			    this->t, this->fAccept.getAverage(), this->pBad.getAverage(), (int)this->mec());
			if (this->file.hasGroundTruth) {
				auto gt = this->compareGroundTruth();
				int hapSize0=this->haplotypes[0].size(),hapSize1=this->haplotypes[1].size();
				assert(hapSize0==hapSize1); // don't multiply by this->haplotypes.size()
				double he = (double)gt / this->haplotypes[0].size();
				printf("  ( Err_vs_truth %d Err_Pct %g%% [%d %d])", (int)gt, 100*he,hapSize0,hapSize1);
				if(fracTime > .99 && he > 0.01) printf("Fail!");
			}
			printf("\n");
			if(fracTime > 0.5 && pBad < 0.01 && MEC <= TARGET_MEC) {
			    this->curIteration = this->maxIterations; // basically done
			    printf("Good enough\n");
			}
		}
	}
	printf("Finished optimizing %d sites using %s cost function\n", (int)this->haplotypes[0].size(), objName[OBJECTIVE]);
}

double Genome::findPbad(double temperature, iteration_t iterations) {
	this->shuffle();

	this->t = temperature;
	this->pBad.total = 0;
	this->pBad.pos = 0;
	this->pBad.sum = 0;

	this->totalBad = 0;
	this->totalBadAccepted = 0;

	for (iteration_t i = 0; i < iterations; ++i) {
		this->iteration();
		// cout << "[simann] temp=" << this->t << ", mec=" << this->mec() << endl;
	}

	cout << "temperature " << temperature << " gives Pbad " << this->pBad.getAverage() << endl;

	return this->pBad.getAverage();
}

#define TARGET_PBAD_START 0.8
#define TARGET_PBAD_END 1e-3
void Genome::autoSchedule(iteration_t iterations) {
	cout << "Finding optimal temperature schedule..." << endl;

	double tInitial = 1;
	while (this->findPbad(tInitial) < TARGET_PBAD_START) tInitial *= 2;
	while (this->findPbad(tInitial) > TARGET_PBAD_START) tInitial /= 2;
	while (this->findPbad(tInitial) < TARGET_PBAD_START) tInitial *= 1.2;
	cout << "tInitial found: " << tInitial << endl;
	double tEnd = tInitial;
	while (this->findPbad(tEnd) > TARGET_PBAD_END) tEnd /= 2;
	while (this->findPbad(tEnd) < TARGET_PBAD_END) tEnd *= 1.2;

	cout << "tInitial = " << tInitial << ", tEnd = " << tEnd << endl;

	this->setParameters(tInitial, tEnd, iterations);
}

void Genome::PbadBuffer::record(double acceptance) {
	size_t next = this->pos == LENGTH - 1 ? 0 : this->pos + 1;
	assert(acceptance >= 0.0 && acceptance <= 1.0);
	if (this->total == LENGTH) {
		this->sum -= this->buffer[next]; // the next one is the first one
	} else {
		this->total++;
	}
	this->pos = next;
	this->buffer[next] = acceptance;
	this->sum += acceptance;
	if(this->sum < 0) this->sum = 0;
}

double Genome::PbadBuffer::getAverage() {
	return this->sum / (double)this->total;
}

void Genome::AcceptBuffer::record(char good) {
	size_t next = this->pos == LENGTH - 1 ? 0 : this->pos + 1;
	assert(good == 0 || good == 1);
	if (this->total == LENGTH) {
		this->sum -= this->buffer[next]; // the next one is the first one
	} else {
		this->total++;
	}
	this->pos = next;
	this->buffer[next] = good;
	this->sum += good;
}

double Genome::AcceptBuffer::getAverage() {
	return this->sum / (double)this->total;
}

dnacnt_t Genome::compareGroundTruth() {
	auto l00 = this->compareGroundTruth(haplotypes[0], file.groundTruth[0]);
	auto l11 = this->compareGroundTruth(haplotypes[1], file.groundTruth[1]);
	auto la = l00 + l11;

	auto l01 = this->compareGroundTruth(haplotypes[0], file.groundTruth[1]);
	auto l10 = this->compareGroundTruth(haplotypes[1], file.groundTruth[0]);
	auto lb = l01 + l10;

	return la < lb ? la : lb;
}

dnacnt_t Genome::compareGroundTruth(const Haplotype& ch, const vector<Allele>& truth) {
	dnacnt_t loss = 0;
	for (size_t i = 0; i < ch.size(); ++i) {
		if (ch.solution[i] != truth[i]) {
			loss++;
		}
	}
	return loss;
}

ostream& operator << (ostream& stream, const Genome& ge) {
	for (size_t i = 0; i < ge.haplotypes.size(); ++i) {
		for (size_t j = 0; j < ge.haplotypes[i].size(); ++j) {
			stream << ge.haplotypes[i].solution[j];
		}
		stream << endl;
	}
	return stream;
}

/*
void Genome::optimize(double temp, double minTemp, double decreaseFactor, int numiters) {
	bool just_reverted = false;
	auto ploidy = this->haplotypes.size();
	// main sim_annealing loop
	auto currentMEC = this->mec();

	while (temp > minTemp) {
		for (int i = 0; i < numiters; i++) {
			// calculate the MEC for the current haplotypes (unless we have just reverted, in which case we do not need to)
			if (!just_reverted) {
				currentMEC = this->mec();
			}
			// choose a pair read to remove and add to another chromosome
			int moveFrom = rand()%ploidy;
			ReadPair * rp = haplotypes[moveFrom].pick(this->randomEngine);
			haplotypes[moveFrom].remove(rp);
			int moveOffset = rand() % ploidy;
			int moveTo = (moveFrom + moveOffset) % ploidy;
			haplotypes[moveTo].add(rp);
			// now, calculate the new MEC score
			auto newMEC = this->mec();
			// calculate the score of the old solution and the new solution
			auto scoreOfCurrent = score(temp, currentMEC);
			auto scoreOfNew = score(temp, newMEC);
			// use the scores to caculate the chance we keep this new solution, if it is bad
			double chanceToKeep = acceptance(scoreOfCurrent, scoreOfNew, temp);
			double randomIndex = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
			//if the solution is bad, AND we do not choose a double below the chance to keep, we revert to the old solution
			if (chanceToKeep <= randomIndex && newMEC > currentMEC) {
				haplotypes[moveTo].remove(rp);
				haplotypes[moveFrom].add(rp);
			}

			cout << "[simann] temp=" << temp << ", mec=" << newMEC << endl;
		}
		temp = temp * decreaseFactor;
	}
}
*/

}
