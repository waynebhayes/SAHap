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

enum _objectives        { OBJ_NONE,   OBJ_MEC,   OBJ_Poisson };
const char *objName[] = {"OBJ_NONE",     "MEC",     "Poisson"};
#ifndef OBJECTIVE
#define OBJECTIVE OBJ_MEC // choices for now are MEC and Poisson
#endif
#if (OBJECTIVE != OBJ_MEC && OBJECTIVE != OBJ_Poisson)
#error "invalid objective"
#endif

#define TARGET_PBAD_START 0.85
#define TARGET_PBAD_END 1e-3
// How is the temperature schedule adjusted to be dynamic?
enum _schedule            { SCHED_NONE,   SCHED_RETREAT,   SCHED_Betz };
const char *schedName[] = {"SCHED_NONE",       "RETREAT",       "Betz"};
#ifndef SCHEDULE
#define SCHEDULE SCHED_RETREAT // choices for now are SCHED_RETREAT and SCHED_Betz
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
	this->range.start = 0;
	this->range.end = this->file.index.size();
	this->total_sites = this->file.index.size();
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

dnacnt_t Genome::pmec() { //Partial MEC
	dnaweight_t out = 0;
	for (size_t i = 0; i < haplotypes.size(); i++) {
		out += haplotypes[i].windowMec();
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

double Genome::totalWindowCoverage() {
	double coverage = 0;
	for (size_t i = 0; i < this->haplotypes.size(); i++)
		coverage += this->haplotypes[i].windowMeanCoverage();
	return coverage;
}

double Genome::siteCostScore() {
	double out = 0;

	for (size_t i = 0; i < haplotypes.size(); i++) {
		// -logT_p(lambda_i, k_i)
#if OBJECTIVE == OBJ_MEC
		out = out + haplotypes[i].windowMec();
#elif OBJECTIVE == OBJ_Poisson
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

	if (isGood) {
		// always accept
	} else if (!isGood && !accept) {
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

double Genome::fracTime() {
    return this->curIteration*1.0/this->maxIterations;
}

void Genome::ResetBuffers() {
    for (iteration_t i = 0; i < this->fAccept.LENGTH; ++i) this->fAccept.buffer[i]=0;
    for (iteration_t i = 0; i < this->pBad.LENGTH; ++i) this->pBad.buffer[i]=0;
    this->pBad.len = this->fAccept.len = 0;
    this->pBad.pos = this->fAccept.pos = 0;
    this->pBad.sum = this->fAccept.sum = 0;
    this->totalBad = 0;
    this->totalBadAccepted = this->totalGood = 0;
}

void Genome::optimize(bool debug) {
	// unsigned int TARGET_MEC = 0;//this->haplotypes[0].size() * this->totalCoverage() * READ_ERROR_RATE;
	// Reset state
	this->t = this->tInitial;
	ResetBuffers();

	unsigned WINDOW_SIZE = 200;
	unsigned INCREMENTS = WINDOW_SIZE / 2;
	double ERROR = 0.02;

	range.end = WINDOW_SIZE;
	
	for (auto& haplotype : this->haplotypes) {
		haplotype.initializeWindow(WINDOW_SIZE, INCREMENTS);
	}

	// Target MEC for the Window
	double PTARGET_MEC = totalWindowCoverage() * ERROR;

	auto start_time = duration_cast<seconds>(system_clock::now().time_since_epoch());
	assert(this->haplotypes.size() == 2); // otherwise need to change a few things below that assume only 0 and 1 exist.
	assert(this->haplotypes[0].size() == this->haplotypes[1].size());
	printf("Performing %ld meta-iterations of %d each using schedule %s,\n",
	    (long)(this->maxIterations/META_ITER), META_ITER, schedName[SCHEDULE]);
	printf("optimizing objective %s across %lu sites with total coverage %g, target MEC %d\n",
	    objName[OBJECTIVE], this->haplotypes[0].size(), this->totalCoverage(), PTARGET_MEC);

	
	int cpuSeconds = 0;
	int tmp = 0;
	
	while (!this->done()) {
		this->t = this->getTemperature(this->curIteration);
		this->iteration();
		this->curIteration++;
		double pBad = this->pBad.getAverage();
		DynamicSchedule(pBad, PTARGET_MEC);
		if (debug && curIteration % REPORT_INTERVAL == 0) {
			auto now_time = duration_cast<seconds>(system_clock::now().time_since_epoch());
			cpuSeconds = (int)(now_time - start_time).count();
			Report(cpuSeconds);
			// if(fracTime() > 0.5 && pBad < 0.01 && mec() <= TARGET_MEC) {
			//     printf("Exiting early because MEC %lu reached target\n", mec());
			//     this->curIteration = this->maxIterations; // basically done
			// }
		}
		if (done()){//} || (pmec() <= PTARGET_MEC)){
			range.start += INCREMENTS;
			range.end += INCREMENTS;
			curIteration = 0;
			tmp = cpuSeconds;

			if (range.start > total_sites - WINDOW_SIZE){
				break;
			}
			
			for (auto& haplotype : haplotypes) {
				haplotype.incrementWindow();
			}
			PTARGET_MEC = totalWindowCoverage() * ERROR;
		}
		if (cpuSeconds  > tmp + 300){ // Avoid getting stuck on one part
			break;
		}
	}
	Report(cpuSeconds, true);
	printf("Finished optimizing %d sites using %s cost function\n", (int)this->haplotypes[0].size(), objName[OBJECTIVE]);
	cout << "MEC: " << (int)mec() << endl;
}



void Genome::Report(int cpuSeconds, bool final) {
    printf("%2dk (%.1f%%,%ds)  T %.3f  fA %.3f  pBad %.3f  MEC %5d", (int)this->curIteration/1000, (100*fracTime()),
	cpuSeconds, this->t, this->fAccept.getAverage(), this->pBad.getAverage(), (int)this->pmec());
    if (this->file.hasGroundTruth) {
	    auto gt = this->compareGroundTruth();
	    int hapSize0=this->haplotypes[0].size(),hapSize1=this->haplotypes[1].size();
	    assert(hapSize0==hapSize1); // don't multiply by this->haplotypes.size()
	    double he = (double)gt / (this->haplotypes[0].size() * haplotypes.size());
	    printf("  ( Err_vs_truth %5d Err_Pct %.2f%% [%d %d])", (int)gt, 100*he,hapSize0,hapSize1);
	    if(final) {
		printf("\nEnding ground truth ");
		if(he > 1.3 * READ_ERROR_RATE) printf("Fail!");
		else printf("Good enough");
	    }
    }
	cout << " " << range.start << "->" << range.end;
    printf("\n");
}

double Genome::findPbad(double temperature, iteration_t iterations) {
	this->shuffle();

	this->t = temperature;
	this->ResetBuffers();

	for (iteration_t i = 0; i < iterations; ++i) {
		this->iteration();
		// cout << "[simann] temp=" << this->t << ", mec=" << this->mec() << endl;
	}

	cout << "temperature " << temperature << " gives Pbad " << this->pBad.getAverage() << endl;

	return this->pBad.getAverage();
}

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
	if (this->len == LENGTH) {
		this->sum -= this->buffer[next]; // the next one is the first one
	} else {
		this->len++;
	}
	this->pos = next;
	this->buffer[next] = acceptance;
	this->sum += acceptance;
	if(this->sum < 0) this->sum = 0;
}

double Genome::PbadBuffer::getAverage() {
	if(this->len == 0) return 0.5;
	return this->sum / (double)this->len;
}

void Genome::AcceptBuffer::record(char good) {
	size_t next = this->pos == LENGTH - 1 ? 0 : this->pos + 1;
	assert(good == 0 || good == 1);
	if (this->len == LENGTH) {
		this->sum -= this->buffer[next]; // the next one is the first one
	} else {
		this->len++;
	}
	this->pos = next;
	this->buffer[next] = good;
	this->sum += good;
}

double Genome::AcceptBuffer::getAverage() {
	if(this->len == 0 ) return 0.5;
	return this->sum / (double)this->len;
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


void Genome::DynamicSchedule(double pBad, double TARGET_MEC)
{
    //printf("pBad %ld %g\n", this->pBad.len, pBad);
#if SCHEDULE==RETREAT
/*
This schedule seems to be about right... up to 1500 SNPs. Then it fails.  Observations:
- There are a lot of magic numbers here. I think it may be possible to simplify. Ideas:
- it looks like maintaing a ration of about 10x as many 1% retreats as 95% retreats; automating that may help.
- the locations in time (fracTime and pBad) do NOT seem very flexible; you can play with what tolerances
you want at those times/pBad's, but if you shift the times or pBads too much, it all falls apart.
- decreasing the threshold (ie., tightening the tolerance) on the 1% retreats can reduce the 95% retreats,
but if you tighten it too much it never gets past 30% fracTime.
- having the first threhold too loose results in way to many 95% retreats. (ie complete restart)
- oddly, increasing the number of meta-iterations on the command line does not always help; then the retreats
are all mis-calibrated, and we never get past 30% fracTime. So again, something to automatically calibrate
them seems necessary.
- there doesn't seem to be much of a difference between optimizing MEC vs Possion. Poisson is a bit slower but
doesn't seem to give any better aswers.
- the dead-ends seem REALLY tight: sometimes you can still have 90% distance from ground truth even with the MEC
being less than a factor for two from the TARGE_MEC!!! That seems *really* strange. Is that even possible?
Perhaps the data is fucked up? It just seems so unlikely the MEC can be almost to the target while we're
still TOTALLY wrong compared to ground truth. For example, optimizing 2100SNP case, target was 614:
598k (93.4%,74s)  T 0.000  fA 0.004  pBad 0.001  MEC 850 ( Err_vs_truth 1955 Err_Pct 93.10% [2100 2100])
599k (93.6%,74s)  T 0.000  fA 0.000  pBad 0.001  MEC 851 ( Err_vs_truth 1954 Err_Pct 93.05% [2100 2100])
600k (93.8%,74s)  T 0.000  fA 0.000  pBad 0.001  MEC 852 ( Err_vs_truth 1956 Err_Pct 93.14% [2100 2100])
601k (93.9%,74s)  T 0.000  fA 0.000  pBad 0.001  MEC 852 ( Err_vs_truth 1956 Err_Pct 93.14% [2100 2100])
602k (94.1%,74s)  T 0.000  fA 0.001  pBad 0.001  MEC 852 ( Err_vs_truth 1954 Err_Pct 93.05% [2100 2100])
Retreat 94% from 94.0625% to 0.0625% because MEC is 852, too big by factor of 1.38762
We're only a factor of 1.387 above the target MEC and we're still 93% wrong compared to the ground truth??? WTF?

- again with ground truth... sometimes we converge to exactly the right ground truth yet the MEC is still large.
Again, how is that possible?  For example, here's a sequence of final lines to a good solution: (target was 268)
49k (61.3%,159s)  T 0.002  fA 0.009  pBad 0.011  MEC   264  ( Err_vs_truth     0 Err_Pct 0.00% [1000 1000])
50k (62.5%,159s)  T 0.002  fA 0.015  pBad 0.015  MEC   272  ( Err_vs_truth     0 Err_Pct 0.00% [1000 1000])
51k (63.7%,159s)  T 0.002  fA 0.011  pBad 0.009  MEC   260  ( Err_vs_truth     0 Err_Pct 0.00% [1000 1000])
How can the MEC continue to decrease while the Err_vs_truth remains exactly zero? No fucking way. There is a
unique solution---all the reads on the correct side. That solution has a fixed MEC. Can't have one change while
the other remain constant.
*/
    double retreat = 0.0;
    static double prev_retreat_frac;
    int num_meta_iters = this->maxIterations/META_ITER;
#define SMALL_RETREAT 0.01 // Let it grow with number of meta-iters? (0.01*(1+2*log(num_meta_iters)))
#define FULL_RETREAT 0.94 // this needs to be less than (1-(REPORT_INTERVAL/2)) from the next line
    if(curIteration % (REPORT_INTERVAL/2) == 0) {
	double factor = (double)pmec()/TARGET_MEC;
	// double factor = (double)pmec()/(TARGET_MEC == 0 ? 0.5 : TARGET_MEC);
	if(fracTime() - prev_retreat_frac > 2*SMALL_RETREAT &&
	    (((fracTime()>0.3||pBad<0.2) && factor > 16) ||   // 14 to 22 seems to work well
	     ((fracTime()>0.5||pBad<0.1) && factor >  8) )){  // quarter to half the above works well?
	    retreat = factor * SMALL_RETREAT / totalCoverage() * log(num_meta_iters);
	}
	if(fracTime()>FULL_RETREAT && factor > 1.3) {//pmec() != 0){
	    retreat = FULL_RETREAT;
	    ResetBuffers();
	}
	if(retreat > 0.0) {
	    cout << "Retreat " << 100*retreat << "% from " << 100 * fracTime();
	    this->curIteration -= retreat * this->maxIterations;
	    assert(this->curIteration>=0);
	    cout << "% to "  << 100 * fracTime() << "% because MEC is " << pmec();
	    cout << ", too big by a factor of " << factor << "(" << TARGET_MEC << 
			", " << totalWindowCoverage() << ")" << endl;
	    prev_retreat_frac = fracTime();
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
}

}
