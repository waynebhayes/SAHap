#include "Genome.hpp"
#include <csignal>
#include <cassert>

#define SAHAP_GENOME_DEBUG 0
#define POISSON 0 // set to zero to use MEC, 1 to use POISSON

using namespace std;

namespace SAHap {

Genome::Genome(InputFile file)
{
	random_device seed;
	this->randomEngine = mt19937(seed());
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

double Genome::siteCostScore() {
	double out = 0;

	for (size_t i = 0; i < haplotypes.size(); i++) {
		// -logT_p(lambda_i, k_i)
#if POISSON
		out = out + haplotypes[i].siteCost();
#else
		out = out + haplotypes[i].mec();
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
	// return this->pbad.getAverage() == 0 || this->t <= 1e-5;
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

	if (!(isGood || oldScore == newScore)) {
		this->totalBad++;
		if (accept) {
			this->totalBadAccepted++;
		}
		this->pbad.record(chanceToKeep);
	}

	/*
	uniform_real_distribution<double> d(0, 3);
	if (d(this->randomEngine) <= 1) {
		this->pbad.record(false);
	} else {
		this->pbad.record(true);
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
	// Reset state
	this->t = this->tInitial;
	this->pbad.total = 0;
	this->pbad.pos = 0;
	this->pbad.sum = 0;

	this->totalBad = 0;
	this->totalBadAccepted = 0;

	auto start_time = duration_cast<seconds>(system_clock::now().time_since_epoch());
	printf("Optimizing %d sites using %s cost function\n", (int)this->haplotypes[0].size(), POISSON ? "Poisson" : "MEC");
	cout << "Sites " << this->haplotypes[0].size() << endl;
	while (!this->done()) {
		this->t = this->getTemperature(this->curIteration);
		this->iteration();
		this->curIteration++;

		if (debug && curIteration % 10000 == 0) {
			auto now_time = duration_cast<seconds>(system_clock::now().time_since_epoch());
			printf("%dk (%.1f%%,%ds)  T %.3g  pBad %.3g  MEC %d", (int)this->curIteration/1000,
			    this->curIteration*100.0/this->maxIterations, (int)(now_time - start_time).count(),
			    this->t, this->pbad.getAverage(), (int)this->mec());
			if (this->file.hasGroundTruth) {
				auto gt = this->compareGroundTruth();
				double he = (double)gt / (this->haplotypes.size() * this->haplotypes[0].size());
				printf("  (Err-vs-truth %d = %g%%)", (int)gt, 100*he);
			}
			printf("\n");
		}
	}
	printf("Finished optimizing %d sites using %s cost function\n", (int)this->haplotypes[0].size(), POISSON ? "Poisson" : "MEC");
}

double Genome::findPbad(double temperature, iteration_t iterations) {
	this->shuffle();

	this->t = temperature;
	this->pbad.total = 0;
	this->pbad.pos = 0;
	this->pbad.sum = 0;

	this->totalBad = 0;
	this->totalBadAccepted = 0;

	for (iteration_t i = 0; i < iterations; ++i) {
		this->iteration();
		// cout << "[simann] temp=" << this->t << ", mec=" << this->mec() << endl;
	}

	cout << "temperature " << temperature << " gives Pbad " << this->pbad.getAverage() << endl;

	return this->pbad.getAverage();
}

void Genome::autoSchedule(iteration_t iterations) {
	cout << "Finding optimal temperature schedule..." << endl;

	double tInitial = 1;
	while (this->findPbad(tInitial, 10000) > .99) {
		tInitial /= 2;
		cout << "down" << endl;
	}
	while (this->findPbad(tInitial, 10000) < .99) {
		tInitial *= 1.2;
		//cout << "up" << endl;
	}

	cout << "tInitial found: " << tInitial << endl;
	iteration_t testiter = 20000;
	double tEnd = tInitial;
	while (true) {
		double pbad = this->findPbad(tEnd, 10000);
		testiter += 10000;

		if (pbad <= pow(10, -10)) {
			break;
		}

		tEnd /= 10;
	}

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
