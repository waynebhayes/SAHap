#include "Genome.hpp"
#include <csignal>

#define SAHAP_GENOME_DEBUG 0

using namespace std;

namespace SAHap {

Genome::Genome(InputFile file)
{
	random_device seed;
	this->randomEngine = mt19937(seed());
	this->file = file;
	this->chromosomes = vector<Chromosome>(this->file.ploidy, Chromosome(this->file.index.size()));
	this->shuffle();
}

Genome::~Genome() {
}

dnaweight_t Genome::mec() {
	dnaweight_t out = 0;
	for (size_t i = 0; i < chromosomes.size(); i++) {
		out = out + chromosomes[i].mec();
	}
	return out;
}

float Genome::score(dnaweight_t mec) {
	float maxMec = this->chromosomes.size() * this->chromosomes[0].size();
	return mec / maxMec;
}

float Genome::score() {
	return this->score(this->mec());
}

void Genome::shuffle() {
	uniform_int_distribution<size_t> distribution(0, this->chromosomes.size() - 1);

	this->maxIterations = this->chromosomes[0].size() * 100;

	if (this->initialized) {
		auto ploidy = this->chromosomes.size();
		auto length = this->chromosomes[0].size();
		this->chromosomes.clear();
		this->chromosomes = vector<Chromosome>(ploidy, Chromosome(length));
	}

	for (auto& r : this->file.reads) {
		this->chromosomes[distribution(this->randomEngine)].add(&r);
	}

	this->initialized = true;

	// for (auto& ch : this->chromosomes) {
	// 	cout << ch << endl;
	// }
}

float Genome::acceptance(dnaweight_t newMec, dnaweight_t curMec) {
	if (newMec < curMec) return 1;
	if (this->t == 0) return 0;
	float energyDiff = this->score(curMec) - this->score(newMec);
	// cout << "Acceptance(" << energyDiff << ") = " << exp(energyDiff / this->t) << endl;
	return exp(energyDiff / this->t);
}

bool Genome::done() {
	// Are we tired yet?
	return this->curIteration >= this->maxIterations;
}

void Genome::setParameters(float tInitial, float tDecay, iteration_t maxIterations) {
	if (this->curIteration) {
		// Probably warn the user in some way that the simulation is in process
	}

	this->tInitial = tInitial;
	this->tDecay = tDecay;
	this->maxIterations = maxIterations;
}

void Genome::setTemperature(float t) {
	this->t = t;
}

void Genome::move() {
	// Perform a random move, saving enough information so we can revert later

	auto ploidy = this->chromosomes.size();
	size_t moveFrom = rand() % ploidy;
	size_t moveTo;

	while (!this->chromosomes[moveFrom].readSize()) {
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

	Read * r = this->chromosomes[moveFrom].pick(this->randomEngine);
#if SAHAP_GENOME_DEBUG
	if (!r) {
		cerr << "DEBUG: Chromosome " << moveFrom << " has no reads remaining" << endl;
		std::raise(SIGINT);
		r = this->chromosomes[moveFrom].pick(this->randomEngine);
	}
#endif
	this->chromosomes[moveTo].add(r);
	this->chromosomes[moveFrom].remove(r);

	this->lastMove.from = moveFrom;
	this->lastMove.to = moveTo;
	this->lastMove.read = r;
}

void Genome::revertMove() {
	const auto& move = this->lastMove;
	this->chromosomes[move.to].remove(move.read);
	this->chromosomes[move.from].add(move.read);
}

void Genome::iteration() {
	// Run an iteration
	auto oldMec = this->mec();
	this->move();
	auto newMec = this->mec();

	uniform_real_distribution<float> distribution(0, 1);
	float chanceToKeep = this->acceptance(newMec, oldMec);
	float randomIndex = distribution(this->randomEngine);

	bool isBad = newMec > oldMec;
	bool accept = randomIndex <= chanceToKeep;

	if (isBad) {
		this->totalBad++;
		if (accept) {
			this->totalBadAccepted++;
		}
		this->pbad.record(accept);
	}

	if (accept) {
		// accept move
	} else {
		// reject move
		this->revertMove();
	}
	/*
	uniform_real_distribution<float> d(0, 3);
	if (d(this->randomEngine) <= 1) {
		this->pbad.record(false);
	} else {
		this->pbad.record(true);
	}
	*/
}

float Genome::getTemperature(iteration_t iteration) {
	float temp = this->tInitial * pow(this->tDecay, iteration);
	cout << iteration << ": " << temp << endl;
	return temp;
}

void Genome::optimize() {
	while (!this->done()) {
		this->t = this->getTemperature(this->curIteration);
		this->iteration();
		this->curIteration++;
		cout << "[simann] temp=" << this->t << ", mec=" << this->mec() << endl;
	}
}

float Genome::findPbad(float temperature) {
	this->shuffle();

	this->t = temperature;
	this->pbad.total = 0;
	this->pbad.pos = 0;
	this->pbad.sum = 0;

	this->totalBad = 0;
	this->totalBadAccepted = 0;

	for (size_t i = 0; i < 10000; ++i) {
		this->iteration();
		// cout << "[simann] temp=" << this->t << ", mec=" << this->mec() << endl;
	}

	return this->pbad.getAverage();
}

void Genome::PbadBuffer::record(bool bad) {
	size_t next = this->pos == LENGTH - 1 ? 0 : this->pos + 1;
	if (this->total == LENGTH) {
		this->sum -= this->buffer[next]; // the next one is the first one
	} else {
		this->total++;
	}
	this->pos = next;
	this->buffer[next] = bad;
	this->sum += bad;
}

float Genome::PbadBuffer::getAverage() {
	return (float)this->sum / this->total;
}

dnacnt_t Genome::compareGroundTruth() {
	auto l00 = this->compareGroundTruth(chromosomes[0], file.groundTruth[0]);
	auto l11 = this->compareGroundTruth(chromosomes[1], file.groundTruth[1]);
	auto la = l00 + l11;
	// cout << "00: " << l00 << ", 11: " << l11 << ", total: " << la << endl;

	auto l01 = this->compareGroundTruth(chromosomes[0], file.groundTruth[1]);
	auto l10 = this->compareGroundTruth(chromosomes[1], file.groundTruth[0]);
	auto lb = l01 + l10;
	// cout << "01: " << l01 << ", 10: " << l10 << ", total: " << lb << endl;

	return la < lb ? la : lb;
}

dnacnt_t Genome::compareGroundTruth(const Chromosome& ch, const vector<Allele>& truth) {
	dnacnt_t loss = 0;
	for (size_t i = 0; i < this->chromosomes[0].size(); ++i) {
		if (ch.solution[i] != truth[i]) {
			loss++;
		}
	}
	return loss;
}

/*
void Genome::optimize(float temp, float minTemp, float decreaseFactor, int numiters) {
	bool just_reverted = false;
	auto ploidy = this->chromosomes.size();
	// main sim_annealing loop
	auto currentMEC = this->mec();

	while (temp > minTemp) {
		for (int i = 0; i < numiters; i++) {
			// calculate the MEC for the current chromosomes (unless we have just reverted, in which case we do not need to)
			if (!just_reverted) {
				currentMEC = this->mec();
			}
			// choose a pair read to remove and add to another chromosome
			int moveFrom = rand()%ploidy;
			ReadPair * rp = chromosomes[moveFrom].pick(this->randomEngine);
			chromosomes[moveFrom].remove(rp);
			int moveOffset = rand() % ploidy;
			int moveTo = (moveFrom + moveOffset) % ploidy;
			chromosomes[moveTo].add(rp);
			// now, calculate the new MEC score
			auto newMEC = this->mec();
			// calculate the score of the old solution and the new solution
			auto scoreOfCurrent = score(temp, currentMEC);
			auto scoreOfNew = score(temp, newMEC);
			// use the scores to caculate the chance we keep this new solution, if it is bad
			float chanceToKeep = acceptance(scoreOfCurrent, scoreOfNew, temp);
			float randomIndex = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
			//if the solution is bad, AND we do not choose a float below the chance to keep, we revert to the old solution
			if (chanceToKeep <= randomIndex && newMEC > currentMEC) {
				chromosomes[moveTo].remove(rp);
				chromosomes[moveFrom].add(rp);
			}

			cout << "[simann] temp=" << temp << ", mec=" << newMEC << endl;
		}
		temp = temp * decreaseFactor;
	}
}
*/

}
