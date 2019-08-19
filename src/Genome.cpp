#include "Genome.hpp"

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

void Genome::move() {
	// Perform a random move, saving enough information so we can revert later

	auto ploidy = this->chromosomes.size();
	size_t moveFrom = rand() % ploidy;
	size_t moveOffset = rand() % ploidy;
	size_t moveTo = (moveFrom + moveOffset) % ploidy;

	Read * r = this->chromosomes[moveFrom].pick(this->randomEngine);
	this->chromosomes[moveFrom].remove(r);
	this->chromosomes[moveTo].add(r);

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

	if (chanceToKeep <= randomIndex) {
		this->revertMove();
	} else {
		this->pbad.record(newMec > oldMec);
		/*
		uniform_real_distribution<float> d(0, 3);
		if (d(this->randomEngine) <= 1) {
			this->pbad.record(false);
		} else {
			this->pbad.record(true);
		}
		*/
	}
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

	for (size_t i = 0; i < 100000; ++i) {
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
