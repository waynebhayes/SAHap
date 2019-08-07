#include "Genome.hpp"

using namespace std;

namespace SAHap {

Genome::Genome(ifstream& file, dnapos_t length, size_t ploidy)
	: chromosomes(vector<Chromosome>(ploidy, Chromosome(length)))
{
	random_device seed;
	this->randomEngine = mt19937(seed());
	this->readPairs = InputReader::read(file);
	this->makeInitialState();
}

Genome::~Genome() {
	for (auto& rp : this->readPairs) {
		delete rp;
	}
}

dnacnt_t Genome::mec() {
	dnacnt_t out = 0;
	for (size_t i = 0; i < chromosomes.size(); i++) {
		out = out + chromosomes[i].mec();
	}
	return out;
}

void Genome::makeInitialState() {
	uniform_int_distribution<size_t> distribution(0, this->chromosomes.size() - 1);

	for (auto& p : this->readPairs) {
		this->chromosomes[distribution(this->randomEngine)].add(p);
	}

	for (auto& ch : this->chromosomes) {
		cout << ch << endl;
	}
}


float Genome::score(float temp, dnacnt_t mec) {
	float max_mec = this->chromosomes.size() * this->chromosomes[0].size();
	float out = mec / max_mec; 
	return out;
}

float Genome::chanceToKeepFunc(float curScore, float newScore, float temp) {
	float p = pow(0.9, (1 - (curScore - newScore))) * temp;
	return p;
}

void Genome::optimize(float temp, float minTemp, float decreaseFactor, int numiters) {
	bool just_reverted = false;
	auto chromosomes_size = this->chromosomes.size();
	// main sim_annealing loop
	auto currentMEC = this->mec();

	while (temp > minTemp) {
		for (int i = 0; i < numiters; i++) {
			// calculate the MEC for the current chromosomes (unless we have just reverted, in which case we do not need to)
			if (!just_reverted) {
				currentMEC = this->mec();
			}
			// choose a pair read to remove and add to another chromosome
			int random_chrom_index = rand()%chromosomes_size;
			ReadPair * to_add = chromosomes[random_chrom_index].pick(this->randomEngine);
			chromosomes[random_chrom_index].remove(to_add);
			int chrom_difference = rand() % chromosomes_size;
			int random_chrom_index_2 = (random_chrom_index + chrom_difference) % chromosomes_size;
			chromosomes[random_chrom_index_2].add(to_add);
			// now, calculate the new MEC score
			auto newMEC = this->mec();
			// calculate the score of the old solution and the new solution
			auto scoreOfCurrent = score(temp, currentMEC);
			auto scoreOfNew = score(temp, newMEC);
			// use the scores to caculate the chance we keep this new solution, if it is bad
			float chanceToKeep = chanceToKeepFunc(scoreOfCurrent, scoreOfNew, temp);
			float randomIndex = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
			//if the solution is bad, AND we do not choose a float below the chance to keep, we revert to the old solution
			if (chanceToKeep <= randomIndex && newMEC > currentMEC) {
				chromosomes[random_chrom_index_2].remove(to_add);
				chromosomes[random_chrom_index].add(to_add);
			}

			cout << "[simann] temp=" << temp << ", mec=" << newMEC << endl;
		}
		temp = temp * decreaseFactor;
	}
}

}
