#include "Genome.hpp"

using namespace std;

namespace SAHap {

Genome::Genome(ifstream& file, dnapos_t length, size_t ploidy)
	: chrom_list(vector<Chromosome>(ploidy, Chromosome(length)))
{
	//use inputReader.hpp to fill the read_pair_list
	read_pair_list = InputReader::read(file);
}

dnacnt_t Genome::get_mec() {
	dnacnt_t out = 0;
	for(size_t i = 0; i < chrom_list.size(); i++){
		out = out + chrom_list[i].mec();
	}
	return out;
}

float Genome::score(float temp, dnacnt_t mec) {
	// FIXME: Implement
	return (float)mec;
}

float Genome::chanceToKeepFunc(float curScore, float newScore, float temp) {
	// FIXME: Implement
	return newScore < curScore ? 1.0 : 0.0;
}

void Genome::sim_ann(float temp, float minTemp, float decreaseFactor, int numiters) {
	bool just_reverted = false;
	auto chrom_list_size = chrom_list.size();
	// iterate over list of read pairs, and add each one to one random chromosome in chrom_list
	for (size_t i = 0; i < read_pair_list.size(); i++) {
		int random_chrom_index = rand() % chrom_list_size;
		chrom_list[random_chrom_index].add(read_pair_list[i]);
	}
	// main sim_annealing loop
	auto currentMEC = this->get_mec();

	while (temp > minTemp) {
		for (int i = 0; i < numiters; i++) {
			// calculate the MEC for the current chromosomes (unless we have just reverted, in which case we do not need to)
			if (!just_reverted) {
				currentMEC = this->get_mec();
			}
			// choose a pair read to remove and add to another chromosome
			int random_chrom_index = rand()%chrom_list_size;
			ReadPair * to_add = chrom_list[random_chrom_index].pick();
			chrom_list[random_chrom_index].remove(to_add);
			int chrom_difference = rand() % chrom_list_size;
			int random_chrom_index_2 = (random_chrom_index + chrom_difference) % chrom_list_size;
			chrom_list[random_chrom_index_2].add(to_add);
			// now, calculate the new MEC score
			auto newMEC = this->get_mec();
			// calculate the score of the old solution and the new solution
			auto scoreOfCurrent = score(temp, currentMEC);
			auto scoreOfNew = score(temp, newMEC);
			// use the scores to caculate the chance we keep this new solution, if it is bad
			float chanceToKeep = chanceToKeepFunc(scoreOfCurrent, scoreOfNew, temp);
			float randomIndex = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
			//if the solution is bad, AND we do not choose a float below the chance to keep, we revert to the old solution
			if (chanceToKeep <= randomIndex && newMEC > currentMEC) {
				chrom_list[random_chrom_index_2].remove(to_add);
				chrom_list[random_chrom_index].add(to_add);
			}
		}
		temp = temp * decreaseFactor;
	}
}

}
