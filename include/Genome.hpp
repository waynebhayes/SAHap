#ifndef GENOME_HPP
#define GENOME_HPP

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "Chromosome.hpp"
#include "InputReader.hpp"
#include "Read.hpp"
#include "types.hpp"

using namespace std;

namespace SAHap {

class Genome {
public:
	Genome(ifstream& filename, dnapos_t length, size_t ploidy);
	dnacnt_t get_mec();
	void sim_ann(float temp, float minTemp, float decreaseFactor, int numiters);

protected:
	vector<Chromosome> chrom_list;
	vector<ReadPair *> read_pair_list;

	float score(float temp, dnacnt_t mec);
	float chanceToKeepFunc(float curScore, float newScore, float temp);
};

}

#endif