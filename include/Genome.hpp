#ifndef GENOME_HPP
#define GENOME_HPP

#include <iostream>
#include <fstream>
#include <random>
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
	~Genome();
	dnacnt_t mec();
	void optimize(float temp, float minTemp, float decreaseFactor, int numiters);

protected:
	vector<Chromosome> chromosomes;
	vector<ReadPair *> readPairs;
	mt19937 randomEngine;

	void makeInitialState();
	float score(float temp, dnacnt_t mec);
	float chanceToKeepFunc(float curScore, float newScore, float temp);
};

}

#endif