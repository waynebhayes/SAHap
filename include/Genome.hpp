#ifndef SAHAP_GENOME_HPP
#define SAHAP_GENOME_HPP

#include <iostream>
#include <fstream>
#include <random>
#include <sstream>
#include <vector>
#include "Chromosome.hpp"
#include "InputReader.hpp"
#include "types.hpp"

using namespace std;

namespace SAHap {

class Genome {
public:
	Genome(InputFile file);
	~Genome();
	dnaweight_t mec();
	float score();
	float score(dnaweight_t mec);

	void shuffle();
	bool done();
	void setParameters(float tInitial, float tDecay, iteration_t maxIterations);
	void setTemperature(float t);
	void move();
	void revertMove();
	void iteration();
	void optimize();
	float findPbad(float temperature);

	// pBad
	struct PbadBuffer {
		static constexpr int LENGTH = 10000;
		bool buffer[LENGTH];
		size_t total = 0;
		size_t pos = 0;
		size_t sum = 0;

		void record(bool bad);
		float getAverage();
	};
	PbadBuffer pbad;
	int totalBad = 0;
	int totalBadAccepted = 0;
	vector<Chromosome> chromosomes;
	void generateOutput();
	dnacnt_t compareGroundTruth();

protected:
	InputFile file;
	mt19937 randomEngine;
	bool initialized = false;

	bool lastMoves[1000];
	size_t lastMovesFront = 0;
	size_t lastMovesBack = 0;

	float t = 100000;
	float tInitial = 100000;
	float tDecay = 0.98;
	iteration_t maxIterations = 0;
	iteration_t curIteration = 0;

	// The last move performed
	struct Move {
		size_t from;
		size_t to;
		Read * read;
	};
	Move lastMove;

	float acceptance(dnaweight_t newMec, dnaweight_t curMec);
	float getTemperature(iteration_t iteration);
	dnacnt_t compareGroundTruth(const Chromosome& ch, const vector<Allele>& truth);
};

}

#endif