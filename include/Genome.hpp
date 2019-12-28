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
	double score();
	double score(dnaweight_t mec);

	void shuffle();
	bool done();
	void setParameters(double tInitial, double tEnd, iteration_t maxIterations);
	void setTemperature(double t);
	void move();
	void revertMove();
	void iteration();
	void optimize(bool debug=false);
	double findPbad(double temperature);

	// pBad
	struct PbadBuffer {
		static constexpr int LENGTH = 10000;
		bool buffer[LENGTH];
		size_t total = 0;
		size_t pos = 0;
		size_t sum = 0;

		void record(bool bad);
		double getAverage();
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

	double t = 100000;
	double tInitial = 100000;
	double tDecay = 0.000001;
	iteration_t maxIterations = 0;
	iteration_t curIteration = 0;

	// The last move performed
	struct Move {
		size_t from;
		size_t to;
		Read * read;
	};
	Move lastMove;

	double acceptance(dnaweight_t newMec, dnaweight_t curMec);
	double getTemperature(iteration_t iteration);
	dnacnt_t compareGroundTruth(const Chromosome& ch, const vector<Allele>& truth);

	friend ostream& operator << (ostream& stream, const Genome& ge);
};

}

#endif
