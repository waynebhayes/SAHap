#ifndef SAHAP_GENOME_HPP
#define SAHAP_GENOME_HPP

#include <iostream>
#include <fstream>
#include <random>
#include <sstream>
#include <vector>
#include <chrono>
#include <iomanip>
#include "Haplotype.hpp"
#include "InputReader.hpp"
#include "types.hpp"

#define META_ITER 10000 // how many iterations per integer on the command line? 1M? 100k?
#define REPORT_INTERVAL (META_ITER/10)

using namespace std;
using namespace std::chrono;

namespace SAHap {

class Genome {
public:
	Genome(InputFile file);
	~Genome();
	dnacnt_t mec();
	double mecScore();
	double siteCostScore();
	double score();
	double score(dnaweight_t mec);
	double totalCoverage();
	double fracTime();

	void shuffle();
	bool done();
	void setParameters(double tInitial, double tEnd, iteration_t maxIterations);
	void setTemperature(double t);
	void move();
	void revertMove();
	void iteration();
	void optimize(bool debug=false);
	void DynamicSchedule(double pBad, int TARGET_MEC);

	void Report(int seconds, bool final=false);
	double findPbad(double temperature, iteration_t iterations = REPORT_INTERVAL);
	void autoSchedule(iteration_t iterations);

	// fAccept = alpha = fraction of moves accepted "recently"
	struct AcceptBuffer {
		static constexpr int LENGTH = REPORT_INTERVAL;
		char buffer[LENGTH];
		size_t len = 0;
		size_t pos = 0;
		int sum = 0.0;

		void record(char good);
		double getAverage();
	};
	AcceptBuffer fAccept;
	int totalGood = 0;

	// pBad
	struct PbadBuffer {
		static constexpr int LENGTH = REPORT_INTERVAL;
		double buffer[LENGTH];
		size_t len = 0;
		size_t pos = 0;
		double sum = 0.0;

		void record(double acceptance);
		double getAverage();
	};
	PbadBuffer pBad;
	int totalBad = 0;
	int totalBadAccepted = 0;
	void ResetBuffers(); // Resets both buffers above


	vector<Haplotype> haplotypes;
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

	double acceptance(double newScore, double curScore);
	double getTemperature(iteration_t iteration);
	dnacnt_t compareGroundTruth(const Haplotype& ch, const vector<Allele>& truth);

	friend ostream& operator << (ostream& stream, const Genome& ge);
};

}

#endif
