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
	dnaweight_t WindowCost(dnapos_t start, dnapos_t end);
	double TotalCost();
	double Score();
	double Score(dnaweight_t mec);
	double MeanCoverage();
	double WindowTotalCoverage();
	double FracTime();

	void Shuffle();
	bool Done();
	void SetParameters(double tInitial, double tEnd, iteration_t maxIterations);
	void SetTemperature(double t);
	void Move();
	void RevertMove();
	void Iteration();
	void Optimize(bool debug);
	void DynamicSchedule(double pBad, double TARGET_MEC);

	void Report(int seconds, bool final=false);
	double FindPbad(double temperature, iteration_t iterations = REPORT_INTERVAL);
	void AutoSchedule(iteration_t iterations);

	// fAccept = alpha = fraction of moves accepted "recently"
	struct AcceptBuffer {
		static constexpr int LENGTH = REPORT_INTERVAL;
		char buffer[LENGTH];
		size_t len = 0;
		size_t pos = 0;
		int sum = 0.0;

		void Record(char good);
		double GetAverage();
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

		void Record(double acceptance);
		double GetAverage();
	};
	PbadBuffer pBad;
	int totalBad = 0;
	int totalBadAccepted = 0;
	void ResetBuffers(); // Resets both buffers above


	vector<Haplotype> haplotypes;
	void GenerateOutput();
	dnacnt_t CompareGroundTruth();

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

	Range range;

	dnapos_t numberOfSites = 0;
	dnacnt_t increments = 0;

	double lastErrorRate = 1;
	double lastCpuTime = 0;

	// The last move performed
	struct Move {
		size_t from;
		size_t to;
		Read * read;
	};
	struct Move lastMove;

	double Acceptance(double newScore, double curScore);
	double GetTemperature(iteration_t iteration);
	dnacnt_t CompareGroundTruth(const Haplotype& ch, const vector<int>& truth);
	
	friend ostream& operator << (ostream& stream, const Genome& ge);

private:
	vector<Range> blocks;

	void CreateBlocks();
	bool Intersects(Range a, Range b);
	Range CombineBlocks(Range a, Range b);
};

}

#endif
