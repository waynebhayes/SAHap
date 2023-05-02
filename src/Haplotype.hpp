#ifndef SAHAP_HAPLOTYPE_HPP
#define SAHAP_HAPLOTYPE_HPP

#include <algorithm>
#include <array>
#include <vector>
#include <unordered_set>
#include <random>
#include "types.hpp"
#include "utils.hpp"

namespace SAHap {

class Haplotype {
public:
	Haplotype(dnapos_t length, unsigned ploidyCount);
	Haplotype(const Haplotype& ch);
	~Haplotype();

	double meanCoverage();

	/**
	 * Returns the size of the solution
	 */
	dnapos_t size() const;

	/**
	 * Returns the number of reads in the Haplotype
	 */
	size_t readSize() const;

	/**
	 * Compute the total MEC
	 */
	double mec();

	/**
	 * Compute partial MEC 
	 */ 
	double mec(dnapos_t start, dnapos_t end);

	/**
	 * Compute current window's MEC
	 */
	double windowMec();

	/**
	 * Compute average coverage of SNPs within the Window
	 */
	double windowMeanCoverage();

	/**
	 * Compute the site-based cost
	 */
	double siteCost();

	/**
	 * Add a Read to a chromosome
	 */
	void add(Read * r);

	/**
	 * Remove a Read to a chromosome
	 */
	void remove(Read * r);

	/**
	 * Randomly pick a Read
	 */
	Read * pick();
	Read * pick(mt19937& engine);

	/**
	 * Print chromosome
	 */
	void print(ostream& stream, bool verbose=false);

	/**
	 * Initializes range of the window and how much it advances each turn
	 */ 
	void initializeWindow(unsigned windowSize, unsigned incrementBy);

	/**
	 * Increments the current window that the program is trying to solve by incrementBy
	 */
	void incrementWindow();

	vector<int> solution;

	friend ostream & operator << (ostream& stream, Haplotype& ch);

	// void print_mec(); // Only for debugging
	void printCoverages();
protected:
	struct VoteInfo {
		dnacnt_t ref_c = 0;
		dnacnt_t alt_c = 0;
		int ref_w = 0;
		int alt_w = 0;

		dnacnt_t& vote(Allele allele);
		int& weight(Allele allele);
	};

	dnapos_t length;
	// vector<VoteInfo> votes;
	vector<vector<int>> weights;
	vector<dnacnt_t> siteCoverages;

	double total_mec = 0; // cached MEC
	double window_mec = 0; // cached current window's MEC
	double isitecost = 0; // cached site-based cost

	unsigned ploidyCount;

	unordered_set<Read *> reads;
	unordered_set<Read *> saved_reads;

	Range window;
	unsigned increment_window_by;

	void findSolution(dnapos_t site);
	void vote(Read& read, bool retract=false);
	void saveReads();
	void pickReads(unsigned overlap);

	bool isInRangeOf(Range r, dnapos_t pos);

	void subtractMECValuesAt(dnapos_t pos);
	void addMECValuesAt(dnapos_t pos);
	void addSite(const Site &s);
	void removeSite(const Site &s);
};

ostream & operator << (ostream& stream, Haplotype& ch);

}

#endif
