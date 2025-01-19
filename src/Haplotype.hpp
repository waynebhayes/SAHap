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

	double MeanCoverage();

	dnapos_t NumSites() const;

	/**
	 * Returns the number of reads in the Haplotype
	 */
	size_t NumReads() const;

	/**
	 * Compute the total Cost across all sites.
	 */
	dnaweight_t TotalCost();

	/**
	 * Compute cost across a window
	 */ 
	dnaweight_t WindowCost(dnapos_t start, dnapos_t end);

	dnaweight_t SiteCost(const Site &s);
	/**
	 * Compute average coverage of SNPs within the Window
	 */
	double WindowTotalCoverage();

	/**
	 * Add a Read to this haplotype
	 */
	void AddRead(Read * r);

	/**
	 * Remove a Read to this haplotype
	 */
	void RemoveRead(Read * r);

	/**
	 * Randomly pick a Read
	 */
	Read * RandomRead();
	Read * RandomRead(mt19937& engine);

	/**
	 * Print chromosome
	 */
	void Print(ostream& stream, bool verbose=false);

	/**
	 * Initializes range of the window and how much it advances each turn
	 */ 
	void InitializeWindow(unsigned windowSize, unsigned incrementBy);

	/**
	 * Increments the current window that the program is trying to solve by incrementBy .... what?
	 */
	void IncrementWindow();

	vector<int> solution; // FIXME: is this a list of reads and which side they're on, or a list of sites with expected letter?

	friend ostream & operator << (ostream& stream, Haplotype& ch);

	// void print_mec(); // Only for debugging
	void PrintCoverages();
protected:
	struct VoteInfo {
		dnacnt_t ref_c = 0;
		dnacnt_t alt_c = 0;
		int ref_w = 0;
		int alt_w = 0;

		dnacnt_t& Vote(Allele allele);
		int& Weight(Allele allele);
	};

	// Why is this needed? Isn't the number of sites a global, the same for all Haplotypes and also for the whole Genome?
	dnapos_t numSites; // number of sites
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

	void FindSolution(dnapos_t site);
	void Vote(Read& read, bool retract=false);
	void SaveReads();
	void PickReads(unsigned overlap);

	bool IsInRangeOf(Range r, dnapos_t pos);

	void SubtractMECValuesAt(dnapos_t pos);
	void AddMECValuesAt(dnapos_t pos);
	void AddSite(const Site &s);
	void RemoveSite(const Site &s);
};

ostream & operator << (ostream& stream, Haplotype& ch);

}

#endif