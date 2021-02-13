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
	Haplotype(dnapos_t length);
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
	 * Compute the MEC
	 */
	double mec();

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

	vector<Allele> solution;

	friend ostream & operator << (ostream& stream, Haplotype& ch);

	double mec(dnapos_t s, dnapos_t e);
	double ave_coverage();
	void save_reads();
	void sep_reads();
	void print_mec();
	bool check(Read *r);
	Range range;
	double pmec = 0;
protected:
	struct VoteInfo {
		dnacnt_t ref_c = 0;
		dnacnt_t alt_c = 0;
		double ref_w = 0;
		double alt_w = 0;

		dnacnt_t& vote(Allele allele);
		double& weight(Allele allele);
	};

	dnapos_t length;
	// vector<VoteInfo> votes;
	vector<array<dnacnt_t, 2>> weights;
	vector<dnacnt_t> siteCoverages;

	double imec = 0; // cached MEC
	double isitecost = 0; // cached site-based cost
	unordered_set<Read *> reads;
	unordered_set<Read *> sreads;

	unordered_set<Read *> archived;

	void tally(dnapos_t site);
	void vote(Read& read, bool retract=false);
};

ostream & operator << (ostream& stream, Haplotype& ch);

}

#endif
