#ifndef SAHAP_CHROMOSOME_HPP
#define SAHAP_CHROMOSOME_HPP

#include <algorithm>
#include <vector>
#include <unordered_set>
#include <random>
#include "DNAChar.hpp"
#include "types.hpp"

namespace SAHap {

class Chromosome {
public:
	Chromosome(dnapos_t length);
	Chromosome(const Chromosome& ch);
	~Chromosome();

	/**
	 * Returns the size of the solution
	 */
	dnapos_t size() const;

	/**
	 * Compute the MEC
	 */
	float mec();

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

	friend ostream & operator << (ostream& stream, Chromosome& ch);

protected:
	struct VoteInfo {
		dnacnt_t ref_c = 0;
		dnacnt_t alt_c = 0;
		float ref_w = 0;
		float alt_w = 0;

		dnacnt_t& vote(Allele allele);
		float& weight(Allele allele);
	};

	dnapos_t length;
	vector<VoteInfo> votes;

	float imec; // cached MEC
	bool mecDirty;
	unordered_set<Read *> reads;

	void tally(dnapos_t site);
	void vote(const Read& read, bool retract=false);
};

ostream & operator << (ostream& stream, Chromosome& ch);

}

#endif