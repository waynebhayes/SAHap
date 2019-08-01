#ifndef CHROMOSOME_HPP
#define CHROMOSOME_HPP

#include <algorithm>
#include <vector>
#include <unordered_set>
#include <random>
#include "DNAChar.hpp"
#include "Read.hpp"
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
	dnacnt_t mec();

	/**
	 * Add a read pair to a chromosome
	 */
	void add(ReadPair * rp);

	/**
	 * Remove a read pair to a chromosome
	 */
	void remove(ReadPair * rp);

	/**
	 * Randomly pick a ReadPair
	 */
	ReadPair * pick();
	ReadPair * pick(mt19937& engine);

	/**
	 * Print chromosome
	 */
	void print(ostream& stream, bool verbose=false);

	DNAChar * solution;

	friend ostream & operator << (ostream& stream, Chromosome& ch);

protected:
	dnapos_t length;
	dnacnt_t ** votes; // [site][letter] = # votes
	dnacnt_t * vsum; // sum of all votes
	dnacnt_t imec; // cached MEC
	bool mecDirty;
	unordered_set<ReadPair *> readPairs;

	void tally(dnapos_t site);
	void vote(const Read& read, bool retract=false);
};

ostream & operator << (ostream& stream, Chromosome& ch);

}

#endif