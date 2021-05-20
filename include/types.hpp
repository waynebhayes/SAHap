#ifndef SAHAP_TYPES_HPP
#define SAHAP_TYPES_HPP

#define READ_ERROR_RATE 0.015 // 0.01 = 1% of letters on a read are incorrect due to sequencing errors

#include <limits.h>
#include "Allele.hpp"
#include <unordered_map>

namespace SAHap {

// Positions/Sites
typedef unsigned long dnapos_t;

// Read counts
typedef unsigned long dnacnt_t;

// Weight
typedef double dnaweight_t;

// Iterations
typedef unsigned long long iteration_t;

struct Site {
	Allele value = Allele::UNKNOWN;
	dnapos_t pos;
	double weight = 1;
};

enum class Zygosity {
	HOMO_REF, // '0'
	HOMO_ALT, // '1'
	HETERO,   // '*'
};

struct Range {
	dnapos_t start;
	dnapos_t end;

	Range(): start(INT_MAX), end(0) {}
	Range(dnapos_t a, dnapos_t b): start(a), end(b) {}
};

struct Read {
	vector<Site> sites;
	double cost = 0; // read-based cost (-log T_p(E_r, k_r))
	
	Range range;
};

struct InputFile {
	dnacnt_t ploidy;
	unordered_map<dnapos_t, dnapos_t> index; // index[matrix pos] = genome pos
	vector<dnapos_t> sites;
	vector<Read> reads;
	vector<Zygosity> zygosity;
	vector<vector<Allele>> groundTruth;
	dnacnt_t groundTruthNotCovered = 0;
	dnacnt_t averageReadLength = 0;
	bool hasZygosity = false; // For MEC/GI and WMEC/GI
	bool hasGroundTruth = false;
};

}

#endif
