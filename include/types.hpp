#ifndef SAHAP_TYPES_HPP
#define SAHAP_TYPES_HPP

#include "Allele.hpp"
#include <unordered_map>

namespace SAHap {

// Positions/Sites
typedef unsigned long dnapos_t;

// Read counts
typedef unsigned long dnacnt_t;

// Weight
typedef float dnaweight_t;

// Iterations
typedef unsigned long long iteration_t;

struct Site {
	Allele value = Allele::UNKNOWN;
	dnapos_t pos;
	float weight = 1;
};

enum class Zygosity {
	HOMO_REF, // '0'
	HOMO_ALT, // '1'
	HETERO,   // '*'
};

struct Read {
	vector<Site> sites;
};

struct InputFile {
	dnacnt_t ploidy;
	unordered_map<dnapos_t, dnapos_t> index; // index[matrix pos] = genome pos
	vector<dnapos_t> sites;
	vector<Read> reads;
	vector<Zygosity> zygosity;
	vector<vector<Allele>> groundTruth;
	dnacnt_t groundTruthNotCovered = 0;
	bool hasZygosity = false; // For MEC/GI and WMEC/GI
	bool hasGroundTruth = false;
};

}

#endif