#ifndef SAHAP_INPUTREADER_HPP
#define SAHAP_INPUTREADER_HPP

#include <algorithm>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <utility>
#include "types.hpp"

using namespace std;

namespace SAHap {

class WIFInputReader {
public:
	static InputFile read(ifstream& file, dnacnt_t ploidy);
	static void readGroundTruth(ifstream& file, InputFile& parsed);
	static Site parseSNP(string snp);
	// Map: actual pos -> matrix pos
	static Read parseRead(unordered_map<dnapos_t, dnapos_t>& index, string line);
};

}

#endif