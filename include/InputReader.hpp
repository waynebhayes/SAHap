#ifndef INPUTREADER_HPP
#define INPUTREADER_HPP

#include <fstream>
#include <sstream>
#include <vector>
#include <utility>
#include "DNAChar.hpp"
#include "Read.hpp"
#include "types.hpp"

using namespace std;

namespace SAHap {

class InputReader {
public:
	static vector<ReadPair *> read(ifstream& file);
	static ReadPair * parseLine(string line);
	static vector<DNAChar> parseSeq(string seq);
};

}

#endif