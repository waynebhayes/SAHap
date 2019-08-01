#ifndef READ_HPP
#define READ_HPP

#include <iostream>
#include <utility>
#include <vector>
#include "DNAChar.hpp"
#include "types.hpp"

using namespace std;

namespace SAHap {

struct Read {
	dnapos_t pos;
	vector<DNAChar> seq;
};

typedef pair<Read, Read> ReadPair;

ostream & operator << (ostream& stream, const Read& r);
ostream & operator << (ostream& stream, const ReadPair& rp);

}

#endif
