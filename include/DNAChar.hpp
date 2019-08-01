#ifndef DNACHAR_HPP
#define DNACHAR_HPP

#include <iostream>
#include <vector>

using namespace std;

namespace SAHap {

enum DNAChar : unsigned char {
	BEGIN,
	A = BEGIN,
	C,
	G,
	T,
	END,
	LENGTH = END,
	UNKNOWN = END,
};

constexpr initializer_list<DNAChar> AllDNAChar = {A, C, G, T};

ostream & operator << (ostream& stream, const DNAChar& l);
ostream & operator << (ostream& stream, const vector<DNAChar>& lv);

}
#endif
