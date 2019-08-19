#ifndef DNACHAR_HPP
#define DNACHAR_HPP

#include <iostream>
#include <vector>

using namespace std;

namespace SAHap {

enum class DNAChar : unsigned char {
	BEGIN,
	A = BEGIN,
	C,
	G,
	T,
	END,
	LENGTH = END,
	UNKNOWN = END,
};

constexpr initializer_list<DNAChar> AllDNAChar = {DNAChar::A, DNAChar::C, DNAChar::G, DNAChar::T};

ostream & operator << (ostream& stream, const DNAChar& l);
ostream & operator << (ostream& stream, const vector<DNAChar>& lv);

}
#endif
