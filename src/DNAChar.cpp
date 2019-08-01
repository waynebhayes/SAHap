#include "DNAChar.hpp"

using namespace std;

namespace SAHap {

ostream & operator << (ostream& stream, const DNAChar& l) {
	const char * letters = "ACGT.";
	stream << letters[l];
	return stream;
}

ostream & operator << (ostream& stream, const vector<DNAChar>& lv) {
	for (const auto& el : lv) {
		stream << el;
	}
	return stream;
}

}