#include "Read.hpp"

using namespace std;

namespace SAHap {

ostream & operator << (ostream& stream, const Read& r) {
	stream << "+" << r.pos << " " << r.seq;
	return stream;
}

ostream & operator << (ostream& stream, const ReadPair& rp) {
	stream << "rp[" << rp.first << ", " << rp.second << "]";
	return stream;
}

}