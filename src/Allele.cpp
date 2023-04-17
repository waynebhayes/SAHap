#include "Allele.hpp"

using namespace std;

namespace SAHap {

size_t allele_i(Allele allele) {
	return (size_t)allele;
}

ostream & operator << (ostream& stream, const Allele& l) {
	if (l == Allele::UNKNOWN) stream << 'X';
	else stream << (int)l;

	return stream;
}

}
