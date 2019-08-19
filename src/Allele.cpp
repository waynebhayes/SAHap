#include "Allele.hpp"

using namespace std;

namespace SAHap {

Allele flip_allele(Allele allele) {
	if (allele == Allele::UNKNOWN) {
		throw "Invalid allele value";
	}

	return allele == Allele::REF ? Allele::ALT : Allele::REF;
}

ostream & operator << (ostream& stream, const Allele& l) {
	if (l == Allele::UNKNOWN) stream << '-';
	else stream << (int)l;

	return stream;
}

}