#include "Allele.hpp"

using namespace std;

namespace SAHap {

size_t allele_i(Allele allele) {
	return (size_t)allele;
}

Allele flip_allele(Allele allele) {
	if (allele == Allele::UNKNOWN) {
		throw "Invalid allele value";
	}

	return allele == Allele::REF ? Allele::ALT : Allele::REF;
}

size_t flip_allele_i(Allele allele) {
	return allele_i(flip_allele(allele));
}

ostream & operator << (ostream& stream, const Allele& l) {
	if (l == Allele::UNKNOWN) stream << '-';
	else stream << (int)l;

	return stream;
}

}