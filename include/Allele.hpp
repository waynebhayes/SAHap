#ifndef ALLELE_HPP
#define ALLELE_HPP

#include <iostream>

using namespace std;

namespace SAHap {

enum class Allele : unsigned char { REF, ALT, UNKNOWN };

size_t allele_i(Allele allele);
Allele flip_allele(Allele allele);
size_t flip_allele_i(Allele allele);
ostream & operator << (ostream& stream, const Allele& l);

}

#endif