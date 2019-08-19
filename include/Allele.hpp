#ifndef ALLELE_HPP
#define ALLELE_HPP

#include <iostream>

using namespace std;

namespace SAHap {

enum class Allele : char { REF, ALT, UNKNOWN };

Allele flip_allele(Allele allele);
ostream & operator << (ostream& stream, const Allele& l);

}

#endif