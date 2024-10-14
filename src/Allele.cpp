#include "Allele.hpp"

using namespace std;

namespace SAHap {

//Expects: allele of type Allele
//Returns: size_t of allele 
size_t allele_i(Allele allele) {
	return (size_t)allele;
}
//Expects: Allele l 
//Return: if known, prints the integer at allele l, otherwise prints 'X'  
//
//Question: what is l?
ostream & operator << (ostream& stream, const Allele& l) {
	if (l == Allele::UNKNOWN) stream << 'X';
	else stream << (int)l;

	return stream;
}

}
