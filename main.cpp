#include <iostream>
#include <fstream>
#include <random>
#include "Genome.hpp"

using namespace SAHap;
using namespace std;

int main(int argc, char *argv[]) {
	if (argc <= 1) {
		cerr << "Usage: " << argv[0] << " [reads]" << endl;
		return 0;
	}

	ifstream file;
	file.open(argv[1]);
	auto parsed = WIFInputReader::read(file, 2);

	try {
		Genome ge(parsed);
		cout << "Loaded" << endl;
		// cout << ge.chromosomes[0].percentAgree() << endl;
		for (float e = -5; e <= 2; e += 0.1) {
			float temp = pow(10, e);
			auto pbad = ge.findPbad(temp);
			cout << temp << " " << pbad << endl;
		}

	} catch (const char* e) {
		cout << e << endl;
	}
	// ge.optimize();

	// /*
	// for (float e = -10; e <= 10; ++e) {
	// 	float temp = pow(10, e);
	// 	cout << temp << ": " << ge.findPbad(temp) << endl;
	// }
	// */
	// /*
	// for (float i = 1; i > 0; i -= 0.1f) {
	// 	cout << "pBad @ " << i << " = " << ge.findPbad(100) << endl;
	// }
	// */

	return 0;
}

