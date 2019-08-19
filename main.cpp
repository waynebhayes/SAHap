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

	try {
		Genome ge(file, 2);
		cout << "Loaded" << endl;
		cout << ge.findPbad(pow(10, 10)) << endl;
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

