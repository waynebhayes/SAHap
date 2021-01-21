#include <iostream>
#include <fstream>
#include <random>
#include <cstdlib>
#include "Genome.hpp"

using namespace SAHap;
using namespace std;

int main(int argc, char *argv[]) {
	if (argc < 2 || argc > 4) {
		cerr << "Usage: " << argv[0] << " <reads> [gt] [millions of iterations = 10]" << endl;
		return 1;
	}

	ifstream file;
	file.open(argv[1]);
	auto parsed = WIFInputReader::read(file, 2);

	if (argc >= 3) {
		ifstream gtruth;
		gtruth.open(argv[2]);

		WIFInputReader::readGroundTruth(gtruth, parsed);
	}

	iteration_t iterations = argc == 4 ? atoi(argv[3]) * META_ITER : META_ITER;

	try {
		Genome ge(parsed);
		ge.start = 0;
		ge.end = ge.haplotypes[0].size();
		ge.autoSchedule(iterations);

		for (unsigned i = 500; i <= ge.haplotypes[0].size() + 100; i+=500) {
			ge.start = i-500;
			ge.end = i;

			try {
				ge.optimize(true, ge.start, ge.end);
				cerr << ge;
			} catch (const char * e) {
				cerr << e << endl;
			}
		}
		cerr << "Total Time: " << ge.total_time << endl;

		/*
		// cout << ge.haplotypes[0].percentAgree() << endl;
		for (float e = -5; e <= 2; e += 0.1) {
			float temp = pow(10, e);
			auto pbad = ge.findPbad(temp);
			cout << temp << " " << pbad << endl;
		}
		*/

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

