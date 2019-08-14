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

	Genome ge(file, 1000, 2);
	// ge.optimize();

	return 0;
}

