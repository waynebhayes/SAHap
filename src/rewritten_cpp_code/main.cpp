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
    auto parsed = WIFInputReader::read(file);

    if (argc > 3) {
        ifstream gtruth;
        gtruth.open(argv[2]);

        WIFInputReader::readGroundTruth(gtruth, parsed);
    }




    iteration_t iterations = argc == 4 ?  atoi(argv[3]) * META_ITER : 10 * META_ITER;

    try {
        Genome ge(parsed);
        ge.AutoSchedule(iterations);
        try {
            ge.Optimize(true);
            cout << ge;
        } catch (const char * e) {
            cerr << e << endl;
        }

    } catch (const char* e) {
        cout << e << endl;
    }

    return 0;
}

