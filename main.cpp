#include <iostream>
#include <fstream>
#include <random>
#include "Chromosome.hpp"
#include "InputReader.hpp"

//#include <utility>
//#include "Read.hpp"
//#include "DNAChar.hpp"

using namespace SAHap;
using namespace std;

int main(int argc, char *argv[]) {
	if (argc <= 1) {
		cerr << "Usage: " << argv[0] << " [reads]" << endl;
		return 0;
	}

	ifstream file;
	file.open(argv[1]);
	auto readPairs = InputReader::read(file);

	// FIXME: Specify m, k
	const int k = 2, m = 1000;
	vector<Chromosome*> ch;
	ch.reserve(k);
	for (int i = 0; i < k; ++i) {
		ch.push_back(new Chromosome(m));
	}

	cout << "Setting up initial state..." << endl;
	random_device seed;
	mt19937 engine(seed());
	uniform_int_distribution<size_t> distribution(0, k - 1);

	for (auto& p : readPairs) {
		ch[distribution(engine)]->add(p);
	}

	for (int i = 0; i < k; ++i) {
		cout << *ch[i] << endl;
	}

	cout << "Cleaning up..." << endl;
	for (int i = 0; i < k; ++i) {
		delete ch[i];
	}

	for (auto& p : readPairs) {
		delete p;
	}

	return 0;

	/*
	// More tests
	Read r1;
	r1.pos = 12;
	r1.seq = vector<DNAChar> {DNAChar::A, DNAChar::G};

	Read r2;
	r2.pos = 36;
	r2.seq = vector<DNAChar> {DNAChar::T, DNAChar::C};

	pair<Read, Read> rp1(r1, r2);

	Read q1;
	q1.pos = 12;
	q1.seq = vector<DNAChar> {DNAChar::A, DNAChar::T};

	pair<Read, Read> rp2(q1, r2);

	cout << r1 << rp1 << endl;

	Chromosome ch(50);
	ch.add(&rp1);
	ch.add(&rp2);

	cout << ch << endl;

	cout << *ch.pick() << endl;

	ch.remove(&rp1);
	cout << ch << endl;

	ch.remove(&rp2);
	cout << ch << endl;
	*/
}

