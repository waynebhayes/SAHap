#include "InputReader.hpp"

using namespace std;

namespace SAHap {

vector<ReadPair *> InputReader::read(ifstream& file) {
	vector<ReadPair *> result;
	string buf;

	while (!file.eof()) {
		getline(file, buf);
		if (buf.find("#") == 0) continue;
		result.push_back(InputReader::parseLine(buf));
	}

	return result;
}

ReadPair * InputReader::parseLine(string line) {
	istringstream iss(line);
	dnapos_t pos1, pos2;
	string seq1, seq2;
	Read r1, r2;

	iss >> pos1 >> seq1 >> pos2 >> seq2;
	r1.pos = pos1;
	r1.seq = InputReader::parseSeq(seq1);
	r2.pos = pos2;
	r2.seq = InputReader::parseSeq(seq2);

	auto result = new ReadPair(r1, r2);

	return result;
}

vector<DNAChar> InputReader::parseSeq(string seq) {
	vector<DNAChar> result;
	result.reserve(seq.length());

	for (const char& c : seq) {
		if (c == 'A' || c == 'a') result.push_back(DNAChar::A);
		else if (c == 'C' || c == 'c') result.push_back(DNAChar::C);
		else if (c == 'G' || c == 'g') result.push_back(DNAChar::G);
		else if (c == 'T' || c == 't') result.push_back(DNAChar::T);
		else throw "Invalid DNA letter";
	}

	return result;
}

}