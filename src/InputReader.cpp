#include "InputReader.hpp"

using namespace std;

namespace SAHap {

InputFile WIFInputReader::read(ifstream& file, dnacnt_t ploidy) {
	InputFile result;
	string buf;

	while (!file.eof()) {
		getline(file, buf);
		if (buf.size() == 0 || buf.find("#") == 0) continue;

		result.reads.push_back(WIFInputReader::parseRead(result.index, buf));
	}

	result.ploidy = ploidy;
	return result;
}

void WIFInputReader::readGroundTruth(ifstream& file, InputFile& parsed) {
	// Obtain sorted list of sites covered
	vector<dnapos_t> sites;
	sites.reserve(parsed.index.size());
	for (const auto& kv : parsed.index) {
		sites.push_back(kv.first);
	}
	sort(sites.begin(), sites.end());
	parsed.sites = sites;

	string buf;
	vector<vector<Allele>> truth;
	truth.reserve(parsed.ploidy);

	while (!file.eof()) {
		getline(file, buf);
		if (buf.size() == 0) continue;
		vector<Allele> ch(parsed.index.size(), Allele::UNKNOWN);
		for (size_t i = 0; i < buf.size(); ++i) {
			char allele = buf[i];
			dnapos_t pos = sites[i];
			dnapos_t matrixPos = parsed.index[pos];
			if (allele == '0') ch[matrixPos] = Allele::REF;
			else if (allele == '1') ch[matrixPos] = Allele::ALT;
			else if (allele == 'X') ch[matrixPos] = Allele::UNKNOWN;
			else throw "Invalid ground truth allele value";
		}
		truth.push_back(ch);
	}

	if (truth.size() != parsed.ploidy) {
		cerr << "Warning: Number of chromosomes in ground truth does not match ploidy! Results will be meaningless." << endl;
	}

	parsed.groundTruth = truth;
	parsed.hasGroundTruth = true;
}

Site WIFInputReader::parseSNP(string snp) {
	istringstream iss(snp);
	Site s;
	iss >> s.pos;

	iss.ignore(256, ' ');
	iss.ignore(256, ' ');

	int weight, value;
	iss >> value >> weight;

	if (s.weight > 0 && s.weight <= 100) {
		s.weight = (float)weight / 100;
	} else {
		throw "Invalid weight value";
	}
	s.weight = 1; // HACK

	if (value == 0) {
		s.value = Allele::REF;
	} else if (value == 1) {
		s.value = Allele::ALT;
	} else {
		cout << value << endl;
		throw "Invalid allele value";
	}
	return s;
}

Read WIFInputReader::parseRead(unordered_map<dnapos_t, dnapos_t>& index, string line) {
	Read result;
	istringstream iss(line);
	string buf;
	while (getline(iss, buf, ':')) {
		istringstream siss(buf);
		char a;
		siss >> a;
		if (a == '#') break;

		Site snp = WIFInputReader::parseSNP(buf);
		if (index.find(snp.pos) == index.end()) {
			// not a known site
			dnapos_t matrixPos = index.size();
			index.insert(make_pair(snp.pos, matrixPos));
			snp.pos = matrixPos;
		} else {
			// known site
			// index[actual pos] = matrix pos
			snp.pos = index[snp.pos];
		}

		result.sites.push_back(snp);
	}

	return result;
}

}