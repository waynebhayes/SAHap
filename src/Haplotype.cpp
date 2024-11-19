#include "Haplotype.hpp"

#define SAHAP_CHROMOSOME_DEBUG_MEC 0 // Not sure what this is meant to debug.... chromosome?
#if SAHAP_CHROMOSOME_DEBUG_MEC
#include <csignal>
#endif
#define SAHAP_CHROMOSOME_ALT_MEC 0
#define SMALL_ENOUGH_TO_IGNORE 1e-10

namespace SAHap {

Haplotype::Haplotype(dnapos_t length, unsigned ploidyCount)
	: length(length), total_mec(0), window_mec(0), isitecost(0)
{
	this->ploidyCount = ploidyCount;
	this->solution = vector<int>(this->length, -1);
	this->weights = vector<vector<int>>(this->length);
	this->siteCoverages = vector<dnacnt_t>(this->length);
	this->window.start = 0;
	this->window.end = this->length;

	for (dnapos_t i = 0; i < this->length; ++i) {
		this->solution[i] = -1;
		this->weights[i] = vector<int>(this->ploidyCount, 0);
		this->siteCoverages[i] = 0;
	}
}

Haplotype::Haplotype(const Haplotype& ch)
	:
		solution(ch.solution),
		length(ch.length),
		weights(ch.weights),
		siteCoverages(ch.siteCoverages),
		total_mec(ch.total_mec),
		window_mec(ch.window_mec),
		ploidyCount(ch.ploidyCount),
		reads(ch.reads),
		saved_reads(ch.saved_reads),
		window(ch.window),
		increment_window_by(ch.increment_window_by)
{
}

Haplotype::~Haplotype() {
}

double Haplotype::meanCoverage() {
    double result = 0.0;
    for (dnapos_t i = 0; i < this->length; ++i) {
		result += this->siteCoverages[i];
    }
    return result/this->length;
}

dnapos_t Haplotype::size() const {
	return this->length;
}

size_t Haplotype::numReads() const {
	return this->reads.size();
}

// Compute the MEC across a window [s,e] for this haplotype ("side")
dnaweight_t Haplotype::windowMEC(dnapos_t s, dnapos_t e) {
    dnaweight_t out = 0;
    assert(e>=s);

    for (unsigned j = 0; j < ploidyCount; j++) {
	for (dnapos_t i = s; i <= e && i < this->length; i++) {
	    if (solution[i] == (int)j) // solution is signed since (-1) is used to mean "undefined"
		continue;
	    if(weights[i][j] < 0 && weights[i][j] > -SMALL_ENOUGH_TO_IGNORE) weights[i][j] = 0;
	    assert(weights[i][j]>=0);
	    out += weights[i][j];
	}
    }
    assert(out>=0);
    return out;
}

dnaweight_t Haplotype::windowCost() {
    cerr << "H:windowCost should choose between Poisson and MEC and call the appropriate function";
    return -1;
}

void Haplotype::saveReads() {
	for (auto r : this->reads) {
		if (r->range.end > this->window.start)
			this->saved_reads.insert(r);
	}
}

void Haplotype::pickReads(unsigned overlap) {
	this->reads.clear();
	
	for (auto r : this->saved_reads) {
		if (r->range.end > this->window.start + overlap && r->range.start <= this->window.end) {
			this->reads.insert(r);
		}
	}
	//cout << "MIN START: " << mx << endl;
	for (auto r : this->reads)
		this->saved_reads.erase(r);
}

void Haplotype::initializeWindow(unsigned windowSize, unsigned incrementBy) {
	this->window.start = 0;
	this->window.end = windowSize;
	this->increment_window_by = incrementBy;
	this->saved_reads = this->reads;

	this->pickReads(0);

	this->window_mec = mec(this->window.start, this->window.end);
}

void Haplotype::incrementWindow() {
	dnapos_t old_end = this->window.end;

	this->window.start += this->increment_window_by;
	this->window.end = min(this->window.end + this->increment_window_by, this->length);

	this->saveReads();
	this->pickReads(old_end - this->window.start);
	
	this->window_mec = mec(this->window.start, this->window.end);
}

double Haplotype::windowTotalCoverage() {
    double result = 0.0;
    for (dnapos_t i = this->window.start; i < this->window.end; ++i) {
		result += this->siteCoverages[i];
    }
    return result;//(this->window.end - this->window.start);
}

void Haplotype::printCoverages() {
	for (auto siteCoverage : siteCoverages) {
		cerr << siteCoverage << " ";
	}
	cerr << endl;
}

// FIXME: most functions returing "double" should probably return dnaweight_t instead.
dnaweight_t Haplotype::siteCost(const Site &s) {
	cerr << "siteCost needs to be implemented\n";
}

void Haplotype::add(Read * r) {
	if (this->reads.find(r) != this->reads.end()) {
		throw "Haplotype already contains read";
	}
	// std::cout << "adding\n";
	this->reads.insert(r);
	this->vote(*r);
}

void Haplotype::remove(Read * r) {
	if (this->reads.find(r) == this->reads.end()) {
		cout << "Offending read is " << r << endl;
		// cout << "Offending read has #sites=" << r->sites.size();
		throw "Haplotype does not contain read";
	}

	this->reads.erase(r);
	this->vote(*r, true);
}

Read * Haplotype::randomRead() {
	long int seed = GetFancySeed(true);
	cout << "Haplotype seed " << seed << endl;
	mt19937 engine(seed);

	return this->randomRead(engine);
}

Read * Haplotype::randomRead(mt19937& engine) {
	if (this->reads.size()==0) return nullptr;

	uniform_int_distribution<size_t> distribution(0, this->reads.size() - 1);
	size_t rd = distribution(engine);
	return *next(begin(this->reads), rd);
}

bool Haplotype::isInRangeOf(Range r, dnapos_t pos) {
	return pos >= r.start && pos <= r.end;
}

void Haplotype::subtractMECValuesAt(dnapos_t pos) {
	for (unsigned i = 0; i < ploidyCount; i++) {
		if ((int)i == solution[pos])
			continue;
		auto mec = weights[pos][i];
		total_mec -= mec;
		assert(total_mec>=0);

		if (isInRangeOf(window, pos))
			window_mec -= mec;
		if(window_mec < 0 && window_mec > -SMALL_ENOUGH_TO_IGNORE) window_mec = 0;
		assert(window_mec>=0);

		if(siteCoverages[i]) isitecost -= -log_poisson_1_cdf(READ_ERROR_RATE * siteCoverages[i], mec);
	}
}

void Haplotype::addMECValuesAt(dnapos_t pos) {
	for (unsigned i = 0; i < ploidyCount; i++) {
		if ((int)i == solution[pos])
			continue;
		auto mec = weights[pos][i];
		total_mec += mec;

		if (isInRangeOf(window, pos))
			window_mec += mec;
		if(window_mec < 0 && window_mec > -SMALL_ENOUGH_TO_IGNORE) window_mec = 0;
		assert(window_mec>=0);

		if(siteCoverages[i]) isitecost += -log_poisson_1_cdf(READ_ERROR_RATE * siteCoverages[i], mec);
	}
}

void Haplotype::addSite(const Site &s) {
	weights[s.pos][s.value] += s.weight;

	if (s.value != solution[s.pos] && weights[s.pos][s.value] > weights[s.pos][solution[s.pos]])
		solution[s.pos] = s.value;
	
	siteCoverages[s.pos]++;
}

void Haplotype::removeSite(const Site &s) {
	weights[s.pos][s.value] -= s.weight;

	if (solution[s.pos] == s.value)
		findSolution(s.pos);

	siteCoverages[s.pos]--;
}

void Haplotype::vote(Read& read, bool retract) {
#if SAHAP_CHROMOSOME_ALT_MEC
	// TODO: Alternative MEC
#else
	for (const Site& site : read.sites) {
		dnapos_t i = site.pos;

		subtractMECValuesAt(i);

		if (!retract) // enter
			addSite(site);
		else // leave
			removeSite(site);

		addMECValuesAt(i);
	}
#endif
}

void Haplotype::findSolution(dnapos_t site) {
	for (unsigned i = 0; i < ploidyCount; i++) 
		if (weights[site][i] > weights[site][solution[site]])
			solution[site] = i;
}

dnacnt_t& Haplotype::VoteInfo::vote(Allele allele) {
	if (allele == Allele::REF) {
		return this->ref_c;
	} else if (allele == Allele::ALT) {
		return this->alt_c;
	}
	throw "vote: Invalid allele value";
}

int& Haplotype::VoteInfo::weight(Allele allele) {
	if (allele == Allele::REF) {
		return this->ref_w;
	} else if (allele == Allele::ALT) {
		return this->alt_w;
	}
	throw "weight: Invalid allele value";
}
	
ostream & operator << (ostream& stream, Haplotype& ch) {
	stream << "ch[";
	stream << "m=" << ch.length << ", ";
	stream << "n=" << ch.reads.size() << ", ";
	stream << "mec=" << ch.mec() << "] ";
	for (dnapos_t i = 0; i < ch.size(); ++i) {
		stream << ch.solution[i];
	}
	return stream;
}

}
