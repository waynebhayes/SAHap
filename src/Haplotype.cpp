#include "Haplotype.hpp"

#define SAHAP_CHROMOSOME_DEBUG_MEC 0 // Not sure what this is meant to debug.... chromosome?
#if SAHAP_CHROMOSOME_DEBUG_MEC
#include <csignal>
#endif
#define SAHAP_CHROMOSOME_ALT_MEC 0
#define SMALL_ENOUGH_TO_IGNORE 1e-10

namespace SAHap {

Haplotype::Haplotype(dnapos_t numSites, unsigned ploidyCount)
	: numSites(numSites), total_mec(0), window_mec(0), isitecost(0)
{
	this->ploidyCount = ploidyCount;
	this->solution = vector<int>(this->numSites, -1);
	this->weights = vector<vector<int>>(this->numSites);
	this->siteCoverages = vector<dnacnt_t>(this->numSites);
	this->window.start = 0;
	this->window.end = this->numSites;

	for (dnapos_t i = 0; i < this->numSites; ++i) {
		this->solution[i] = -1;
		this->weights[i] = vector<int>(this->ploidyCount, 0);
		this->siteCoverages[i] = 0;
	}
}

Haplotype::Haplotype(const Haplotype& ch)
	:
		solution(ch.solution),
		numSites(ch.numSites),
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

double Haplotype::MeanCoverage() {
    double result = 0.0;
    for (dnapos_t i = 0; i < this->numSites; ++i) {
		result += this->siteCoverages[i];
    }
    return result/this->numSites;
}

dnapos_t Haplotype::NumSites() const {
	return this->numSites;
}

size_t Haplotype::NumReads() const {
	return this->reads.size();
}

// Compute the MEC across a window [s,e] for this haplotype ("side")
dnaweight_t Haplotype::WindowCost(dnapos_t s, dnapos_t e) {
    dnaweight_t out = 0;
    assert(e>=s);

    for (unsigned j = 0; j < ploidyCount; j++) {
	for (dnapos_t i = s; i <= e && i < this->numSites; i++) {
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

dnaweight_t Haplotype::TotalCost() {
    return WindowCost(0, this->numSites);
}

void Haplotype::SaveReads() {
	for (auto r : this->reads) {
		if (r->range.end > this->window.start)
			this->saved_reads.insert(r);
	}
}

void Haplotype::PickReads(unsigned overlap) {
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

void Haplotype::InitializeWindow(unsigned windowSize, unsigned incrementBy) {
	this->window.start = 0;
	this->window.end = windowSize;
	this->increment_window_by = incrementBy;
	this->saved_reads = this->reads;

	this->PickReads(0);

	this->window_mec = WindowCost(this->window.start, this->window.end);
}

void Haplotype::IncrementWindow() {
	dnapos_t old_end = this->window.end;

	this->window.start += this->increment_window_by;
	this->window.end = min(this->window.end + this->increment_window_by, this->numSites);

	this->SaveReads();
	this->PickReads(old_end - this->window.start);
	
	this->window_mec = WindowCost(this->window.start, this->window.end);
}

double Haplotype::WindowTotalCoverage() {
    double result = 0.0;
    for (dnapos_t i = this->window.start; i < this->window.end; ++i) {
		result += this->siteCoverages[i];
    }
    return result;//(this->window.end - this->window.start);
}

void Haplotype::PrintCoverages() {
	for (auto siteCoverage : siteCoverages) {
		cerr << siteCoverage << " ";
	}
	cerr << endl;
}

// FIXME: most functions returing "double" should probably return dnaweight_t instead.
dnaweight_t Haplotype::SiteCost(const Site &s) {
	cerr << "siteCost needs to be implemented\n";
	return -1;
}

void Haplotype::AddRead(Read * r) {
	if (this->reads.find(r) != this->reads.end()) {
		throw "Haplotype already contains read";
	}
	// std::cout << "adding\n";
	this->reads.insert(r);
	this->Vote(*r);
}

void Haplotype::RemoveRead(Read * r) {
	if (this->reads.find(r) == this->reads.end()) {
		cout << "Offending read is " << r << endl;
		// cout << "Offending read has #sites=" << r->sites.size();
		throw "Haplotype does not contain read";
	}

	this->reads.erase(r);
	this->Vote(*r, true);
}

Read * Haplotype::RandomRead() {
	long int seed = GetFancySeed(true);
	cout << "Haplotype seed " << seed << endl;
	mt19937 engine(seed);

	return this->RandomRead(engine);
}

Read * Haplotype::RandomRead(mt19937& engine) {
	if (this->reads.size()==0) return nullptr;

	uniform_int_distribution<size_t> distribution(0, this->reads.size() - 1);
	size_t rd = distribution(engine);
	return *next(begin(this->reads), rd);
}

bool Haplotype::IsInRangeOf(Range r, dnapos_t pos) {
	return pos >= r.start && pos <= r.end;
}

void Haplotype::SubtractMECValuesAt(dnapos_t pos) {
	for (unsigned i = 0; i < ploidyCount; i++) {
		if ((int)i == solution[pos])
			continue;
		auto mec = weights[pos][i];
		total_mec -= mec;
		assert(total_mec>=0);

		if (IsInRangeOf(window, pos))
			window_mec -= mec;
		if(window_mec < 0 && window_mec > -SMALL_ENOUGH_TO_IGNORE) window_mec = 0;
		assert(window_mec>=0);

		if(siteCoverages[i]) isitecost -= -log_poisson_1_cdf(READ_ERROR_RATE * siteCoverages[i], mec);
	}
}

void Haplotype::AddMECValuesAt(dnapos_t pos) {
	for (unsigned i = 0; i < ploidyCount; i++) {
		if ((int)i == solution[pos])
			continue;
		auto mec = weights[pos][i];
		total_mec += mec;

		if (IsInRangeOf(window, pos))
			window_mec += mec;
		if(window_mec < 0 && window_mec > -SMALL_ENOUGH_TO_IGNORE) window_mec = 0;
		assert(window_mec>=0);

		if(siteCoverages[i]) isitecost += -log_poisson_1_cdf(READ_ERROR_RATE * siteCoverages[i], mec);
	}
}

void Haplotype::AddSite(const Site &s) {
	weights[s.pos][s.value] += s.weight;

	if (s.value != solution[s.pos] && weights[s.pos][s.value] > weights[s.pos][solution[s.pos]])
		solution[s.pos] = s.value;
	
	siteCoverages[s.pos]++;
}

void Haplotype::RemoveSite(const Site &s) {
	weights[s.pos][s.value] -= s.weight;

	if (solution[s.pos] == s.value)
		FindSolution(s.pos);

	siteCoverages[s.pos]--;
}

void Haplotype::Vote(Read& read, bool retract) {
#if SAHAP_CHROMOSOME_ALT_MEC
	// TODO: Alternative MEC
#else
	for (const Site& site : read.sites) {
		dnapos_t i = site.pos;

		SubtractMECValuesAt(i);

		if (!retract) // enter
			AddSite(site);
		else // leave
			RemoveSite(site);

		AddMECValuesAt(i);
	}
#endif
}

void Haplotype::FindSolution(dnapos_t site) {
	for (unsigned i = 0; i < ploidyCount; i++) 
		if (weights[site][i] > weights[site][solution[site]])
			solution[site] = i;
}

dnacnt_t& Haplotype::VoteInfo::Vote(Allele allele) {
	if (allele == Allele::REF) {
		return this->ref_c;
	} else if (allele == Allele::ALT) {
		return this->alt_c;
	}
	throw "vote: Invalid allele value";
}

int& Haplotype::VoteInfo::Weight(Allele allele) {
	if (allele == Allele::REF) {
		return this->ref_w;
	} else if (allele == Allele::ALT) {
		return this->alt_w;
	}
	throw "weight: Invalid allele value";
}
	
ostream & operator << (ostream& stream, Haplotype& ch) {
	stream << "ch[";
	stream << "m=" << ch.numSites << ", ";
	stream << "n=" << ch.reads.size() << ", ";
	stream << "mec=" << ch.TotalCost() << "] ";
	for (dnapos_t i = 0; i < ch.NumSites(); ++i) {
		stream << ch.solution[i];
	}
	return stream;
}

}