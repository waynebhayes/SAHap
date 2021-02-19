#include "Haplotype.hpp"

#define SAHAP_CHROMOSOME_DEBUG_MEC 0
#if SAHAP_CHROMOSOME_DEBUG_MEC
#include <csignal>
#endif
#define SAHAP_CHROMOSOME_ALT_MEC 0

namespace SAHap {

Haplotype::Haplotype(dnapos_t length)
	: length(length), total_mec(0), window_mec(0), isitecost(0)
{
	this->solution = vector<Allele>(this->length, Allele::UNKNOWN);
	this->weights = vector<array<dnacnt_t, 2>>(this->length);
	this->siteCoverages = vector<dnacnt_t>(this->length);
	this->window.start = 0;
	this->window.end = this->length;

	for (dnapos_t i = 0; i < this->length; ++i) {
		this->solution[i] = Allele::UNKNOWN;
		this->weights[i][0] = 0;
		this->weights[i][1] = 0;
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

size_t Haplotype::readSize() const {
	return this->reads.size();
}

double Haplotype::mec() {
#if SAHAP_CHROMOSOME_DEBUG_MEC
	double imec = 0;
	for (dnapos_t i = 0; i < this->length; ++i) {
		this->tally(i);
		auto majority = this->solution[i];
		if (majority != Allele::UNKNOWN) {
			imec += this->weights[i][flip_allele_i(majority)];
		}
	}
	if (imec != this->total_mec) {
		cerr << "DEBUG: Bad MEC value " << this->total_mec << ", should be " << imec << endl;
		raise(SIGINT);
	}
#endif

	return this->total_mec;
}

double Haplotype::mec(dnapos_t s, dnapos_t e) {
	double out = 0;

	for (auto i = s; i <= e && i < this->length; i++) {
		if (this->solution[i] == Allele::UNKNOWN)
			out += this->weights[i][0] + this->weights[i][1];
		else
			out += this->weights[i][flip_allele_i(this->solution[i])];
	}

	return out;
}

double Haplotype::windowMec() {
	return this->window_mec;
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
		int scope = (int)min(r->range.end, this->window.end) - max(r->range.start, this->window.start + overlap);
		if (scope > 0) {// Maybe we should just pick reads that overlaps with the current window by more than 1 SNPs
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

	this->window.start += increment_window_by;
	this->window.end += increment_window_by;

	this->saveReads();
	this->pickReads(old_end - this->window.start);
	
	this->window_mec = mec(this->window.start, this->window.end);
}

void Haplotype::print_mec() {
	for (dnapos_t i = 0; i < this->length; i++) {
		if (this->solution[i] == Allele::UNKNOWN)
			cerr << this->weights[i][0] + this->weights[i][1];
		else 
			cerr << this->weights[i][flip_allele_i(this->solution[i])];
	}
	cerr << endl;
}

double Haplotype::siteCost() {
	return this->isitecost;
}

void Haplotype::add(Read * r) {
	if (this->reads.find(r) != this->reads.end()) {
		throw "Haplotype already contains read";
	}

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

Read * Haplotype::pick() {
	long int seed = GetFancySeed(true);
	cout << "Haplotype seed " << seed << endl;
	mt19937 engine(seed);

	return this->pick(engine);
}

Read * Haplotype::pick(mt19937& engine) {
	if (!this->reads.size()) return nullptr;

	uniform_int_distribution<size_t> distribution(0, this->reads.size() - 1);
	size_t rd = distribution(engine);
	return *next(begin(this->reads), rd);
}

void Haplotype::vote(Read& read, bool retract) {
#if SAHAP_CHROMOSOME_ALT_MEC
	// TODO: Alternative MEC
#else
	for (const Site& site : read.sites) {
		dnapos_t i = site.pos;
		Allele allele = site.value;
		Allele majority = this->solution[i];

		if (majority != Allele::UNKNOWN) {
			auto mec = this->weights[i][flip_allele_i(this->solution[i])];
			this->total_mec -= mec;
			if (i >= this->window.start && i <= this->window.end)
				this->window_mec -= mec;
			// FIXME: we have only 1 bit to specify the "main" letter or *THREE* altertanes, so
			// they are not equally probable. Need to account for this lopsidedness.
			this->isitecost -= -log_poisson_1_cdf(READ_ERROR_RATE * this->siteCoverages[i], mec);
		}

		if (!retract) {
			// enter
			this->weights[i][allele_i(allele)] += site.weight;

			if (
				allele != majority &&
				(
					majority == Allele::UNKNOWN ||
					this->weights[i][allele_i(allele)] > this->weights[i][flip_allele_i(allele)]
				)
			) {
				majority = allele;
				this->solution[i] = allele;
			}

			this->siteCoverages[i]++;
		} else {
			// leave
			this->weights[i][allele_i(allele)] -= site.weight;

			if (allele == majority) {
				this->tally(i);
				majority = this->solution[i];
			}

			this->siteCoverages[i]--;
		}

		if (majority != Allele::UNKNOWN) {
			auto mec = this->weights[i][flip_allele_i(this->solution[i])];
			this->total_mec += mec;
			if (i >= this->window.start && i <= this->window.end)
				this->window_mec += mec;
			// FIXME: we have only 1 bit to specify the "main" letter or *THREE* altertanes, so
			// they are not equally probable. Need to account for this lopsidedness.
			this->isitecost += -log_poisson_1_cdf(READ_ERROR_RATE * this->siteCoverages[i], mec);
		}
	}
#endif
}

void Haplotype::tally(dnapos_t site) {
	auto weights = this->weights[site];
	if (weights[allele_i(Allele::REF)] == 0 && weights[allele_i(Allele::ALT)] == 0) {
		this->solution[site] = Allele::UNKNOWN;
	} else if (weights[allele_i(Allele::REF)] >= weights[allele_i(Allele::ALT)]) {
		this->solution[site] = Allele::REF;
	} else {
		this->solution[site] = Allele::ALT;
	}
}

dnacnt_t& Haplotype::VoteInfo::vote(Allele allele) {
	if (allele == Allele::REF) {
		return this->ref_c;
	} else if (allele == Allele::ALT) {
		return this->alt_c;
	}
	throw "vote: Invalid allele value";
}

double& Haplotype::VoteInfo::weight(Allele allele) {
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
