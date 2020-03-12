#include "Haplotype.hpp"

#define SAHAP_CHROMOSOME_DEBUG_MEC 0
#if SAHAP_CHROMOSOME_DEBUG_MEC
#include <csignal>
#endif
#define SAHAP_CHROMOSOME_ALT_MEC 0

namespace SAHap {

Haplotype::Haplotype(dnapos_t length)
	: length(length), imec(0), isitecost(0)
{
	this->solution = vector<Allele>(this->length, Allele::UNKNOWN);
	this->weights = vector<array<dnacnt_t, 2>>(this->length);
	this->siteCoverages = vector<dnacnt_t>(this->length);

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
		imec(ch.imec),
		reads(ch.reads)
{
}

Haplotype::~Haplotype() {
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
	if (imec != this->imec) {
		cerr << "DEBUG: Bad MEC value " << this->imec << ", should be " << imec << endl;
		raise(SIGINT);
	}
#endif

	return this->imec;
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
	random_device seed;
	mt19937 engine(seed());

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
			this->imec -= mec;
			// FIXME
			this->isitecost -= -log_poisson_1_cdf(0.15 * this->siteCoverages[i], mec);
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
			this->imec += mec;
			// FIXME
			this->isitecost += -log_poisson_1_cdf(0.15 * this->siteCoverages[i], mec);
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
