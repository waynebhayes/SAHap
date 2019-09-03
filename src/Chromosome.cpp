#include "Chromosome.hpp"

#define SAHAP_CHROMOSOME_DEBUG_MEC 0
#if SAHAP_CHROMOSOME_DEBUG_MEC
#include <csignal>
#endif

namespace SAHap {

Chromosome::Chromosome(dnapos_t length)
	: length(length), imec(0)
{
	this->solution = vector<Allele>(this->length, Allele::UNKNOWN);
	this->weights = vector<array<dnacnt_t, 2>>(this->length);

	for (dnapos_t i = 0; i < this->length; ++i) {
		this->solution[i] = Allele::UNKNOWN;
		this->weights[i][0] = 0;
		this->weights[i][1] = 0;
	}
}

Chromosome::Chromosome(const Chromosome& ch)
	:
		solution(ch.solution),
		length(ch.length),
		weights(ch.weights),
		imec(ch.imec),
		reads(ch.reads)
{
}

Chromosome::~Chromosome() {
}

dnapos_t Chromosome::size() const {
	return this->length;
}

size_t Chromosome::readSize() const {
	return this->reads.size();
}

float Chromosome::mec() {
#if SAHAP_CHROMOSOME_DEBUG_MEC
	float imec = 0;
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

void Chromosome::add(Read * r) {
	if (this->reads.find(r) != this->reads.end()) {
		throw "Chromosome already contains read";
	}

	this->reads.insert(r);
	this->vote(*r);
}

void Chromosome::remove(Read * r) {
	if (this->reads.find(r) == this->reads.end()) {
		cout << "Offending read is " << r << endl;
		// cout << "Offending read has #sites=" << r->sites.size();
		throw "Chromosome does not contain read";
	}

	this->reads.erase(r);
	this->vote(*r, true);
}

Read * Chromosome::pick() {
	random_device seed;
	mt19937 engine(seed());

	return this->pick(engine);
}

Read * Chromosome::pick(mt19937& engine) {
	if (!this->reads.size()) return nullptr;

	uniform_int_distribution<size_t> distribution(0, this->reads.size() - 1);
	size_t rd = distribution(engine);
	return *next(begin(this->reads), rd);
}

void Chromosome::vote(const Read& read, bool retract) {
	for (const Site& site : read.sites) {
		dnapos_t i = site.pos;
		Allele allele = site.value;
		Allele majority = this->solution[i];

		if (majority != Allele::UNKNOWN) {
			this->imec -= this->weights[i][flip_allele_i(this->solution[i])];
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
		} else {
			// leave
			this->weights[i][allele_i(allele)] -= site.weight;

			if (allele == majority) {
				this->tally(i);
				majority = this->solution[i];
			}
		}

		if (majority != Allele::UNKNOWN) {
			this->imec += this->weights[i][flip_allele_i(this->solution[i])];
		}
	}
}

void Chromosome::tally(dnapos_t site) {
	auto weights = this->weights[site];
	if (weights[allele_i(Allele::REF)] == 0 && weights[allele_i(Allele::ALT)] == 0) {
		this->solution[site] = Allele::UNKNOWN;
	} else if (weights[allele_i(Allele::REF)] >= weights[allele_i(Allele::ALT)]) {
		this->solution[site] = Allele::REF;
	} else {
		this->solution[site] = Allele::ALT;
	}
}

dnacnt_t& Chromosome::VoteInfo::vote(Allele allele) {
	if (allele == Allele::REF) {
		return this->ref_c;
	} else if (allele == Allele::ALT) {
		return this->alt_c;
	}
	throw "vote: Invalid allele value";
}

float& Chromosome::VoteInfo::weight(Allele allele) {
	if (allele == Allele::REF) {
		return this->ref_w;
	} else if (allele == Allele::ALT) {
		return this->alt_w;
	}
	throw "weight: Invalid allele value";
}
	
ostream & operator << (ostream& stream, Chromosome& ch) {
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
