#include "Chromosome.hpp"

namespace SAHap {

Chromosome::Chromosome(dnapos_t length)
	: length(length), imec(0), mecDirty(true)
{
	this->solution = vector<Allele>(this->length, Allele::UNKNOWN);
	this->votes = vector<VoteInfo>(this->length);

	for (dnapos_t i = 0; i < this->length; ++i) {
	}
}

Chromosome::Chromosome(const Chromosome& ch)
	:
		solution(ch.solution),
		length(ch.length),
		votes(ch.votes),
		imec(ch.imec),
		mecDirty(ch.mecDirty),
		reads(ch.reads)
{
}

Chromosome::~Chromosome() {
}

dnapos_t Chromosome::size() const {
	return this->length;
}

float Chromosome::mec() {
	if (this->mecDirty) {
		this->imec = 0;
		for (dnapos_t i = 0; i < this->length; ++i) {
			auto majority = this->solution[i];
			if (majority != Allele::UNKNOWN) {
				this->imec += this->votes[i].weight(majority);
			}
		}
		this->mecDirty = false;
	}

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
			this->imec -= this->votes[i].weight(flip_allele(this->solution[i]));
		}

		if (!retract) {
			this->votes[i].vote(allele)++;
			this->votes[i].weight(allele) += site.weight;

			if (
				allele != majority &&
				(
					majority == Allele::UNKNOWN ||
					this->votes[i].vote(allele) > this->votes[i].vote(flip_allele(allele))
				)
			) {
				majority = allele;
				this->solution[i] = allele;
			}
		} else {
			this->votes[i].vote(allele)--;
			this->votes[i].weight(allele) -= site.weight;

			if (allele == majority) {
				this->tally(i);
				majority = this->solution[i];
			}
		}

		if (majority != Allele::UNKNOWN) {
			this->imec += this->votes[i].weight(flip_allele(this->solution[i]));
		}
	}
}

void Chromosome::tally(dnapos_t site) {
	auto vote = this->votes[site];
	if (vote.ref_c == 0 && vote.alt_c == 0) {
		this->solution[site] = Allele::UNKNOWN;
	} else if (vote.ref_c >= vote.alt_c) {
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