#include "Chromosome.hpp"

namespace SAHap {

Chromosome::Chromosome(dnapos_t length)
	: length(length), imec(0), mecDirty(true)
{
	this->solution = new DNAChar[this->length];
	this->votes = new dnacnt_t*[this->length];
	this->vsum = new dnacnt_t[this->length];

	for (dnapos_t i = 0; i < this->length; ++i) {
		this->solution[i] = DNAChar::UNKNOWN;
		this->votes[i] = new dnacnt_t[DNAChar::LENGTH];
		this->vsum[i] = 0;

		for (auto letter : AllDNAChar) {
			this->votes[i][letter] = 0;
		}
	}
}

Chromosome::Chromosome(const Chromosome& ch)
	: length(ch.length), imec(ch.imec), mecDirty(ch.mecDirty), readPairs(ch.readPairs)
{
	this->solution = new DNAChar[this->length];
	this->votes = new dnacnt_t*[this->length];
	this->vsum = new dnacnt_t[this->length];

	for (dnapos_t i = 0; i < this->length; ++i) {
		this->votes[i] = new dnacnt_t[DNAChar::LENGTH];
		copy(ch.votes[i], ch.votes[i] + DNAChar::LENGTH, this->votes[i]);
	}

	copy(ch.solution, ch.solution + ch.length, this->solution);
	copy(ch.vsum, ch.vsum + ch.length, this->vsum);
}

Chromosome::~Chromosome() {
	for (dnapos_t i = 0; i < this->length; ++i) {
		delete[] this->votes[i];
	}
	delete[] this->votes;
	delete[] this->vsum;
	delete[] this->solution;
}

dnapos_t Chromosome::size() const {
	return this->length;
}

dnacnt_t Chromosome::mec() {
	if (this->mecDirty) {
		this->imec = 0;
		for (dnapos_t i = 0; i < this->length; ++i) {
			auto majority = this->solution[i];
			if (majority != DNAChar::UNKNOWN) {
				this->imec += this->vsum[i] - this->votes[i][majority];
			}
		}
		this->mecDirty = false;
	}

	return this->imec;
}

void Chromosome::add(ReadPair * rp) {
	if (this->readPairs.find(rp) != this->readPairs.end()) {
		throw "Chromosome already contains read pair";
	}

	this->readPairs.insert(rp);
	this->vote(rp->first);
	this->vote(rp->second);
}

void Chromosome::remove(ReadPair * rp) {
	if (this->readPairs.find(rp) == this->readPairs.end()) {
		throw "Chromosome does not contain read pair";
	}

	this->readPairs.erase(rp);
	this->vote(rp->first, true);
	this->vote(rp->second, true);
}

ReadPair * Chromosome::pick() {
	random_device seed;
	mt19937 engine(seed());

	return this->pick(engine);
}

ReadPair * Chromosome::pick(mt19937& engine) {
	if (!this->readPairs.size()) return nullptr;

	uniform_int_distribution<size_t> distribution(0, this->readPairs.size() - 1);
	size_t rd = distribution(engine);
	return *next(begin(this->readPairs), rd);
}

void Chromosome::vote(const Read& read, bool retract) {
	dnapos_t rs = read.seq.size();

	for (dnapos_t i = 0; i < rs; ++i) {
		dnapos_t site = read.pos + i;
		DNAChar letter = read.seq[i];
		DNAChar majority = this->solution[site];

		if (majority != DNAChar::UNKNOWN) {
			this->imec -= this->vsum[site] - this->votes[site][majority];
		}

		if (!retract) {
			this->votes[site][letter]++;
			this->vsum[site]++;
			if (
				letter != majority &&
				(
					majority == DNAChar::UNKNOWN ||
					this->votes[site][letter] > this->votes[site][majority]
				)
			) {
				this->solution[site] = letter;
				majority = letter;
			}
		} else {
			this->votes[site][letter]--;
			this->vsum[site]--;

			if (letter == majority) {
				this->tally(site);
				majority = this->solution[site];
			}
		}

		if (majority != DNAChar::UNKNOWN) {
			this->imec += this->vsum[site] - this->votes[site][majority];
		}
	}
}

void Chromosome::tally(dnapos_t site) {
	auto majority = DNAChar::UNKNOWN;
	dnacnt_t max = 0;
	dnacnt_t sum = 0;

	for (auto l : AllDNAChar) {
		sum += this->votes[site][l];
		if (this->votes[site][l] > max) {
			max = this->votes[site][l];
			majority = l;
		}
	}

	this->vsum[site] = sum;
	this->solution[site] = majority;
}

ostream & operator << (ostream& stream, Chromosome& ch) {
	stream << "ch[";
	stream << "m=" << ch.length << ", ";
	stream << "n=" << ch.readPairs.size() << ", ";
	stream << "mec=" << ch.mec() << "] ";
	for (dnapos_t i = 0; i < ch.size(); ++i) {
		stream << ch.solution[i];
	}
	return stream;
}

}