#include "Genome.hpp"
#include "utils.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <iostream>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <string>

#define SAHAP_GENOME_DEBUG 0

enum _objectives        { OBJ_NONE,   OBJ_MEC,   OBJ_Poisson };
const char *objName[] = {"OBJ_NONE",     "MEC",     "Poisson"};
#ifndef OBJECTIVE
#define OBJECTIVE OBJ_MEC // choices for now are MEC and Poisson
#endif
#if (OBJECTIVE != OBJ_MEC && OBJECTIVE != OBJ_Poisson)
#error "invalid objective"
#endif

#define TARGET_PBAD_START 0.85
#define TARGET_PBAD_END 1e-3
// How is the temperature schedule adjusted to be dynamic?
enum _schedule            { SCHED_NONE,   SCHED_RETREAT,   SCHED_Betz };
const char *schedName[] = {"SCHED_NONE",       "RETREAT",       "Betz"};
#ifndef SCHEDULE
#define SCHEDULE SCHED_RETREAT // choices for now are SCHED_RETREAT and SCHED_Betz
#endif
#if (SCHEDULE != RETREAT && SCHEDULE != Betz)
#error "invalid schedule"
#endif

int maxCoverageAssumption = 120;


using namespace std;

namespace SAHap {


Genome::Genome(InputFile file)
{
	long int seed = GetFancySeed(true);
	cout << "Genome seed " << seed << endl;
	this->randomEngine = mt19937(seed);
	this->file = file;
	this->haplotypes = vector<Haplotype>(this->file.ploidy, Haplotype(this->file.index.size(), file.ploidy));
	this->range.start = 0;
	this->range.end = this->file.index.size();
	this->numberOfSites = this->file.index.size();
	this->increments = file.averageReadLength;
	this->Shuffle();
}

Genome::~Genome() {
}


dnaweight_t Genome::TotalCost() {
	dnaweight_t out = 0;
	auto ploidy = this->haplotypes.size();
	for (size_t i = 0; i < ploidy; i++) {
		out = out + haplotypes[i].TotalCost();
	}
	double maxMec = ploidy * this->haplotypes[0].NumSites() * (maxCoverageAssumption/ploidy);
	return out; // maybe normalize by macMec?
}


double Genome::MeanCoverage() {
    double coverage = 0;
    for (size_t i = 0; i < haplotypes.size(); i++)
		coverage += this->haplotypes[i].MeanCoverage();
    return coverage;
}

//Expected: Nothing
//Returns: sum of mean coverages from the haplotype vector within a specific window
double Genome::WindowTotalCoverage() {
	double coverage = 0;
	for (size_t i = 0; i < this->haplotypes.size(); i++)
		coverage += this->haplotypes[i].WindowTotalCoverage();
	return coverage;
}

//Expected: nothing
//Returns: partial minimum error correction (is partial because it only calculates MEC for the window)
dnaweight_t Genome::WindowCost(dnapos_t start, dnapos_t end) {
	dnaweight_t out = 0;
	auto ploidy = this->haplotypes.size();
	for (size_t i = 0; i < ploidy; i++) {
		assert(haplotypes[i].WindowCost(start,end)>=0);
		out += haplotypes[i].WindowCost(start,end);
		assert(out>=0);
	}
	return out;
}


//Expected: nothing
//Returns: percentage of mec over maxMEC
//Question: what is maxMEC?
double Genome::Score(dnaweight_t mec) {
	double maxMec = this->haplotypes.size() * this->haplotypes[0].NumSites();
	return mec / maxMec;
}


//Expected: nothing
//Returns: nothing, however it shuffles the reads 
void Genome::Shuffle() {
	uniform_int_distribution<size_t> distribution(0, this->haplotypes.size() - 1);

	this->maxIterations = this->haplotypes[0].NumSites() * 100;

	if (this->initialized) {
		auto ploidy = this->haplotypes.size();
		auto length = this->haplotypes[0].NumSites();
		this->haplotypes.clear();
		this->haplotypes = vector<Haplotype>(ploidy, Haplotype(length, ploidy));
	}

	for (auto& r : this->file.reads) {
		this->haplotypes[distribution(this->randomEngine)].AddRead(&r);
	}

	this->initialized = true;

	// for (auto& ch : this->haplotypes) {
	// 	cout << ch << endl;
	// }
}

//Expected: updated score and the current (previous) score
//Returns: number that states whether to accept a move or not
double Genome::Acceptance(double newScore, double curScore) {
	if (newScore < curScore) return 1;
	if (this->t == 0) return 0;
	double energyDiff = curScore - newScore;
	// cout << "Acceptance(" << energyDiff << ") = " << exp(energyDiff / this->t) << endl;

	// cout << "newScore: " << newScore << ", energyDiff: " << energyDiff << ", acceptance: " << exp(energyDiff / this->t) << endl;
	return exp(energyDiff / this->t);
}

//Expected: nothing
//Returns: true or false based on if the program is done running
bool Genome::Done() {
	return this->curIteration >= this->maxIterations;
}

void Genome::SetParameters(double tInitial, double tEnd, iteration_t maxIterations) {
	if (this->curIteration) {
		// Probably warn the user in some way that the simulation is in process
	}

	this->tInitial = tInitial;
	// this->tEnd = tEnd;
	this->tDecay = -log(tEnd / this->tInitial);
	this->maxIterations = maxIterations;

	cout << "decay is " << this->tDecay << endl;
}

void Genome::SetTemperature(double t) {
	this->t = t;
}

void Genome::Move() {
	// Perform a random move, saving enough information so we can revert later

	auto ploidy = this->haplotypes.size();
	size_t moveFrom = rand() % ploidy;
	size_t moveTo;
	// FIXME1: why are we picking the source haplotype FIRST when it's possible there's no
	// reads in it? This while() loop checks for this condition and keeps looping until it finds a haplotype
	// with reads in it... but it's probably better to pick a READ from the entire universe first, the FIND
	// it's haploptype, then choose ANOTHER haplotype to move it to.
	// In any case the current code *should* work... but it's ugly.
	// It's also bad in a WINDOWING system because if the read we pick is one we shouldn't touch, we have to
	// come back all the way to here to start again... which means TWO nested while loops. Very bad.
	// FIXME1: so, pick the read FIRST.
	unsigned numTries = 0;
	while (!this->haplotypes[moveFrom].NumReads()) { // chose moveFrom (keep choosing until we find a non-empty haplotype)
		if (ploidy == 2) {
			moveFrom = !moveFrom;
		} else {
			moveFrom = rand() % ploidy;
		}
		assert(numTries++ < 10*ploidy); // this should be MORE than enough to NEVER iterate forever
	}
	// now choose where to move *to* (someplace other than from)
	if (ploidy == 2) {
		moveTo = !moveFrom;
	} else {
		size_t moveOffset = rand() % (ploidy - 1);
		moveTo = (moveFrom + moveOffset + 1) % ploidy;
	}
	assert(moveFrom != moveTo);
	// Move this "pick the read" up above FIXME1, but choose among the entire universe of reads, not just those on
	// "moveFrom". Instead, pick the read, then set moveFrom to it's current haplotype, then choose moveTo as above.
	Read * r = this->haplotypes[moveFrom].RandomRead(this->randomEngine);
	// HERE is where you put your code to check if this read shouldn't be touched because it has substantial
	// overlap with the previous window. (And you can make that decision even before looking at it's haplotype.

#if SAHAP_GENOME_DEBUG
	printf("move read %p from %d to %d\n", r, moveFrom, moveTo);
	if (!r) {
		cerr << "DEBUG: Haplotype " << moveFrom << " has no reads remaining" << endl;
		std::raise(SIGINT);
		r = this->haplotypes[moveFrom].RandomRead(this->randomEngine);
	}
#endif

	// std::cout << "From: " << moveFrom << ", To: " << moveTo << std::endl;

	this->haplotypes[moveTo].AddRead(r);
	this->haplotypes[moveFrom].RemoveRead(r);

	this->lastMove.from = moveFrom;
	this->lastMove.to = moveTo;
	this->lastMove.read = r;
}

void Genome::RevertMove() {
	const auto& move = this->lastMove;
	this->haplotypes[move.to].RemoveRead(move.read);
	this->haplotypes[move.from].AddRead(move.read);
}

void Genome::Iteration() {
	auto oldScore = this->TotalCost();
	// FIXME: these lines recomputes ALL the sites?? It should only incrementally compute the old and new scores at the sites touched by this read! Inefficient!
	this->Move();
	auto newScore = this->TotalCost();

	uniform_real_distribution<double> distribution(0, 1);
	double chanceToKeep = this->Acceptance(newScore, oldScore);
	double randomIndex = distribution(this->randomEngine);

	bool isGood = newScore < oldScore;
	bool accept = randomIndex <= chanceToKeep;
	assert(!isGood || accept);

	if (isGood) {
		// always accept
	} else if (!accept) {
		// reject
		this->RevertMove();
	}

	this->fAccept.Record(isGood);

	if (!(isGood || oldScore == newScore)) {
		this->totalBad++;
		if (accept) {
			this->totalBadAccepted++;
		}
		this->pBad.Record(chanceToKeep);
	}

#if 0
	uniform_real_distribution<double> d(0, 3);
	if (d(this->randomEngine) <= 1) {
		this->pBad.Record(false);
	} else {
		this->pBad.Record(true);
	}
#endif
}

double Genome::GetTemperature(iteration_t iteration) {
	double s = iteration / (double)this->maxIterations;
	double temp = (double)this->tInitial * exp(-this->tDecay * s);

	// cout << iteration << ": " << temp << endl;
	return temp;
}

double Genome::FracTime() {
    return this->curIteration*1.0/this->maxIterations;
}

void Genome::ResetBuffers() {
    for (iteration_t i = 0; i < this->fAccept.LENGTH; ++i) this->fAccept.buffer[i]=0;
    for (iteration_t i = 0; i < this->pBad.LENGTH; ++i) this->pBad.buffer[i]=0;
    this->pBad.len = this->fAccept.len = 0;
    this->pBad.pos = this->fAccept.pos = 0;
    this->pBad.sum = this->fAccept.sum = 0;
    this->totalBad = 0;
    this->totalBadAccepted = this->totalGood = 0;
}

void Genome::Optimize(bool debug) {
	// unsigned int TARGET_MEC = 0;//this->haplotypes[0].NumSites() * this->totalCoverage() * READ_ERROR_RATE;
	// Reset state
	this->t = this->tInitial;
	ResetBuffers();

	unsigned WINDOW_SIZE = increments * 2;
	double ERROR = READ_ERROR_RATE;
	double add = 0.0001; // FIXME: WTF is this?
	range.end = WINDOW_SIZE;
	
	for (auto& haplotype : this->haplotypes) {
		haplotype.InitializeWindow(WINDOW_SIZE, increments);
	}

	// Target MEC for the Window
	double PTARGET_MEC = WindowTotalCoverage() * ERROR;

	auto start_time = duration_cast<seconds>(system_clock::now().time_since_epoch());
	// assert(this->haplotypes.size() == 2); // FIXME need to change a few things below that assume only 0 and 1 exist.
	assert(this->haplotypes[0].NumSites() == this->haplotypes[1].NumSites());
	printf("Performing %ld meta-iterations of %d each using schedule %s,\n",
	    (long)(this->maxIterations/META_ITER), META_ITER, schedName[SCHEDULE]);
	printf("optimizing objective %s across %lu sites with total coverage %g, target MEC %g\n",
	    objName[OBJECTIVE], this->haplotypes[0].NumSites(), this->MeanCoverage(), PTARGET_MEC);

	
	int cpuSeconds = 0;
	int tmp = 0;
	
	while (range.start + WINDOW_SIZE <= numberOfSites) {
		this->t = this->GetTemperature(this->curIteration);
		this->Iteration();
		this->curIteration++;
		double pBad = this->pBad.GetAverage();
		iteration_t prev = curIteration;
		DynamicSchedule(pBad, PTARGET_MEC);
		if (debug && curIteration % REPORT_INTERVAL == 0) {
			auto now_time = duration_cast<seconds>(system_clock::now().time_since_epoch());
			cpuSeconds = (int)(now_time - start_time).count();
			Report(cpuSeconds);
			// if(FracTime() > 0.5 && pBad < 0.01 && mec() <= TARGET_MEC) {
			//     printf("Exiting early because MEC %lu reached target\n", mec());
			//     this->curIteration = this->maxIterations; // basically done
			// }
		}
		if (this->Done()){//} || (pmec() <= PTARGET_MEC)){
			range.start += increments;
			range.end += increments;
			curIteration = 0;
			tmp = cpuSeconds;
			
			for (auto& haplotype : haplotypes) {
				haplotype.IncrementWindow();
			}

			add = 0;
			PTARGET_MEC = WindowTotalCoverage() * ERROR;
		} else if (prev > curIteration) {
			add += 0.0005; // FIXME: WTF is this magic number?
			PTARGET_MEC = WindowTotalCoverage() * (ERROR + add);
		}

		if (cpuSeconds  > tmp + 50){ // For Debugging to break out if program can't find optimal solution
			break;
		}
	}
	Report(cpuSeconds, true);
	printf("Finished optimizing %d sites using %s cost function\n", (int)this->haplotypes[0].NumSites(), objName[OBJECTIVE]);
	cout << "MEC: " << TotalCost() << endl;

	// for (auto h : haplotypes) {
	// 	// h.printCoverages();
	// 	h.Print_mec();
	// }
	CreateBlocks();
}

void Genome::CreateBlocks() {
	for (Read r : file.reads) {
		if (blocks.empty() || !Intersects(blocks.back(), r.range))
			blocks.push_back(r.range);
		else
			blocks.back() = CombineBlocks(blocks.back(), r.range);
	}
}

bool Genome::Intersects(Range a, Range b) {
	return min(a.end, b.end) >= max(a.start, b.start);
}

Range Genome::CombineBlocks(Range a, Range b) {
	Range ret;
	ret.start = min(a.start, b.start);
	ret.end = max(a.end, b.end);
	return ret;
}

void Genome::Report(int cpuSeconds, bool final) {
    assert(this->TotalCost() >= 0);
    printf("%2dk (%.1f%%,%ds)  T %.3f  fA %.3f  pBad %.4f  MEC %.2f", (int)this->curIteration/1000, (100*FracTime()),
	cpuSeconds, this->t, this->fAccept.GetAverage(), this->pBad.GetAverage(), this->TotalCost());
    if (this->file.hasGroundTruth && ((cpuSeconds-lastCpuTime>1 || lastErrorRate > .25) || final || this->curIteration>=this->maxIterations)) {
	    auto gt = this->CompareGroundTruth();
	    auto ploidy = this->haplotypes.size();
	    assert(ploidy==2); // the below code only works with 2 haplotypes
	    int hapSize0=this->haplotypes[0].NumSites(),hapSize1=this->haplotypes[1].NumSites();
	    assert(hapSize0==hapSize1); // don't multiply by this->haplotypes.size()
	    double he = (double)gt / (this->haplotypes[0].NumSites() * haplotypes.size());
		this->lastErrorRate = he;
		this->lastCpuTime = cpuSeconds;
	    printf("  ( Err_vs_truth %5d Err_Pct %.2f%% [%d %d])", (int)gt, 100*he,hapSize0,hapSize1);
	    if(final) {
		printf("\nEnding ground truth ");
		if(he > 1.3 * READ_ERROR_RATE) printf("Fail!");
		else printf("Good enough");
	    }
    }
	cout << " " << range.start << "->" << range.end;
    printf("\n");
}

double Genome::FindPbad(double temperature, iteration_t iterations) {
	this->Shuffle();

	this->t = temperature;
	this->ResetBuffers();

	double sumPbad = 0.0, prevSumPbad;
	int i = 0;
	do {
		this->Iteration();
		// cout << "[simann] temp=" << this->t << ", mec=" << this->TotalCost() << endl;
		double newPbad = this->pBad.GetAverage();
		prevSumPbad = sumPbad;
		sumPbad += newPbad;
		++i;
		// below, we compute the previous and current average and demand they agree to some precision
	} while (i < 30 || fabs(prevSumPbad/(i-1) / (sumPbad/i) - 1) > 1e-3); // stabilize to at least this relative precision

	cout << "temperature " << temperature << " gives Pbad " << this->pBad.GetAverage() << " after " << i << " iterations\n";

	return this->pBad.GetAverage();
}

void Genome::AutoSchedule(iteration_t iterations) {
	cout << "Finding optimal temperature schedule..." << endl;

	double tInitial = 1;
	while (this->FindPbad(tInitial) < TARGET_PBAD_START) tInitial *= 2;
	while (this->FindPbad(tInitial) > TARGET_PBAD_START) tInitial /= 2;
	while (this->FindPbad(tInitial) < TARGET_PBAD_START) tInitial *= 1.2;
	cout << "tInitial found: " << tInitial << endl;
	double tEnd = tInitial;
	while (this->FindPbad(tEnd) > TARGET_PBAD_END) tEnd /= 2;
	while (this->FindPbad(tEnd) < TARGET_PBAD_END) tEnd *= 1.2;

	cout << "tInitial = " << tInitial << ", tEnd = " << tEnd << endl;

	this->SetParameters(tInitial, tEnd, iterations);
}

void Genome::PbadBuffer::Record(double acceptance) {
	size_t next = this->pos == LENGTH - 1 ? 0 : this->pos + 1;
	assert(acceptance >= 0.0 && acceptance <= 1.0);
	if (this->len == LENGTH) {
		this->sum -= this->buffer[next]; // the next one is the first one
	} else {
		this->len++;
	}
	this->pos = next;
	this->buffer[next] = acceptance;
	this->sum += acceptance;
	if(this->sum < 0) this->sum = 0;
}

double Genome::PbadBuffer::GetAverage() {
	if(this->len == 0) return 0.5;
	return this->sum / (double)this->len;
}

void Genome::AcceptBuffer::Record(char good) {
	size_t next = this->pos == LENGTH - 1 ? 0 : this->pos + 1;
	assert(good == 0 || good == 1);
	if (this->len == LENGTH) {
		this->sum -= this->buffer[next]; // the next one is the first one
	} else {
		this->len++;
	}
	this->pos = next;
	this->buffer[next] = good;
	this->sum += good;
}

double Genome::AcceptBuffer::GetAverage() {
	if(this->len == 0 ) return 0.5;
	return this->sum / (double)this->len;
}

dnacnt_t Genome::CompareGroundTruth() {

	if(haplotypes.size()<=0){
		return -1;
	}

	
	auto lowest = this->CompareGroundTruth(haplotypes[0], file.groundTruth[0]);
	for(unsigned i = 1; i<haplotypes.size();i++){
		lowest += this->CompareGroundTruth(haplotypes[i], file.groundTruth[i]);
	}
	vector<int> a(haplotypes.size());
	for(unsigned i =0;i<haplotypes.size();i++){
		a[i] = i;
	}

	do{
		auto temp = this->CompareGroundTruth(haplotypes[0], file.groundTruth[a[0]]);
		for(unsigned i = 1; i<haplotypes.size();i++){
			temp += this->CompareGroundTruth(haplotypes[i], file.groundTruth[a[i]]);
		}

		if(lowest > temp){
			lowest = temp;
		}
	}while(next_permutation(a.begin(), a.end()));


	return lowest;
}

dnacnt_t Genome::CompareGroundTruth(const Haplotype& ch, const vector<int>& truth) {
	dnacnt_t loss = 0;
	for (size_t i = 0; i < ch.NumSites(); ++i) {
		if (ch.solution[i] != truth[i]) {
			loss++;
		}
	}
	return loss;
}

ostream& operator << (ostream& stream, const Genome& ge) {
	int blockNum = 1;
	for (Range r : ge.blocks) {
		stream << "BLOCK " << blockNum++ << endl;
		for (size_t i = 0; i < ge.haplotypes.size(); ++i) {
			for (size_t j = 0; j < ge.haplotypes[i].NumSites(); ++j) {
				if (j >= r.start && j <= r.end) {
					if (ge.haplotypes[i].solution[j] < 0)
						stream << 'X';
					else
						stream << ge.haplotypes[i].solution[j];
				}
				else
					stream << '-';
			}
			stream << endl;
		}
	}
	return stream;
}

void Genome::DynamicSchedule(double pBad, double TARGET_MEC)
{
    //printf("pBad %ld %g\n", this->pBad.len, pBad);
#if SCHEDULE==RETREAT
/*
This schedule seems to be about right... up to 1500 SNPs. Then it fails.  Observations:
- There are a lot of magic numbers here. I think it may be possible to simplify. Ideas:
- it looks like maintaing a ration of about 10x as many 1% retreats as 95% retreats; automating that may help.
- the locations in time (fracTime and pBad) do NOT seem very flexible; you can play with what tolerances
you want at those times/pBad's, but if you shift the times or pBads too much, it all falls apart.
- decreasing the threshold (ie., tightening the tolerance) on the 1% retreats can reduce the 95% retreats,
but if you tighten it too much it never gets past 30% fracTime.
- having the first threhold too loose results in way to many 95% retreats. (ie complete restart)
- oddly, increasing the number of meta-iterations on the command line does not always help; then the retreats
are all mis-calibrated, and we never get past 30% fracTime. So again, something to automatically calibrate
them seems necessary.
- there doesn't seem to be much of a difference between optimizing MEC vs Possion. Poisson is a bit slower but
doesn't seem to give any better aswers.
- the dead-ends seem REALLY tight: sometimes you can still have 90% distance from ground truth even with the MEC
being less than a factor for two from the TARGE_MEC!!! That seems *really* strange. Is that even possible?
Perhaps the data is fucked up? It just seems so unlikely the MEC can be almost to the target while we're
still TOTALLY wrong compared to ground truth. For example, optimizing 2100SNP case, target was 614:
598k (93.4%,74s)  T 0.000  fA 0.004  pBad 0.001  MEC 850 ( Err_vs_truth 1955 Err_Pct 93.10% [2100 2100])
599k (93.6%,74s)  T 0.000  fA 0.000  pBad 0.001  MEC 851 ( Err_vs_truth 1954 Err_Pct 93.05% [2100 2100])
600k (93.8%,74s)  T 0.000  fA 0.000  pBad 0.001  MEC 852 ( Err_vs_truth 1956 Err_Pct 93.14% [2100 2100])
601k (93.9%,74s)  T 0.000  fA 0.000  pBad 0.001  MEC 852 ( Err_vs_truth 1956 Err_Pct 93.14% [2100 2100])
602k (94.1%,74s)  T 0.000  fA 0.001  pBad 0.001  MEC 852 ( Err_vs_truth 1954 Err_Pct 93.05% [2100 2100])
Retreat 94% from 94.0625% to 0.0625% because MEC is 852, too big by factor of 1.38762
We're only a factor of 1.387 above the target MEC and we're still 93% wrong compared to the ground truth??? WTF?

- again with ground truth... sometimes we converge to exactly the right ground truth yet the MEC is still large.
Again, how is that possible?  For example, here's a sequence of final lines to a good solution: (target was 268)
49k (61.3%,159s)  T 0.002  fA 0.009  pBad 0.011  MEC   264  ( Err_vs_truth     0 Err_Pct 0.00% [1000 1000])
50k (62.5%,159s)  T 0.002  fA 0.015  pBad 0.015  MEC   272  ( Err_vs_truth     0 Err_Pct 0.00% [1000 1000])
51k (63.7%,159s)  T 0.002  fA 0.011  pBad 0.009  MEC   260  ( Err_vs_truth     0 Err_Pct 0.00% [1000 1000])
How can the MEC continue to decrease while the Err_vs_truth remains exactly zero? No fucking way. There is a
unique solution---all the reads on the correct side. That solution has a fixed MEC. Can't have one change while
the other remain constant.
*/
    double retreat = 0.0;
    static double prev_retreat_frac;
    int num_meta_iters = this->maxIterations/META_ITER;
#define SMALL_RETREAT 0.01 // Let it grow with number of meta-iters? (0.01*(1+2*log(num_meta_iters)))
#define FULL_RETREAT 0.94 // this needs to be less than (1-(REPORT_INTERVAL/2)) from the next line
    if(curIteration % (REPORT_INTERVAL/2) == 0) {
	double factor = (double)TotalCost()/TARGET_MEC;
	// double factor = (double)pmec()/(TARGET_MEC == 0 ? 0.5 : TARGET_MEC);
	if(FracTime() - prev_retreat_frac > 2*SMALL_RETREAT &&
	    (((FracTime()>0.3||pBad<0.2) && factor > 16) ||   // 14 to 22 seems to work well
	     ((FracTime()>0.5||pBad<0.1) && factor >  8) )){  // quarter to half the above works well?
	    retreat = factor * SMALL_RETREAT / MeanCoverage() * log(num_meta_iters);
	}
	if(FracTime()>FULL_RETREAT && factor > 1.3) {//pmec() != 0){
	    retreat = FULL_RETREAT;
	    ResetBuffers();
	}
	if(retreat > 0.0) {
	    cout << "Retreat " << 100*retreat << "% from " << 100 * FracTime();
	    this->curIteration -= retreat * this->maxIterations;
	    assert(this->curIteration>=0);
	    cout << "% to "  << 100 * FracTime() << "% because MEC is " << TotalCost();
	    cout << ", too big by a factor of " << factor << "(" << TARGET_MEC << 
			", " << WindowTotalCoverage() << ")" << endl;
	    prev_retreat_frac = FracTime();
	}
    }
#elif SCHEDULE==Betz
    static double computedTdecay, LOWER=atof(getenv("LOWER")),
	ACCEPT_TARGET=atof(getenv("ACCEPT_TARGET")), PBAD_TARGET=atof(getenv("PBAD_TARGET"));
    if(!computedTdecay){
	computedTdecay = this->tDecay;
	printf("Betz values: LOWER %g ACCEPT_TARGET %g PBAD_TARGET %g\n",LOWER,ACCEPT_TARGET,PBAD_TARGET);
    }
    if(this->curIteration > this->fAccept.LENGTH) {
	this->tDecay = computedTdecay * (LOWER +
	    min(fabs(ACCEPT_TARGET-this->fAccept.GetAverage()), fabs(PBAD_TARGET-this->pBad.GetAverage())));
    }
#endif
}

}