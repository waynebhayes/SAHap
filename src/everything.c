#include "misc.h"
#include "sets.h"
#include "rand48.h"

#define KNOWN_TRUTH 1
#define GENOME_LEN 1000
#define READ_LEN 10 // constant for now, for simplicity
#define COVERAGE 20
#define PLOIDY 2
#define NUM_LETTERS 2 // number of possible letters in the solution
#define NUM_READS (GENOME_LEN / READ_LEN * COVERAGE)

typedef struct _read {
    int start, top; // this read runs from start to (top-1), ie., top is NOT included
    int id; // starts at 0 and ends at (numReads-1).
    char *val; // n=(top-start), val[n] will be allocated
    char hapSol; // which haplotype are we currently assigned to? (initially -1)
#if KNOWN_TRUTH
    char trueHap;
#endif
} READ;

typedef struct _site {
    int location;
    SET *readsThatTouch; // Note this could be a list for a HAPLOTYPE, or for a whole GENOME!
} SITE;

typedef struct _haplotype {
    int id, totalMEC;
    SET *readSet; // currently assigned reads
    char *sol; // the current solution from start to GENOME_LEN-1, according to the majority vote, per-site
    char *MEC; // the current error at each site (ie., number of reads that disagree with the solution)
} HAPLOTYPE;

typedef struct _genome {
    int len; // how long is the genome?
    int totalMEC;
    short *MEC; // the current error at each site (per-site, summed across haplotypes)
} GENOME;


static READ _read[NUM_READS]; // the global list of reads
static SITE _site[GENOME_LEN]; // we know which reads touch each site

// per-site functions

char HaplotypeComputeSiteMajority(HAPLOTYPE *h, SITE *site, int *pNumTouch) {
    unsigned readsThatTouch[NUM_READS]; // outrageously over-sized array of reads that touch this site
    unsigned numReads = SetToArray(readsThatTouch, site->readsThatTouch);
    int count[NUM_LETTERS]; // keep count to find majority later
    unsigned i, r;
    *pNumTouch=0;
    for(i=0;i<NUM_LETTERS;i++) count[i]=0;
    for(i=0; i<numReads; i++) {
	READ *r = &_read[readsThatTouch[i]]; // get pointer to each read (globally) that touches this site
	assert(0 <= r->hapSol && r->hapSol < PLOIDY);
	if(r->hapSol == h->id) { // would be more efficient to have per-hap readsThatTouch, but this'll do for now
	    ++*pNumTouch;
	    assert(r->start <= site->location && site->location < r->top);
	    int readLoc = site->location - r->start;
	    ++count[r->val[readLoc]];
	}
    }
    if(*pNumTouch==0) return -1; // haplotype has no reads at this site
    int maxCount=-1, maxVal=-1;
    for(i=0;i<NUM_LETTERS;i++) if(count[i] > maxCount) {maxCount=count[i]; maxVal=i;}
    assert(maxCount>=0 && maxVal>=0);
    return maxVal;
}


int HaplotypeComputeSiteMEC(HAPLOTYPE *h, SITE *site) {
    int numTouch;
    char majority = HaplotypeComputeSiteMajority(h, site, &numTouch);
    if(majority<0) return 0; // no reads touch this site at all, so MEC is 0
    unsigned readsThatTouch[NUM_READS]; // outrageously over-sized array of reads that touch this site
    unsigned numReads = SetToArray(readsThatTouch, site->readsThatTouch);
    int i, MEC=0;
    for(i=0; i<numReads; i++) {
	READ *r = &_read[readsThatTouch[i]]; // get pointer to each read (globally) that touches this site
	assert(0 <= r->hapSol && r->hapSol < PLOIDY);
	if(r->hapSol == h->id) { // would be more efficient to have per-hap readsThatTouch, but this'll do for now
	    assert(r->start <= site->location && site->location < r->top);
	    int readLoc = site->location - r->start;
	    if(r->val[readLoc] != majority) ++MEC;
	}
    }
    return MEC;
}


int main(int argc, char *argv)
{
    srand48(time(NULL)+getpid());
    // Generate the true haplotypes
    HAPLOTYPE trueHap[PLOIDY];
    HAPLOTYPE hapSol[PLOIDY]; // the current "solution"

    for(int h=0; h<PLOIDY; h++) {
	// the true haplotype
	trueHap[h].id=h;
	trueHap[h].readSet = SetAlloc(NUM_READS);
	trueHap[h].sol = Calloc(GENOME_LEN,sizeof(trueHap[h].sol[0]));
	trueHap[h].MEC = Calloc(GENOME_LEN,sizeof(trueHap[h].MEC[0]));
	for(int i=0; i<GENOME_LEN; i++)
	    trueHap[h].sol[i] = NUM_LETTERS*drand48();

	// the "solution" (computed) haplotype
	hapSol[h].id=h;
	hapSol[h].readSet = SetAlloc(NUM_READS);
	hapSol[h].sol = Calloc(GENOME_LEN,sizeof(trueHap[h].sol[0]));
	hapSol[h].MEC = Calloc(GENOME_LEN,sizeof(trueHap[h].MEC[0]));
    }

    // allocate the sites and initialize to EMPTY the set of reads that touch each site
    for(int i=0; i<GENOME_LEN; i++) {
	_site[i].location = i;
	_site[i].readsThatTouch = SetAlloc(NUM_READS);
    }

    // Generate the reads
    for(int r=0; r<NUM_READS; r++) {
	_read[r].id=r;
	_read[r].start = drand48()*(GENOME_LEN-READ_LEN);
	_read[r].top = _read[r].start+READ_LEN;
	assert(_read[r].top <= GENOME_LEN);
	_read[r].val = Calloc(READ_LEN, sizeof(_read[r].val[0]));
	int hap = PLOIDY*drand48(); // the TRUE haplotype this read is from
#if KNOWN_TRUTH
        _read[r].trueHap = hap;
#endif
	SetAdd(trueHap[hap].readSet, r); // record this read's true haplotype
	for(int i=0; i<READ_LEN; i++) {
	    int location = _read[r].start + i;
	    _read[r].val[i] = trueHap[hap].sol[location];
	    SetAdd(_site[location].readsThatTouch, r); // record the list of sites this read touches.
	}
	_read[r].hapSol = PLOIDY*drand48(); // record this read's RANDOMLY assigned haplotype
	SetAdd(hapSol[_read[r].hapSol].readSet, r);
    }

    printf("Genome len %d, coverage %d, numReads %d, numHap %d\n", GENOME_LEN, COVERAGE, NUM_READS, PLOIDY);
    for(int h=0; h<PLOIDY; h++) printf("Haplotype %d currently has %d reads\n", h, SetCardinality(trueHap[h].readSet));

    for(int i=0; i<GENOME_LEN; i++) {
	printf("site %d touches %d reads\n", i, SetCardinality(_site[i].readsThatTouch));
	for(int j=0;j<PLOIDY;j++) {
	    int numTouch, majority;
	    majority = HaplotypeComputeSiteMajority(&trueHap[j], &_site[i], &numTouch),
	    printf("\thap[%d] touches %d reads, is majority %d, with MEC %d\n", j, numTouch, majority,
		HaplotypeComputeSiteMEC(&trueHap[j], &_site[i]));
	}
    }
}
