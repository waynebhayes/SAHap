#include "misc.h"
#include "sets.h"
#include "rand48.h"

#define VERBOSE 1
#define KNOWN_TRUTH 1
#define GENOME_LEN 500
#define READ_LEN 5 // constant for now, for simplicity
#define COVERAGE 30
#define PLOIDY 2
#define NUM_LETTERS 2 // number of possible letters in the solution
#define NUM_READS (GENOME_LEN / READ_LEN * COVERAGE)

typedef struct _read {
    int id; // the "id" is just the subscript in the _read[] array.
    int start, top; // this read runs from start to (top-1), ie., top is NOT included
    char *let; // the array of READ_LEN letters
    char hap; // which haplotype are we currently assigned to? (initially -1)
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
    assert(0<=h->id && h->id<PLOIDY);
    unsigned readsThatTouch[NUM_READS]; // outrageously over-sized array of reads that touch this site
    unsigned numReads = SetToArray(readsThatTouch, site->readsThatTouch);
    int count[NUM_LETTERS]; // keep count to find majority later
    unsigned i, r;
    *pNumTouch=0;
    for(i=0;i<NUM_LETTERS;i++) count[i]=0;
    for(i=0; i<numReads; i++) {
	READ *r = &_read[readsThatTouch[i]]; // get pointer to each read (globally) that touches this site
	assert(0 <= r->hap && r->hap < PLOIDY);
	if(r->hap == h->id) { // would be more efficient to have per-hap readsThatTouch, but this'll do for now
	    ++*pNumTouch;
	    assert(r->start <= site->location && site->location < r->top);
	    int readLoc = site->location - r->start;
	    int letter = r->let[readLoc];
	    assert(0 <= letter && letter < NUM_LETTERS);
	    ++count[letter];
	}
    }
    if(*pNumTouch==0) return -1; // haplotype has no reads at this site
    int maxCount=-1, maxLet=-1;
    for(i=0;i<NUM_LETTERS;i++) if(count[i] >= maxCount) {maxCount=count[i]; maxLet=i;}
    assert(0<=maxCount && maxCount<=*pNumTouch && 0<=maxLet && maxLet<NUM_LETTERS);
    return maxLet;
}


int HaplotypeComputeSiteMEC(HAPLOTYPE *h, SITE *site, int *pNumTouch) {
    char majority = HaplotypeComputeSiteMajority(h, site, pNumTouch);
    if(majority<0) return 0; // no reads touch this site at all, so MEC is 0
    unsigned readsThatTouch[NUM_READS]; // outrageously over-sized array of reads that touch this site
    unsigned numReads = SetToArray(readsThatTouch, site->readsThatTouch);
    int i, MEC=0;
    for(i=0; i<numReads; i++) {
	READ *r = &_read[readsThatTouch[i]]; // get pointer to each read (globally) that touches this site
	assert(0 <= r->hap && r->hap < PLOIDY);
	if(r->hap == h->id) { // would be more efficient to have per-hap readsThatTouch, but this'll do for now
	    assert(r->start <= site->location && site->location < r->top);
	    int readLoc = site->location - r->start;
	    int letter = r->let[readLoc];
	    if(letter != majority) ++MEC;
	}
    }
    return MEC;
}

// Create the ground truth genome + haplotypes (though the haplotypes will be erased after this function call).
// Then generate the reads: each comes from a random (true) haplotype at a certain location.
static void CreateReads(void) {
    HAPLOTYPE trueHap[PLOIDY]; // this will go away after the function call
    assert(PLOIDY==2);
    for(int h=0; h<PLOIDY; h++)
	trueHap[h].sol = Calloc(GENOME_LEN,sizeof(trueHap[h].sol[0]));
    for(int h=0; h<PLOIDY/2; h++) {
	for(int i=0; i<GENOME_LEN; i++) {
	    trueHap[h].sol[i] = NUM_LETTERS*drand48();
	    trueHap[h+1].sol[i] = (trueHap[h].sol[i]+1) % NUM_LETTERS;
	}
    }

    // allocate the sites and initialize to EMPTY the set of reads that touch each site
    for(int i=0; i<GENOME_LEN; i++) {
	_site[i].location = i;
	_site[i].readsThatTouch = SetAlloc(NUM_READS); // initially empty set
    }

    // Generate the reads
    for(int r=0; r<NUM_READS; r++) {
	_read[r].id=r;
	_read[r].start = drand48()*(GENOME_LEN-READ_LEN);
	assert(_read[r].start >=0);
	_read[r].top = _read[r].start+READ_LEN;
	assert(_read[r].top <= GENOME_LEN);
	_read[r].let = Calloc(READ_LEN, sizeof(_read[r].let[0]));
	int hap = PLOIDY*drand48(); // the TRUE haplotype this read is from
#if KNOWN_TRUTH
        _read[r].trueHap = hap;
#endif
#if VERBOSE>1
	printf("R%d[%d,%d)<-H%d[",r,_read[r].start,_read[r].top,hap);
#endif
	for(int i=0; i<READ_LEN; i++) {
	    int location = _read[r].start + i;
	    _read[r].let[i] = trueHap[hap].sol[location];
#if VERBOSE>1
	printf("%d",_read[r].let[i]);
#endif
	    SetAdd(_site[location].readsThatTouch, r); // record the list of sites this read touches.
	}
#if VERBOSE>1
	printf("]\n");
#endif
    }

    // Now nuke the truth, we're not allowed to know it.
    for(int h=0; h<PLOIDY; h++) {
	Free(trueHap[h].sol);
    }
}

int main(int argc, char *argv)
{
    srand48(time(NULL)+getpid());

    CreateReads(); // reads "created" from the true (but unknown here) haplotype set

    HAPLOTYPE hapSol[PLOIDY]; // These are the constructed haplotypes, initially empty.
    for(int h=0; h<PLOIDY; h++) {
	hapSol[h].id=h;
	hapSol[h].totalMEC=0;
	hapSol[h].readSet = SetAlloc(NUM_READS);
	hapSol[h].sol = Calloc(GENOME_LEN,sizeof(hapSol[h].sol[0])); // computed from the most common element
	hapSol[h].MEC = Calloc(GENOME_LEN,sizeof(hapSol[h].MEC[0])); // number of differences from most common element
    }

    // Initially, assign each read to a random haplotype
    for(int r=0; r<NUM_READS; r++) {
	int hap = _read[r].hap = PLOIDY*drand48();
	SetAdd(hapSol[hap].readSet, r);
    }

    // Initialize the solutions and MEC values
    for(int i=0; i<GENOME_LEN; i++) {
	for(int h=0; h<PLOIDY; h++) {
	    int numTouch, majority, MEC;
	    hapSol[h].sol[i] = HaplotypeComputeSiteMajority(&hapSol[h], &_site[i], &numTouch);
	    hapSol[h].MEC[i] = HaplotypeComputeSiteMEC(&hapSol[h], &_site[i], &numTouch);
	}
    }
#if VERBOSE
    printf("Genome len %d, coverage %d, numReads %d, numHap %d\n", GENOME_LEN, COVERAGE, NUM_READS, PLOIDY);
    for(int h=0; h<PLOIDY; h++) printf("Haplotype %d currently has %d reads\n", h, SetCardinality(hapSol[h].readSet));
#endif
    for(int i=0; i<GENOME_LEN; i++) {
	printf("site %d touches %d reads\n", i, SetCardinality(_site[i].readsThatTouch));
	for(int h=0;h<PLOIDY;h++) {
	    int numTouch, majority = HaplotypeComputeSiteMajority(&hapSol[h], &_site[i], &numTouch);
	    assert(majority == hapSol[h].sol[i]);
	    printf("\thap[%d] touches %d reads, is majority %d, with MEC %d\n", h, numTouch, majority, hapSol[h].MEC[i]);
	}
    }
}
