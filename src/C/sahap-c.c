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
    HAPLOTYPE haps[PLOIDY];
    int totalMEC;
    short *MEC; // the current error at each site (per-site, summed across haplotypes)
} GENOME;


static READ _read[NUM_READS]; // the global list of reads
static SITE _site[GENOME_LEN]; // we know which reads touch each site


// Computes the whole-genome MEC at a particular site by summing across the per-Haplotype MECs.
// On the way we also compute and assign the solution at this site for each haplotpye
// It MIGHT have made more sense to have a function that explicitly computes the MEC only for 1 haplotype,
// but I don't want to change that right now.
int ComputeSiteMEC(GENOME *G, SITE *site) {
    unsigned readsThatTouch[NUM_READS]; // outrageously over-sized array of reads that touch this site
    unsigned numReads = SetToArray(readsThatTouch, site->readsThatTouch);
    int count[PLOIDY][NUM_LETTERS], not[PLOIDY][NUM_LETTERS]; // count of letter[h][X] and disagreeing count
    unsigned h, i, j, r, numTouch;
    for(h=0;h<PLOIDY;h++) for(i=0;i<NUM_LETTERS;i++) count[h][i]=not[h][i]=0;
    for(i=0; i<numReads; i++) {
	READ *r = &_read[readsThatTouch[i]]; // get pointer to each read (globally) that touches this site
	assert(0 <= r->hap && r->hap < PLOIDY);
	assert(r->start <= site->location && site->location < r->top);
	int readLoc = site->location - r->start;
	int letter = r->let[readLoc];
	assert(0 <= letter && letter < NUM_LETTERS);
	++count[r->hap][letter];
	for(j=0;j<NUM_LETTERS;j++) if(j!=letter) ++not[r->hap][j];
    }
    int siteMEC = 0;
    for(h=0;h<PLOIDY;h++) {
	int maxCount=-1, maxLet=-1;
	for(i=0;i<NUM_LETTERS;i++) if(count[h][i] >= maxCount) {maxCount=count[h][i]; maxLet=i;}
	G->haps[h].sol[site->location] = maxLet;
	G->haps[h].MEC[site->location] = not[h][maxLet];
	siteMEC += G->haps[h].MEC[site->location];
    }
    return siteMEC;
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
	_read[r].start = drand48()*(GENOME_LEN-READ_LEN+1);
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

void Report(GENOME *G) {
#if VERBOSE>1
    for(int h=0; h<PLOIDY; h++) printf("Haplotype %d currently has %d reads\n", h, SetCardinality(G->haps[h].readSet));
#endif
    int genomeMEC=0;
    for(int i=0; i<GENOME_LEN; i++) {
	genomeMEC += ComputeSiteMEC(G, &_site[i]);
#if VERBOSE>1
	printf("site %d touches %d reads\n", i, SetCardinality(_site[i].readsThatTouch));
	for(int h=0;h<PLOIDY;h++) {
	    static SET *intersect;
	    if(!intersect) intersect=SetAlloc(NUM_READS); else SetReset(intersect);
	    SetIntersect(intersect, G->haps[h].readSet, _site[i].readsThatTouch);
	    printf("\thap[%d] touches %d reads, is majority %d, with MEC %d\n", h, SetCardinality(intersect),
		G->haps[h].sol[i], G->haps[h].MEC[i]);
	}
#endif
    }
    printf("\ngenomeMEC %d", genomeMEC);
}

void FlipRead(READ *r) {

}

void HillClimb(GENOME *G) {
    int maxMEC=0, iter=0, stagnant=0;
    // Note: the "stagnant" variable isn't needed in SA, it's only needed in Hill Climbing, because if we've
    // tried 1000 moves without any improvement, it's probably time to give up because we're at a local minimum.
    while(stagnant<1000) { // for SA, this "while" will be replace with a loop over the temperature range estimated first by SA
	++iter;
	Report(G);
	int r=drand48()*NUM_READS; // pick a read at random
	int hap=_read[r].hap;
	assert(hap==0 || hap==1); // C code limitation: the negation below only works if PLOIDY==2
#if VERBOSE
	printf(" Read %d is current in H%d...", r, hap);
#endif
	int beforeMEC=0, afterMEC=0;
	for(int i=_read[r].start; i<_read[r].top;i++) beforeMEC += ComputeSiteMEC(G, &_site[i]);
	// Now flip it and recompute the MEC along its sites (assumes PLOIDY==1)
	SetDelete(G->haps[hap].readSet, r);
	SetAdd  (G->haps[!hap].readSet, r);
	_read[r].hap=!hap;
	for(int i=_read[r].start; i<_read[r].top;i++) afterMEC += ComputeSiteMEC(G, &_site[i]);
	if(beforeMEC>maxMEC || afterMEC>maxMEC) maxMEC=MAX(beforeMEC,afterMEC);
#if VERBOSE
	printf("before %d, after %d...", beforeMEC, afterMEC);
#endif
	// Here is where SA might choose to accept a bad move rather than rejecting all bad moves.
	if(afterMEC < beforeMEC) {
	    stagnant=0;
#if VERBOSE
	    printf("accept ") ; // accept the move (do nothing, it's already moved)
#endif
	}
	else { // reject
	    ++stagnant;
#if VERBOSE
	    printf("reject(%d) ", stagnant);
#endif
	    afterMEC=0;
	    _read[r].hap=hap;
	    SetDelete(G->haps[!hap].readSet, r);
	    SetAdd    (G->haps[hap].readSet, r);
	    for(int i=_read[r].start; i<_read[r].top;i++) afterMEC += ComputeSiteMEC(G, &_site[i]);
	    assert(afterMEC == beforeMEC);
	}
    }
    printf("\nHill Climbing stagnated\n");
}

int main(int argc, char *argv)
{
    srand48(time(NULL)+getpid());

    CreateReads(); // reads "created" from the true (but unknown here) haplotype set

    GENOME G;

    for(int h=0; h<PLOIDY; h++) {
	G.haps[h].id=h;
	G.haps[h].totalMEC=0;
	G.haps[h].readSet = SetAlloc(NUM_READS);
	G.haps[h].sol = Calloc(GENOME_LEN,sizeof(G.haps[h].sol[0])); // computed from the most common element
	G.haps[h].MEC = Calloc(GENOME_LEN,sizeof(G.haps[h].MEC[0])); // number of differences from most common element
    }

    // Initially, assign each read to a random haplotype
    for(int r=0; r<NUM_READS; r++) {
	int hap = _read[r].hap = PLOIDY*drand48();
	SetAdd(G.haps[hap].readSet, r);
    }

    printf("Genome len %d, coverage %d, %d reads of length %d, numHap %d\n", GENOME_LEN, COVERAGE, NUM_READS, READ_LEN, PLOIDY);
    HillClimb(&G);
}
