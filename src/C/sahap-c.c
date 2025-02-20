#include "misc.h"
#include "sets.h"
#include "rand48.h"

#define VERBOSE 1
#define MAX_NUM_SITES 2000000
#define MAX_READ_LEN 400 // constant for now, for simplicity
#define COVERAGE 20
#define PLOIDY 2
#define NUM_LETTERS 2 // number of possible letters in the solution
#define MAX_NUM_READS 15000

typedef struct _site {
    int globalLoc; // NEVER use this global location except when reading in a WIF file--use the WhichSite instead
    SET *readsThatTouch; // Note this could be a list for a HAPLOTYPE, or for a whole GENOME!
} SITE;

typedef struct _haplotype {
    int id, totalMEC;
    SET *readSet; // currently assigned reads
    char *sol; // the current solution from start to MAX_NUM_SITES-1, according to the majority vote, per-site
    char *MEC; // the current error at each site (ie., number of reads that disagree with the solution)
} HAPLOTYPE;

typedef struct _genome {
    HAPLOTYPE haps[PLOIDY];
    int totalMEC;
    short *MEC; // the current error at each site (per-site, summed across haplotypes)
} GENOME;


typedef struct _read {
    int id; // the "id" is just the subscript in the _read[] array.
    int firstSite, numSites; // this read runs from start to (top-1), ie., top is NOT included
    char *let; // the array of nucleotides (or bits)
    char hap, trueHap; // currently assigned haplotype, and ground truth (-1 if unknown)
} READ;

static READ _read[MAX_NUM_READS]; // the global list of reads
static SITE _site[MAX_NUM_SITES]; // we know which reads touch each site
static int _numReads; // number of reads actually read in

static int _numSites;

// given a global genome location, which site is it? -1 means "not found"
static int WhichSite(int loc) {
    assert(loc>0); // we can't handle the true location being zero
    int i; // FIXME: linear search for now, maybe binary search when this gets big...
    for(i=0; i<_numSites;i++) if(_site[i].globalLoc==loc) return i;
    return -1;
}


void ReadWIF(char filename[]) {
    FILE *fp = fopen(filename, "r");
    int readSite=0, siteLoc, conf;
    char nuc, bit;
    char siteLocStr[20]; // integer locations in a genome have at most... 10? characters
    while(fscanf(fp, "%s ", siteLocStr)==1) {
	siteLoc = atoi(siteLocStr);
	if(siteLoc == 0) { // we've reached the end of a read
	    assert(siteLocStr[0] != '0'); // oops, the siteLoc REALLY IS zero... FIXME
	    char junk[1000]; char *foo = fgets(junk, sizeof(junk), fp); assert(foo); // read through the comment
	    if(VERBOSE>1) printf("Finished reading Read %d (%d sites)\n", _numReads, readSite);
	    _read[_numReads].numSites = readSite;
	    assert(_numReads < MAX_NUM_READS);
	    _numReads++; readSite=0;
	}
	else {
	    int whichSite = WhichSite(siteLoc);
	    if(whichSite == -1) { // it's a new siteLoc we haven't yet seen
		assert(_numSites < MAX_NUM_SITES);
		whichSite=_numSites;
		_site[_numSites].globalLoc = siteLoc;
		_numSites++;
	    }
	    int n = fscanf(fp, "%c %c %d : ", &nuc, &bit, &conf);
	    if(VERBOSE>1) printf("n=%d, siteLoc=%d, nuc=%c, bit=%c conf=%d\n", n, siteLoc, nuc, bit, conf);
	    if(readSite==0) { // this is the first site of a new read
		if(VERBOSE>1) printf("Found beginning of read %d (line %d), siteLoc %d\n", _numReads, _numReads+1, siteLoc);
		_read[_numReads].id = _numReads;
		_read[_numReads].firstSite = whichSite;
		_read[_numReads].let = Malloc(MAX_READ_LEN);
		_read[_numReads].trueHap = -1;
	    }
	    assert(readSite < MAX_READ_LEN);
	    _read[_numReads].let[readSite++]=bit-'0';
	}
    }
    printf("Read %d reads across %d sites\n", _numReads, _numSites);
    fclose(fp);
}

// Create the ground truth genome + haplotypes (though the haplotypes will be erased after this function call).
// Then generate the reads: each comes from a random (true) haplotype at a certain location.
static void CreateRandomReads(void) {
    HAPLOTYPE trueHap[PLOIDY]; // this will go away after the function call
    for(int h=0; h<PLOIDY; h++)
	trueHap[h].sol = Calloc(MAX_NUM_SITES,sizeof(trueHap[h].sol[0]));
    // When PLOIDY==2, the letters are ALWAYS opposites at the same site... but I don't know what happens if PLOIDY>2
    assert(PLOIDY==2);
    for(int h=0; h<PLOIDY/2; h++) {
	for(int i=0; i<MAX_NUM_SITES; i++) {
	    trueHap[h].sol[i] = NUM_LETTERS*drand48();
	    trueHap[h+1].sol[i] = (trueHap[h].sol[i]+1) % NUM_LETTERS;
	}
    }

    // allocate the sites and initialize to EMPTY the set of reads that touch each site
    _numSites = MAX_NUM_SITES;
    for(int i=0; i<_numSites; i++) _site[i].globalLoc = i;

    // Generate the reads
    _numReads = MAX_NUM_READS;
    for(int r=0; r<_numReads; r++) {
	_read[r].id=r;
	_read[r].firstSite = drand48()*(MAX_NUM_SITES-MAX_READ_LEN+1);
	assert(_read[r].firstSite >= 0);
	_read[r].numSites = MAX_READ_LEN;
	_read[r].let = Calloc(MAX_READ_LEN, sizeof(_read[r].let[0]));
	int hap = PLOIDY*drand48(); // the TRUE haplotype this read is from
        _read[r].trueHap = hap;
	if(VERBOSE>1) printf("R%d[%d=%d]<-H%d[",r,_read[r].firstSite,_read[r].numSites,hap);
	for(int i=0; i<MAX_READ_LEN; i++) {
	    int site = _read[r].firstSite + i;
	    _read[r].let[i] = trueHap[hap].sol[site];
	    if(VERBOSE>1) printf("%d",_read[r].let[i]);
	}
	if(VERBOSE>1) printf("]\n");
    }
}

static void InitializeSystem(void) {
    // allocate the sites and initialize to EMPTY the set of reads that touch each site
    int i, r;
    for(i=0; i<_numSites; i++) _site[i].readsThatTouch = SetAlloc(_numReads); // initially empty set
    for(r=0; r<_numReads; r++) {
	_read[r].hap = PLOIDY*drand48(); // initialized to random
	for(i=0; i<_read[r].numSites; i++) SetAdd(_site[_read[r].firstSite+i].readsThatTouch, r);
    }
}

// Computes the whole-genome MEC at a particular site by summing across the per-Haplotype MECs.
// On the way we also compute and assign the solution at this site for each haplotpye
// It MIGHT have made more sense to have a function that explicitly computes the MEC only for 1 haplotype,
// but I don't want to change that right now.
int ComputeSiteMEC(GENOME *G, SITE *site) {
    unsigned readsThatTouch[MAX_NUM_READS]; // outrageously over-sized array of reads that touch this site
    unsigned numReads = SetToArray(readsThatTouch, site->readsThatTouch);
    int count[PLOIDY][NUM_LETTERS], not[PLOIDY][NUM_LETTERS]; // count of letter[h][X] and disagreeing count
    unsigned h, i, j, r, numTouch;
    for(h=0;h<PLOIDY;h++) for(i=0;i<NUM_LETTERS;i++) count[h][i]=not[h][i]=0;
    for(i=0; i<numReads; i++) {
	READ *r = &_read[readsThatTouch[i]]; // get pointer to each read (globally) that touches this site
	assert(0 <= r->hap && r->hap < PLOIDY);
	int siteIndex = site - _site;
	int readLoc = siteIndex - r->firstSite; assert(0 <= readLoc && readLoc < r->numSites);
	int letter = r->let[readLoc];
	assert(0 <= letter && letter < NUM_LETTERS);
	++count[r->hap][letter];
	for(j=0;j<NUM_LETTERS;j++) if(j!=letter) ++not[r->hap][j];
    }
    int siteMEC = 0;
    for(h=0;h<PLOIDY;h++) {
	int maxCount=-1, maxLet=-1;
	for(i=0;i<NUM_LETTERS;i++) if(count[h][i] >= maxCount) {maxCount=count[h][i]; maxLet=i;}
	G->haps[h].sol[site-_site] = maxLet;
	G->haps[h].MEC[site-_site] = not[h][maxLet];
	siteMEC += G->haps[h].MEC[site-_site];
    }
    return siteMEC;
}


void Report(int iter, GENOME *G) {
    if(VERBOSE) {
	printf("iter %d; Read counts:", iter);
	for(int h=0; h<PLOIDY; h++)
	    printf(" H%d=%d", h, SetCardinality(G->haps[h].readSet));
    }

    int genomeMEC=0;
    for(int i=0; i<_numSites; i++) {
	genomeMEC += ComputeSiteMEC(G, &_site[i]);
	if(VERBOSE>1) {
	    printf("site %d touches %d reads\n", i, SetCardinality(_site[i].readsThatTouch));
	    for(int h=0;h<PLOIDY;h++) {
		static SET *intersect;
		if(!intersect) intersect=SetAlloc(MAX_NUM_READS); else SetReset(intersect);
		SetIntersect(intersect, G->haps[h].readSet, _site[i].readsThatTouch);
		printf("\thap[%d] touches %d reads, is majority %d, with MEC %d\n", h, SetCardinality(intersect),
		    G->haps[h].sol[i], G->haps[h].MEC[i]);
	    }
	}
    }
    printf(" total genomeMEC %d\n", genomeMEC);
}

void FlipRead(READ *r) {

}

#define MAX_TRIES 1000
void HillClimb(GENOME *G) {
    int maxMEC=0, iter=0, numTries=MAX_TRIES;
    // Note: the "stagnant" variable isn't needed in SA, it's only needed in Hill Climbing, because if we've
    // tried 1000 moves without any improvement, it's probably time to give up because we're at a local minimum.
    while(numTries) { // for SA, this "while" will be replace with a loop over the temperature range estimated first by SA
	if(iter % 1000 == 0) Report(iter, G);
	++iter;

	int r=drand48()*MAX_NUM_READS; // pick a read at random
	foint f; f.i = r; // pass this foint f into the "move" function below. It will look like this:
	// foint Move(foint f) { int r = f.i; // this is the read to move....
	// compute score (before move)
	int beforeMEC=0, afterMEC=0;
	for(int i=0; i<_read[r].numSites;i++) beforeMEC += ComputeSiteMEC(G, &_site[_read[r].firstSite+i]);

	// compute new location
	int hap=_read[r].hap, hapShift = drand48()*(PLOIDY-1)+1, newHap = (hap+hapShift)%PLOIDY;
	assert(0 <= newHap && newHap < PLOIDY && newHap != hap);
	if(VERBOSE>1) printf(" iter %d Read %d is current in H%d, about to move it to H%d...", iter, r, hap, newHap);
	// Now ACTUALLY do the move and recompute the MEC along its sites (assumes PLOIDY==1)
	SetDelete(G->haps[hap].readSet, r);
	SetAdd  (G->haps[newHap].readSet, r);
	_read[r].hap=newHap;

	// compute score (after move--exact same calculation as "before" above)
	for(int i=0; i<_read[r].numSites;i++) afterMEC += ComputeSiteMEC(G, &_site[_read[r].firstSite+i]);

	// SA code will make this decision, NOT you
	if(beforeMEC>maxMEC || afterMEC>maxMEC) maxMEC=MAX(beforeMEC,afterMEC);
	if(VERBOSE>1) printf("before %d, after %d...", beforeMEC, afterMEC);
	// Here is where SA might tell you to accept or reject: your AcceptReject function will look at the first argument
	// (a Boolean) and simply do the accept, or reject, as instructed.
	if(afterMEC < beforeMEC) { // do this if told to accept
	    numTries=MAX_TRIES;
	    if(VERBOSE>1) printf("accept ") ; // accept the move (do nothing, it's already moved)
	}
	else { // do this if told to reject
	    numTries--;
	    if(VERBOSE>1) printf("reject(%d) ", numTries);
	    afterMEC=0;
	    _read[r].hap=hap;
	    SetDelete(G->haps[!hap].readSet, r);
	    SetAdd    (G->haps[hap].readSet, r);
	    for(int i=0; i<_read[r].numSites;i++) afterMEC += ComputeSiteMEC(G, &_site[_read[r].firstSite+i]);
	    assert(afterMEC == beforeMEC);
	}
    }
    printf("\nHill Climbing stagnated\n");
}

int main(int argc, char *argv[])
{
    srand48(time(NULL)+getpid());

    if(argc>1) ReadWIF(argv[1]);
    else CreateRandomReads(); // reads "created" from the true (but unknown here) haplotype set

    InitializeSystem();

    GENOME G;

    for(int h=0; h<PLOIDY; h++) {
	G.haps[h].id=h;
	G.haps[h].totalMEC=0;
	G.haps[h].readSet = SetAlloc(MAX_NUM_READS);
	G.haps[h].sol = Calloc(MAX_NUM_SITES,sizeof(G.haps[h].sol[0])); // computed from the most common element
	G.haps[h].MEC = Calloc(MAX_NUM_SITES,sizeof(G.haps[h].MEC[0])); // number of differences from most common element
    }

    // Initially, assign each read to a random haplotype
    for(int r=0; r<MAX_NUM_READS; r++) {
	int hap = _read[r].hap = PLOIDY*drand48();
	SetAdd(G.haps[hap].readSet, r);
    }

    printf("Global numSites %d, coverage %d, %d reads of length %d, numHap %d\n", MAX_NUM_SITES, COVERAGE, MAX_NUM_READS, MAX_READ_LEN, PLOIDY);
    HillClimb(&G);
}
