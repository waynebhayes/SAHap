#ifndef __SAHAP_H__

#define K 2 // ploidy (number of chromosomes)
extern int m; // number of SNP positions
extern char *S[K]; // current best solution
typedef struct _readPlace READ_PLACEMENT;

typedef struct _ReadPair {
    // CONSTANTS: these don't change once the read is created.
    int gap; // estimated gap between the two ends (may not be used)
    int pos[2]; // the starting positions [0,..m) of the two reads; assert(pos[1] > pos[0]+len[0]);
    int len[2]; // lengths of the two reads
    char *seq[2]; // the actual sequences, of length len[0] and len[1]
    // Above are all constants. Below is our one variable.
    // GOAL VARIABLE: our goal is to figure out the correct value of "which", meaning which chromosome does this come from?
    char which; // which chromosome are we currently associated with?? 0 through K-1
} READ_PAIR;

