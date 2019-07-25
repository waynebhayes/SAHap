#ifndef __SAHAP_H__

extern int Slen;
extern char *S; // current best solution
typedef struct _readPlace READ_PLACEMENT;

typedef struct _ReadPair {
    char *seq[2]; // two sequences of characters, ie., ACGT's
    int gap, len[2]; // estimated gap between, and precise lengths of, the two sequences
    READ_PLACEMENT *currentPlace;
} READ_PAIR;

typedef struct _readPlace {
    READ_PAIR *R; // the read we are referring to

    // Orientation: it must be +1 or -1, meaning the read is placed forward or backward along the sequence
    int orientation;

    unsigned int loc[2]; // loc is always positive. However...
    // if orientation is +1:
	// loc[0] is the left-most location in S[] of R->seq[0][0], and
	// loc[1] is the left-most location in S[] of R->seq[1][0].
    // if orientation is -1:
	// R->seq[0][i] is at S[Slen-1-loc[0]-i].
	// R->seq[1][i] is at S[Slen-1-loc[1]-i].
    // Note also that the objective function should probably have some element of preferring that
    // ABS(loc[1]-loc[0]) is close to R->gap, so for example penalize by ABS(ABS(loc[1]-loc[0])-gap),
    // although that may be too much.
} READ_PLACEMENT;

#endif

