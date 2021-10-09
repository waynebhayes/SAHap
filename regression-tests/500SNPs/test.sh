#!/bin/bash
USAGE=' '
REG_DIR=${REG_DIR:?"needs to be set"}
SAHAP_ITERS=${SAHAP_ITERS:=1}

# Temporary Filename + Directory (both, you can use either, note they'll have different random stuff in the XXXXXX part)
TMP=`mktemp /tmp/$BASENAME.XXXXXX`
TMPDIR=`mktemp -d /tmp/$BASENAME.XXXXXX`
trap "/bin/rm -rf $TMP $TMPDIR; exit" 0 1 2 3 15 # call trap "" N to remove the trap for signal N

warn(){ echo "WARNING: $@" >&2;}
newlines(){ awk '{for(i=1; i<=NF;i++)print $i}' "$@";}

PARALLEL="awk '{printf "'"(%s)\n",$0}'"' | bash" # in case there's no parallel
if [ -x ./parallel ]; then
    CORES=${CORES:?"needs to be set to the number of cores (cpus) you want to use in parallel"}
    echo | awk 'BEGIN{for(i=0;i<2;i++)printf "hello %d\n",i}' > $TMPDIR/correct
    echo | awk 'BEGIN{for(i=0;i<2;i++)printf "sleep %d; echo hello %d\n",i,i}' | ./parallel 2 > $TMPDIR/parallel
    if cmp $TMPDIR/*; then
	PARALLEL="./parallel $CORES"
    else
	warn "./parallel 2 gave wrong output:"
	set -x; diff $TMPDIR/*; set +x
    fi
fi
if echo "$PARALLEL" | fgrep -q bash; then
    warn "No parallel available; using only 1 core"
    CORES=1
fi

TRIES=0
NEED='MEC Poisson'
echo "Running tests in parallel with $CORES cores"
while [ "$NEED" != "" -a $TRIES -lt 10 ]; do
    echo "Before try #$TRIES, still need success on '$NEED'"
    for OBJ in $NEED; do
	echo "time ./sahap.$OBJ data/500SNPs_30x/Model_14.wif data/500SNPs_30x/Model_14.txt $SAHAP_ITERS > $REG_DIR/$OBJ.out 2>$REG_DIR/$OBJ.err; cat $REG_DIR/$OBJ.err | tee -a $REG_DIR/$OBJ.out"
    done | eval $PARALLEL
    (( NUM_FAILS+=$? ))

    NUM_FAILS=0
    for OBJ in $NEED; do
	echo --- $REG_DIR/$OBJ.out ---
	fgrep '(100.0' $REG_DIR/$OBJ.out | fgrep '%' |
	    awk 'function MIN(a,b){if(a<=b)return a; else return b;}
	    BEGIN{
		# Best possible, so far as I know:
		# 10k (0.0%,1s)  T 0.259  fA 0.000  pBad 0.0000  MEC 4.88  ( Err_vs_truth    13 Err_Pct 1.34% [486 486])
		best["MEC"]=4.27;
		best["Err_vs_truth"]=4;
		best["Err_Pct"]=0.41;
	    }
	    /\(100\.0*%/{
		#print;
		for(i=1;i<NF;i++)if($i in best){
		    got=$(i+1);
		    if(!($i in GOT))GOT[$i]=got;
		    else GOT[$i]=MIN(GOT[$i],got);
		}
	    }
	    END {
		for(m in best) {
		    printf "%s best %g got \"%s\": ", m,best[m],GOT[m];
		    if(!(m in GOT)) print "FAIL: NO RESULT"
		    else if(GOT[m]<=1*best[m]) {++NUM_GOOD; print "BEST!!";}
		    else if(GOT[m]<=2*best[m]) {++NUM_GOOD; print "GOOD";}
		    else if(GOT[m]<=3*best[m]) {++NUM_GOOD; print "FAIR";}
		    else if(GOT[m]<=4*best[m]) {++NUM_GOOD; print "BORDERLINE";}
		    else print "FAIL"
		}
		exit(length(best)-NUM_GOOD);
	    }'
	NEW_FAILS=$?
	if [ $NEW_FAILS -eq 0 ]; then
	    NEED=`echo $NEED | newlines | grep -v $OBJ`
	else
	    (( NUM_FAILS+=$NEW_FAILS ))
	fi
    done
    echo "After try #$TRIES, NUM_FAILS is $NUM_FAILS, still looking for success on '$NEED'"
    (( ++TRIES ))
done
echo "After $TRIES tries, NUM_FAILS is $NUM_FAILS; still looking for success on '$NEED'"
exit $NUM_FAILS
