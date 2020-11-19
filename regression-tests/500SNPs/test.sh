#!/bin/bash
USAGE=' '
REG_DIR=${REG_DIR:?"needs to be set"}

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
    echo "Before try #$TRIES, need '$NEED'"
    for OBJ in $NEED; do
	echo "/bin/time ./sahap.$OBJ data/500SNPs_30x/Model_14.wif data/500SNPs_30x/Model_14.txt 1 > $REG_DIR/$OBJ.out 2>&1"
    done | eval $PARALLEL

    NUM_FAILS=0
    for OBJ in $NEED; do
	echo --- $REG_DIR/$OBJ.out ---
	fgrep '(100.0000%' $REG_DIR/$OBJ.out |
	    awk 'BEGIN{
		best["MEC"]=71;
		best["Err_vs_truth"]=4;
		best["Err_Pct"]=0.411523
	    }
	    /\(100\.0*%/{
		for(i=1;i<NF;i++)if($i in best){
		    got=$(i+1); printf "%s best %d got %s: ", $i,best[$i],got;
		    if(got>2*best[$i]){
			print "FAIL";
			++NUM_FAILS;
		    }
		    else print "GOOD"
		}
	    }
	    END {
		exit(NUM_FAILS);
	    }'
	NEW_FAILS=$?
	if [ $NEW_FAILS -eq 0 ]; then
	    NEED=`echo $NEED | newlines | grep -v $OBJ`
	else
	    (( NUM_FAILS+=$NEW_FAILS ))
	fi
    done
    echo "After try #$TRIES, NUM_FAILS is $NUM_FAILS, still need '$NEED'"
    (( ++TRIES ))
done
echo "Needed a total of $TRIES attempts to get all regressions to pass"
[ "$TRIES" -gt 1 ] && warn "Needed more than one try??"
exit $NUM_FAILS
