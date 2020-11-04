#!/bin/bash
USAGE=' '
source ~/bin/misc.sh
REG_DIR=${REG_DIR:?"needs to be set"}

PARALLEL="awk '{printf "'"(%s)\n",$0}'"' | bash" # in case there's no parallel
if which parallel >/dev/null; then
    CORES=${CORES:?"needs to be set to the number of cores (cpus) you want to use in parallel"}
    echo | awk 'BEGIN{for(i=0;i<2;i++)printf "hello %d\n",i}' > $TMPDIR/correct
    echo | awk 'BEGIN{for(i=0;i<2;i++)printf "sleep %d; echo hello %d\n",i,i}' | parallel 2 > $TMPDIR/parallel
    if cmp $TMPDIR/*; then
	PARALLEL="parallel $CORES"
    fi
fi
if echo "$PARALLEL" | fgrep -q bash; then CORES=1; fi

TRIES=0
NEED='MEC Poisson'
echo "Running tests in parallel with $CORES cores"
while [ "$NEED" != "" -a $TRIES -lt 10 ]; do
    echo "Before try #$TRIES, need '$NEED'"
    for OBJ in $NEED; do
	echo "./sahap.$OBJ data/500SNPs_30x/Model_14.wif data/500SNPs_30x/Model_14.txt 1 > $REG_DIR/$OBJ.out"
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
exit $NUM_FAILS
