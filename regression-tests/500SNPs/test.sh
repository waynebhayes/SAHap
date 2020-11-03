#!/bin/sh
REG_DIR=${REG_DIR:?"needs to be set"}

<<<<<<< HEAD
TRIES=0
NUM_FAILS=1
while [ $NUM_FAILS -ne 0 -a $TRIES -lt 10 ]; do
    for OBJ in MEC Poisson; do
	echo "./sahap.$OBJ data/500SNPs_30x/Model_14.wif data/500SNPs_30x/Model_14.txt 1 > $REG_DIR/$OBJ.out"
    done | tee /dev/tty | parallel $CORES
=======
for OBJ in MEC Poisson; do
    echo "./sahap.$OBJ data/500SNPs_30x/Model_14.wif data/500SNPs_30x/Model_14.txt 1 > $REG_DIR/$OBJ.out"
done | tee /dev/tty | parallel $CORES
>>>>>>> aeccd44aa969a7c45afca0e337169a81d14ce9ac

    NUM_FAILS=0
    for OBJ in MEC Poisson; do
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
	(( NUM_FAILS+=$NEW_FAILS ))
    done
    echo "Try #$TRIES, NUM_FAILS is $NUM_FAILS"
    (( ++TRIES ))
done
exit $NUM_FAILS
