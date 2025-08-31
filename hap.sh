#!/bin/bash
################## SKELETON: DO NOT TOUCH THESE 2 LINES
EXEDIR=`dirname "$0"`; BASENAME=`basename "$0" .sh`; TAB='	'; NL='
'
#################### ADD YOUR USAGE MESSAGE HERE, and the rest of your code after END OF SKELETON ##################
USAGE="USAGE: $BASENAME [-c] foo.wif
PURPOSE: script to perform... stuff... on a haplotype WIF file
OPTIONS:
    -v: use the site variant (3rd entry of the 4-entry site specifier) rather than the genomic letter (default)"

################## SKELETON: DO NOT TOUCH CODE HERE
# check that you really did add a usage message above
USAGE=${USAGE:?"$0 should have a USAGE message before sourcing skel.sh"}
die(){ echo "$USAGE${NL}FATAL ERROR in $BASENAME:" "$@" >&2; exit 1; }
[ "$BASENAME" == skel ] && die "$0 is a skeleton Bourne Shell script; your scripts should source it, not run it"
echo "$BASENAME" | grep "[ $TAB]" && die "Shell script names really REALLY shouldn't contain spaces or tabs"
[ $BASENAME == "$BASENAME" ] || die "something weird with filename in '$BASENAME'"
warn(){ (echo "WARNING: $@")>&2; }
not(){ if eval "$@"; then return 1; else return 0; fi; }
newlines(){ awk '{for(i=1; i<=NF;i++)print $i}' "$@"; }
parse(){ awk "BEGIN{print $*}" </dev/null; }
which(){ echo "$PATH" | tr : "$NL" | awk '!seen[$0]{print}{++seen[$0]}' | while read d; do eval /bin/ls $d/$N; done 2>/dev/null | newlines; }

unset TMPDIR

[ "$MYTMP" ] || export MYTMP="/tmp"
export TMPDIR=${TMPDIR:-`mktemp -d $MYTMP/$BASENAME.XXXXXX`}
 trap "/bin/rm -rf $TMPDIR; exit" 0 1 2 3 15 # call trap "" N to remove the trap for signal N

#################### END OF SKELETON, ADD YOUR CODE BELOW THIS LINE

# General idea: for each site, list the letters recorded from all reads that touch that site, sorted by site+letter.
# This will place agreeing reads at that site all in one place (kinda like k-mers).
# Then gather them by line-by-line as (site, letter, {set of agreeing reads}) (called slRs)
# Then for each slR, record all the read pairs that agree with each other.
# Once finishing reading all the slRs, we have an "agree count" for every pair of reads.
# Then sort all the read pairs by agree count, and greedily build the two haplotype sets.

LETTER_COL=1 # set this to 2 if we only get the digit (variation) of the site rather than the actual genomic letter
case "$1" in
-v) LETTER_COL=2; shift ;;
-*) die "unknown option '$1'";;
esac

[ $# -eq 1 ] || die "expecting exactly 1 argument, a WIF file"
[ -r "$1" ] || die "1st arg must be a file"

# First, separate the WIF file into scaffolds, because there is no way to disambiguate haplotypes between scaffolds
hawk -F: 'BEGIN{scaff=0; fmt="'$TMPDIR'/scaff%06d.wif"; outFile=sprintf(fmt,scaff)}
    {sub(" : #.*","")} # nuke the useless colon and comment
    NR==1{furthest=1*$NF}
    NR>1{if(1*$1 > furthest) # the first site of this read is beyond the reach of all reads so far
	    outFile=sprintf(fmt,++scaff);
    }
    {
	print > outFile; reach=1*$NF; if(reach>furthest)furthest=reach
    }' "$@"

cd $TMPDIR || die "couldn't cd into '$TMPDIR'"
echo There are `ls | wc -l` scaffolds

for scaff in scaff*.wif;
do
    echo "Scaffold $scaff";
    # a WIF file is one read per line
    sed 's/ : #.*$//' "$scaff" | # remove the final colon and comment text
	awk '{for(c=1;c<NF;c+=5)print $c,$(c+'$LETTER_COL'),FNR}' | # site, letter, and read, in the order seen in the WIF file
	sort | tee $TMPDIR/slr.txt | # now sorted by site as a STRING--worry about sorting them numerically later
	gawk '{ # prevSL = prev site+letter
	    if(prevSL==$1" "$2) printf " %d", $3; # print read number of the site if site+letter agree with previous line
	    else {
		sameSL=0; if(start++) print ""; # print newline except first time
		prevSL=$1" "$2; printf "%s %s\t%s",$1,$2,$3 # otherwise print new site, letter, and read
	    }
	}' | sort -n | tee $TMPDIR/slRs.txt | # sorted by "site letter {set of reads that have this letter at this site}"
	hawk 'function SetLine(res,   i){delete res; for(i=3;i<=NF;i++)res[$i]=1}
	function PrintGroups(G,  g,r,p) {
	    p=PROCINFO["sorted_in"];
	    PROCINFO["sorted_in"]="@ind_num_asc";
	    for(g in G){
	    printf "G[%s]:", g;
		for(r in G[g]) printf " %s", r;
		print ""
	    }
	    PROCINFO["sorted_in"]=p;
	}
	{printf "%s", $0} # print every line
	NR==1{ # first line is special: initialize H[1] (entire Haplotype 1)
	    numGroups=1; H[1][0]=1;
	    SetLine(H[1]);
	    print ""
	}
	NR>1{
	    SetLine(L); # L becomes the set of reads that see this letter at this site
	    MakeEmptySet(numMatches);
	    for(g in H) {
		# number of matches = overlap between the reads seeing this letter at this site, and the haplotypes-so-far
		nm = SetIntersect(res,H[g],L);
		if(nm) numMatches[g] = nm;
	    }
	    nm=0; delete bm;
	    PROCINFO["sorted_in"]="@val_num_desc"; # force the next "for" loop to iterate through numMatches largest-to-smallest
	    for(g in numMatches) bm[++nm]=g; # bm[1] is the haplotype with the greatest number of reads that match this site.
	    ASSERT(length(numMatches) == length(bm));
	    if(length(numMatches)==0) { # this must be a new group
		print "\tNo match to existing groups:"; PrintGroups(H);
		printf "New group %d\n", ++numGroups;
		H[numGroups][0]=1;
		SetCopy(H[numGroups],L);
	    } else if(length(numMatches)==1) { # this line is a continuation of exactly one existing group
		SetUnion(res, H[bm[1]], L); # extend the Haplotype
		SetCopy(H[bm[1]], res);
		print "\tMatches only group", bm[1]
	    } else {
		printf "\tDo nothing because there is more than one match:";
		for(i=1;i<=length(bm);i++) printf " %d has %d;",bm[i],numMatches[bm[i]];
		print ""
		PrintGroups(H);

        ambigCount++
        ambigReads[$3] = 1
        }
	}
	END{
	    print "FINAL GROUPS"
        PrintGroups(H);
        print ""
        print "AMBIGUOUS READS"

        PROCINFO["sorted_in"]="@ind_num_asc";

        count = 0
        for (read in ambigReads) {
            count++
            printf "Ambiguous[%d]: %s\n", count, read
        }

        print count " total ambiguous cases."
	}'
done

exit

 *************************** OLD, PROBABLY USELESS CODE BELOW THIS LINE **********************************
    hawk '{for(i=3;i<NF;i++) for(j=i+1;j<=NF;j++) ++agree[MIN($i,$j)][MAX($i,$j)]} # accumulate #sites at which reads agree
	END{for(i in agree)for(j in agree[i]) print agree[i][j],i,j}' |
    sort -nr | tee $TMPDIR/nrr.txt | # countAgree r1 r2
# Now go through the nrr file twice--once to get the length of each read, and second time to do merging.
hawk 'BEGIN{MakeEmptySet(G)}
    ARGIND==1{r1=$2;r2=$3; ++len[r1]; ++len[r2]}
    ARGIND==2{
	printf "* %s *", $0 > "/dev/stderr";
	r1=$2;r2=$3; ++len[r1]; ++len[r2];
	if(!(r1 in G)&&!(r2 in G)) { # neither are in a group
	    G[r1]=G[r2]=++ng; A[ng]=$1; set[ng][r1]=set[ng][r2]=1; # A[group]=totalAgreeCount
	    printf "\tboth assigned to group %d [%s %s]\n",ng, G[r1],G[r2] > "/dev/stderr"
	}
	else if((r1 in G)&&(r2 in G)) { # both are in a group
	    gr1=G[r1]; gr2=G[r2];
	    printf "G[%d]=%d |%d|, G[%d]=%d |%d|;",r1,gr1,A[gr1], r2,gr2, A[gr2] > "/dev/stderr"
	    ASSERT(set[gr1][r1] && set[gr2][r2], set[gr1][r1]" && "set[gr2][r2]);
	    if(G[r1]==G[r2]) print " all good" > "/dev/stderr"
	    else if($1 < (len[r1]+len[r2])/'$HAP') print "not enough agreement to merge" > "/dev/stderr"
	    else { # merge the groups *IF* they agree more than half their length
		# merge the groups into the smaller numbered group
		src=MAX(gr1,gr2); dest=MIN(gr1,gr2);
		printf "merging set %d into set %d\n", src, dest > "/dev/stderr"
		A[dest]+=A[src];
		for(r in set[src]) {set[dest][r]=1; G[r]=dest}
		delete set[src]; delete A[src];
	    }
	} # at this point EXACTLY one is in a group
	else if(r1 in G){ASSERT(!(r2 in G),r2" oops "G[r2]); gr1=G[r1]; G[r2]=gr1; A[gr1]+=$1; set[gr1][r2]=1; printf "G[%d]=%d\n",r2,gr1 > "/dev/stderr"}
	else if(r2 in G){ASSERT(!(r1 in G),r1" oops "G[r1]); gr2=G[r2]; G[r1]=gr2; A[gr2]+=$1; set[gr2][r1]=1; printf "G[%d]=%d\n",r1,gr2 > "/dev/stderr"}
	else ASSERT(0, "should not get here");
    }
    END{ printf "\nFinal groupings:\n";
	for(g in set) {
	    gFile = sprintf("'$TMPDIR'/g%d",g);
	    printf "%d reads, set %s =\t{", length(set[g]), g;
	    for(r in set[g]){printf " %s",r; printf " %d$\n", r > gFile}
	    print " }"
	}
    }' $TMPDIR/nrr.txt $TMPDIR/nrr.txt | sort -nr

ls -S $TMPDIR/g* | while read g; do
    ls -l $g | awk '{printf "%d ", $5}'; ggrep -f $g $TMPDIR/slr.txt | awk '{print $2}' | uniq -c | awk '$1>1{printf "%s",$2}END{print "\n"}'
done
