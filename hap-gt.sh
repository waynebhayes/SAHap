#!/bin/bash
################## SKELETON: DO NOT TOUCH THESE 2 LINES
EXEDIR=`dirname "$0"`; BASENAME=`basename "$0" .sh`; TAB='	'; NL='
'
#################### ADD YOUR USAGE MESSAGE HERE, and the rest of your code after END OF SKELETON ##################
USAGE="USAGE: $BASENAME <base-filename>
PURPOSE: given the basename of a WIF file and the associated ground truth .TXT file, evaluate the GT"

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

[ "$MYTMP" ] || export MYTMP="/tmp"
export TMPDIR=${TMPDIR:-`mktemp -d $MYTMP/$BASENAME.XXXXXX`}
 trap "/bin/rm -rf $TMPDIR; exit" 0 1 2 3 15 # call trap "" N to remove the trap for signal N

#################### END OF SKELETON, ADD YOUR CODE BELOW THIS LINE

F="$1"
[ -f $F.txt ] || die "can't find ground truth file $F.txt"
[ -f $F.wif ] || die "can't find WIF file $F.wif"

sed 's/ : #.*//' $F.wif | tr : "$NL" | awk '{print $1}' | sort -u | sort -n > $TMPDIR/sites

echo "Ambiguous reads are those in which the winning haplotpype have fewer than twice the matches of the losing one, according the ground truth in the .txt file. Here they are sorted by the line number they appeared in the WIF file"
hawk 'ARGIND==1{for(i=1;i<=length($0);i++)L[FNR][i]=substr($0,i,1);next}
    ARGIND==2{siteline[$1]=FNR;next}
    {
	sub(" : #.*","");
	delete m;
	for(i=0;i<NF/5;i++) for(k=1;k<=2;k++) if($(5*i+3)==L[k][siteline[$(5*i+1)]])++m[k];
	printf "%3d %3d\t%d\n", MIN(m[1],m[2]), MAX(m[1],m[2]), FNR
    }' $F.txt $TMPDIR/sites $F.wif |
    sort -n |
    awk '$2<2*$1' | # find reads where the winner matches fewer than twice the loser
    sort -k 3n | # print the ambiguous reads sorted by the line number of the WIF file they appeared on
    cut -f2
