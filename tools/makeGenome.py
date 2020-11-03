import random
import sys

'''
makeHaplotypeSequence only creates the base truth.
Example input: (2, 5)
Example output: [[atgtt],[tgacc]]
'''
def makeHalpotypeSequence(ploidy, haplotypeLength):
    letters = ["a","t","g","c"]
    genomes = [[] for _ in range(0, ploidy)]
    for _ in range(0, haplotypeLength):
        toAdd=[]
        allEqual=True
        toAdd.append(random.choice(letters))
        for i in range(1,ploidy):
            toAdd.append(random.choice(letters))
            allEqual = (toAdd[i]==toAdd[i-1])
        if allEqual:
            toAdd[random.randrange(ploidy-1)]=random.choice([i for i in letters if i!=toAdd[ploidy-2]])
        for i in range(0,ploidy):
            genomes[i].append(toAdd[i])
    return genomes

if __name__ == "__main__":
    script, ploidy, haplotypeLength = sys.argv
    seq1 = makeHalpotypeSequence(int(ploidy), int(haplotypeLength))
    for seq in seq1:
        for i in seq:
            print(i,end="")
        print("")