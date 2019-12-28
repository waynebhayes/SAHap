import random
import sys

def makeSeqFromFile(filename):
    out = []
    with open(filename) as f:
        for line in f:
            out.append([l for l in line.strip()])
    return out

'''
makeReadIndicies makes a list of 2-tuples, showing the beginning and ending of each read.
This function ensures that the reads are contingent.
'''
def makeReadIndices(seq, numReads, endread, readgap):
    '''
    endread should be a small number (<5), which is the number of reads on either side 
    readgap should be 2-tuples, with the max and min gap sizes (INCLUDING the max!)
    '''
    minReadLen = (readgap[0]+(2*endread))
    maxOverlap = int(((minReadLen*numReads)-len(seq[0]))/numReads)
    assert(int(maxOverlap)>0),"The minimum read is too small, and may fail to form contiguous blocks: try increasing the minimum readgap."
    readIndices = []
    reads=[]
    nextRead=(0,random.randint(readgap[0],readgap[1])+(2*endread))
    readIndices.append(nextRead)
    readsMade = 0
    while readsMade<numReads:
        overlap = random.randint(1,maxOverlap)
        nextReadStart=nextRead[1]-overlap
        nextReadEnd=nextReadStart+random.randint(readgap[0],readgap[1])+(2*endread)
        if(nextReadEnd>=len(seq[0])-1):
            nextReadEnd = len(seq[0])-1
            nextReadStart = nextReadEnd-(random.randint(readgap[0],readgap[1])+(2*endread))
            nextRead = (nextReadStart,nextReadEnd)
            readIndices.append(nextRead)
            readsMade+=1
            break
        nextRead = (nextReadStart,nextReadEnd)
        readIndices.append(nextRead)
        readsMade+=1
    readsLeftToMake = numReads - readsMade
    for i in range(0,readsLeftToMake):
        readSize = random.randint(readgap[0], readgap[1]) + (2*endread)
        nextReadStart = random.randrange(0,len(seq[0])-readSize)
        '''
        error above?
        '''
        nextReadEnd = nextReadStart+readSize
        nextRead = (nextReadStart,nextReadEnd)
        readIndices.append(nextRead)
    return readIndices

'''
makeReads takes the indices and the error rate, and produces the actual error.
'''
def makeReads(readIndicies, seq, endread, err):
    letters = ["a","t","g","c"]
    out = []
    for startI, endI in readIndicies:
        k = random.randrange(0,len(seq))
        frontEndRead = seq[k][startI:startI+endread]
        backEndRead = seq[k][endI-endread+1:endI+1]
        for i in range(0,len(frontEndRead)):
            chance = random.uniform(0,1)
            if chance<err:
                frontEndRead[i] = random.choice([l for l in letters if l!=frontEndRead[i]])
        for i in range(0,len(backEndRead)):
            chance = random.uniform(0,1)
            if chance<err:
                backEndRead[i] = random.choice([l for l in letters if l!=backEndRead[i]])
        out.append((startI, endI-endread, frontEndRead, backEndRead))
    return out

if __name__ == "__main__":
    script,filename,numReads,endread,gapmin,gapmax,err = sys.argv
    seq1 = makeSeqFromFile(filename)
    readsI1 = makeReadIndices(seq1, int(numReads), int(endread), (int(gapmin),int(gapmax)))
    reads1 = makeReads(readsI1, seq1, int(endread), float(err))
    for num1, num2, leftRead,rightRead in reads1:
        print(num1, end = "\t")
        for i in leftRead:
            print(i,end="")
        print("\t",end="")
        print(num2, end = "\t")
        for i in rightRead:
            print(i,end="")
        print("")
        
        
        
        
        
    