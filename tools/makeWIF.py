from __future__ import print_function

import sys
import numpy as np 

# Create WIF files from scratch

if len(sys.argv) != 5:
    print("Usage: {} <Ploidy Count> <Target SNPs> <Target Coverage> <Read Error Rate>".format(sys.argv[0]))
    sys.exit(1)

ploidy = int(sys.argv[1])
snp_count = int(sys.argv[2])
target_coverage =  int(sys.argv[3])
read_errors = float(sys.argv[4])

if read_errors >= 1:
    read_errors /= 100.0

haplotypes = [[] for i in range(ploidy)]

for i in range(snp_count):
    haplist = list(range(0, ploidy))
    np.random.shuffle(haplist)
    for h in range(0, ploidy):
        haplotypes[h].append(haplist[h])

#WIF FILE
cur_coverage = 0.0
min_snp = snp_count
max_snp = 0
reads = {}
while cur_coverage < target_coverage:
    K = np.random.poisson(50)
    L = np.random.randint(snp_count)
    while L+K/2 >= snp_count or L-K/2 < 0:
        L = np.random.randint(snp_count)

    if reads.get(L-K/2, -1) == -1:
        reads[L-K/2] = []
    reads[L-K/2].append(L+K/2)

    cur_coverage += float(K) / snp_count
    min_snp = min(min_snp, L-K/2)
    max_snp = max(max_snp, L+K/2)

ks = list(reads.keys())
ks.sort()


temp_wif = open('data/temp/temp.wif', 'w')

for start in ks:
    for end in reads[start]:
        H = np.random.randint(2)
        err_pos = -1
        if np.random.random_sample() <= read_errors:
            err_pos = np.random.randint(end + 1)
        for i in range(int(start),int(end) + 1):
            if i == err_pos:
                haplotypes[H][i] = int(not haplotypes[H][i])
            temp_wif.write(f'{i} X {haplotypes[H][i]} 61 : ')
            if i == err_pos:
                haplotypes[H][i] = int(not haplotypes[H][i])
        temp_wif.write("# 60 : NA")
        
temp_wif.close()


# 4492 C 1 61 : 14636 T 0 61 :
#GROUND TRUTH

tempGT_txt = open('data/temp/tempGT.txt', 'w')

for haplotype in haplotypes:
    for i in haplotype[int(min_snp):int(max_snp)+1]:
        tempGT_txt.write(str(i))
    tempGT_txt.write("\n")

tempGT_txt.close()

print("data files located in 'SAHap/data/temp/'")