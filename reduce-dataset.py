from __future__ import print_function

import sys
import math
import numpy as np

if len(sys.argv) != 6:
    print("Usage: {} <target SNPs> <wif file> <ground truth file> <coverage> <read error %>".format(sys.argv[0]))
    sys.exit(1)

target_snps = int(sys.argv[1])
wif_file = sys.argv[2]
gt_file = sys.argv[3]
target_coverage = int(sys.argv[4])
r_error_rate = float(sys.argv[5])

wif = open(wif_file, 'r')
gt = open(gt_file, 'r')
sites = set()
letters = {}

for line in wif:
    sil = line.split('#')[0].strip().split(':')
    for s in sil:
        q = s.strip().split(' ')
        if len(q) != 4: continue

        sites.add(int(q[0]))
        
        if letters.get(int(q[0]), -1) == -1:
              letters[int(q[0])] = set()
        letters[int(q[0])].add(q[1])

for k in letters.keys():
    letters[k] = list(letters[k])
    if len(letters[k]) < 2:
        letters[k].append('X')

sites = list(sorted(sites))

startidx = np.random.randint(0, len(sites) - target_snps)
chosen = sites[startidx:startidx + target_snps]

haplotypes = []

for line in gt:
    haplotypes.append(line)

haplotypes[0] = haplotypes[0][startidx:startidx + target_snps]
haplotypes[1] = haplotypes[1][startidx:startidx + target_snps]

cur_coverage = 0
min_snp = target_snps
max_snp = 0
reads = {}
while cur_coverage < target_coverage:
    K = np.random.poisson(50)
    L = np.random.randint(target_snps)
    while L+K/2 >= target_snps or L-K/2 < 0:
        L = np.random.randint(target_snps)
    
    if reads.get(L-K/2, -1) == -1:
        reads[L-K/2] = []
    reads[L-K/2].append(L+K/2)

    cur_coverage += float(K) / target_snps
    min_snp = min(min_snp, L-K/2)
    max_snp = max(max_snp, L+K/2)

ks = list(reads.keys())
ks.sort()
for start in ks:
    for ends in reads[start]:
        H = np.random.randint(2)
        for i in range(start,ends + 1):
            print("{} X {} 61".format(i, haplotypes[H][i]), end=" : ")
        print("# 60 : NA")

# 4492 C 1 61 : 14636 T 0 61 :

#GROUND TRUTH
for h in haplotypes:
    print(h[min_snp:max_snp+1], file=sys.stderr)
