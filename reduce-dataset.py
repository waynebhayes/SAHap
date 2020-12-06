import sys
import random
import math

if len(sys.argv) != 6:
    print("Usage: {} <target SNPs> <wif file> <ground truth file> <coverage %> <read error %>".format(sys.argv[0]))
    sys.exit(1)

target = int(sys.argv[1])
wif_file = sys.argv[2]
gt_file = sys.argv[3]
target_coverage = float(sys.argv[4])
r_error_rate = float(sys.argv[5])

if target_coverage > 1:
    target_coverage /= 100
if r_error_rate > 1:
    r_error_rate /= 100

wif = open(wif_file, 'r')
gt = open(gt_file, 'r')
sites = set()

for line in wif:
    sil = line.split('#')[0].strip().split(':')
    for s in sil:
        q = s.strip().split(' ')
        if len(q) != 4: continue

        sites.add(int(q[0]))

sites = list(sorted(sites))

site_map = {}
for i in range(len(sites)):
    site_map[sites[i]] = i

if (len(sites) < target):
    target = len(sites)
    #print("Target SNPs was greater than # of sites reducing target to current # of sites")

startidx = random.randint(0, len(sites) - target)
chosen = sites[startidx:startidx + target]

R = []

let_insites = {}
for c in chosen:
    let_insites[c] = {'A', 'C', 'T', 'G'}

wif.seek(0)
for line in wif:
    sil = line.split('#')[0].strip().split(':')
    l = []
    for s in sil:
        q = s.strip().split(' ')
        if len(q) != 4: continue

        if int(q[0]) in chosen:
            l.append(s.strip())
            let_insites[int(q[0])].discard(q[1])

    if l:
        R.append(' : '.join(l))

random.shuffle(R)

LETTERS = ['A', 'C', 'T', 'G']
NR = []
chosen = set()

for i in range(0, math.ceil(len(R) * target_coverage)):
    sil = R[i].strip().split(':')
    res = []
    l = []
    for s in sil:
        q = s.strip().split(' ')
        if len(q) != 4: continue

        if random.random() < r_error_rate:
            for let in let_insites[int(q[0])]:
                q[1] = let
                break;
            #ind = (LETTERS.index(q[1]) + 1) % 4
            #q[1] = LETTERS[ind]
        chosen.add(int(q[0]))
        l.append(' '.join(q))

    if l:
        print(' : '.join(l))#, ' # 60 : NA']))

chosen = list(sorted(chosen))

for line in gt:
    #print(line.strip()[startidx:startidx+target], file=sys.stderr)
    l = line.strip()
    for s in chosen:
        print(l[site_map[s]], end="", file=sys.stderr)
    print("", file=sys.stderr)

