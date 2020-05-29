import sys
import random

if len(sys.argv) < 4:
    print("Usage: {} <target SNPs> <wif file> <ground truth file>".format(sys.argv[0]))
    sys.exit(1)

target = int(sys.argv[1])
wif_file = sys.argv[2]
gt_file = sys.argv[3]

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

startidx = random.randint(0, len(sites) - target)
chosen = sites[startidx:startidx + target]

wif.seek(0)
for line in wif:
    sil = line.split('#')[0].strip().split(':')
    l = []
    for s in sil:
        q = s.strip().split(' ')
        if len(q) != 4: continue

        if int(q[0]) in chosen:
            l.append(s)

    if l:
        print(' : '.join(l))

for line in gt:
    print(line.strip()[startidx:startidx+target], file=sys.stderr)
