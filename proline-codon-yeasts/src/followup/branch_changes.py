"""Acceptor-stem/N73 changes on the target branch, per tRNA class: target genomes vs the sisters' majority.
usage: python branch_changes.py trnas.tsv. Also counts C1:G72 among UAG tRNAs of all 49 genomes."""
import sys, csv
from collections import defaultdict, Counter
sys.path.insert(0, 'src/followup')  # run from the folder root
from origin import SIS, WICK, SPIN, CUN, ac

by = defaultdict(list)
for r in csv.DictReader(open(sys.argv[1]), delimiter='\t'):
    by[r['acc']].append(r)


def key(r):
    return 'UAG' if ac(r) in CUN else ac(r)


def acc8(r):
    return r['stem5'] + r['stem3'] + r['N73']


for grp in (WICK, SPIN):
    sis = SIS[grp[0]]
    print('== targets', grp)
    n_cls = n_chg = 0
    for k in sorted({key(r) for a in grp for r in by[a]}):
        ss = [acc8(r) for a in sis for r in by[a] if key(r) == k]
        ts = [acc8(r) for a in grp for r in by[a] if key(r) == k]
        if not ss or not ts or any(len(x) != 15 for x in ss + ts):
            continue
        cons = ''.join(Counter(x[i] for x in ss).most_common(1)[0][0] for i in range(15))
        # positions changed in every target copy vs sisters' consensus
        lab = [f'{i + 1}' for i in range(7)] + [f'{66 + i}' for i in range(7)] + ['73']
        chg = [f'{lab[i]}:{cons[i]}>{ts[0][i]}' for i in range(15) if all(t[i] != cons[i] for t in ts)
               and all(t[i] == ts[0][i] for t in ts)]
        n_cls += 1
        n_chg += bool(chg)
        if chg:
            print(f'  {k:4s} changed in all targets: {chg}')
    print(f'  classes compared {n_cls}, classes with a shared change {n_chg}')
accs = open('data/tree/accs.txt').read().split()
c = Counter((r['stem5'][0] + ':' + r['stem3'][-1]) for a in accs for r in by[a] if ac(r) in CUN)
print('1:72 of UAG tRNAs, 49 genomes:', dict(c))
print('genomes with a C1:G72 UAG tRNA:', [a for a in accs for r in by[a] if ac(r) in CUN and r['stem5'][0] + r['stem3'][-1] == 'CG'])
