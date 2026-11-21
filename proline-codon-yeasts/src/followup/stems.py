"""Side-by-side acceptor stem, N73, anticodon loop and variable-region length of UAG and Pro tRNAs, targets vs sisters.
usage: python stems.py trnas.tsv. V = nt between anticodon-stem 3' strand end and T-stem start (T arm taken as 17 nt)."""
import sys, csv
from collections import defaultdict
sys.path.insert(0, 'src/followup')  # run from the folder root
from origin import SIS, WICK, SPIN, CUN, ac

by = defaultdict(list)
for r in csv.DictReader(open(sys.argv[1]), delimiter='\t'):
    by[r['acc']].append(r)


def row(acc, r, tag):
    s, o, n = r['seq'].upper(), int(r['ac_off']) - 1, int(r['ac_len'])
    i3 = s.rfind(r['stem3'])
    v = i3 - 17 - (o + n + 7)
    loop = s[o - 2:o + n + 2]
    print(f"{tag:6s} {acc} {ac(r):5s} len {len(s):3d}  5'{r['stem5']} 3'{r['stem3']} 1:72 {r['stem5'][0]}:{r['stem3'][-1]} "
          f"N73 {r['N73']} loop {loop} V {v}")


for grp in (WICK, SPIN):
    for acc in grp + SIS[grp[0]]:
        for r in by[acc]:
            if ac(r) in CUN:
                row(acc, r, 'TARGET' if acc in grp else 'sister')
    for acc in grp + SIS[grp[0]]:
        for r in by[acc]:
            if ac(r) == 'TGG':
                row(acc, r, 'Pro')
    print()
