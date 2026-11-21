"""T0b: candidate tRNA loci missed by ARAGORN (>=5 query hits, <30 nt overlap with a call): ARAGORN on locus +-50 nt
(linear), and the subject triplet aligned to the query anticodon (34-36) when the HSP spans it.
usage: python t0b.py trnas.tsv hits_aln.tsv targets.fa"""
import sys, csv
from collections import defaultdict
import pyaragorn

calls = list(csv.DictReader(open(sys.argv[1]), delimiter='\t'))
qinfo = {f"{r['acc']}|{r['anticodon']}|{i}": r for i, r in enumerate(
    [r for r in calls if r['acc'] in set(open('data/tree/accs.txt').read().split())])}
iv = defaultdict(list)
for r in calls:
    iv[r['acc']].append((int(r['begin']), int(r['end'])))
seq, name = {}, None
for l in open(sys.argv[3]):
    if l[0] == '>':
        name = l[1:].split()[0]; seq[name] = []
    else:
        seq[name].append(l.strip())
seq = {k: ''.join(v) for k, v in seq.items()}
loci = defaultdict(list)
for l in open(sys.argv[2]):
    q, s, pid, ln, qs, qe, ss, se, e, qa, sa = l.rstrip('\n').split('\t')
    a, b = sorted((int(ss), int(se)))
    if max([min(b, y) - max(a, x) for x, y in iv[s]] + [-1]) >= 30:
        continue
    loci[(s, a // 200)].append((q, int(qs), int(qe), a, b, int(ss) > int(se), float(e), qa, sa))
finder = pyaragorn.RNAFinder(4, tmrna=False, linear=True)
for (s, _), hs in sorted(loci.items()):
    if len(hs) < 5:
        continue
    lo, hi = min(h[3] for h in hs), max(h[4] for h in hs)
    trip = defaultdict(int)
    for q, qs, qe, a, b, rev, e, qa, sa in hs:
        r = qinfo[q]
        o = int(r['ac_off'])  # 1-based start of anticodon in query
        if qs <= o and o + 2 <= qe:
            i, k, out = qs, 0, ''
            for x, y in zip(qa, sa):
                if x != '-':
                    if o <= i <= o + 2:
                        out += y
                    i += 1
            trip[(r['anticodon'], out)] += 1
    region = seq[s][max(0, lo - 51):hi + 50]
    ar = [(g.amino_acid, g.anticodon, g.begin + max(0, lo - 51), g.strand, round(g.energy, 1))
          for g in finder.find_rna(region) if g.type == 'tRNA']
    qtypes = defaultdict(int)
    for h in hs:
        qtypes[h[0].split('|')[1]] += 1
    print(f'{s} {lo}-{hi} hits {len(hs)} by query anticodon {dict(qtypes)} best e {min(h[6] for h in hs):.1e} '
          f'max HSP {max(h[4] - h[3] + 1 for h in hs)} nt | aligned to query 34-36: {dict(trip)} | ARAGORN: {ar}')
