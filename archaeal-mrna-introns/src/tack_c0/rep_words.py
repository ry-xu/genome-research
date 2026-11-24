"""A5 words (RULES.md A5): replication of the c0 splice in the PRJNA992765 Methanosuratincola cultures.
usage: rep_words.py RNA_RUN,RUN,... METAG_RUN [READCHECK.tsv]   (counts from out/rep/runs, index out/rep)
Prints per-homolog words pooled over the RNA runs, decoys, metaG, the A3 read check, the A5 verdict, and the tRNA-intron
diagnostic split by read orientation (sense / antisense to the feature, mate 1 orientation as in strand_check.sh)."""
import sys, csv, os
from collections import defaultdict
D = '/Users/ryan/lab/work/tack-rna/out/rep'
rna = sys.argv[1].split(','); mg_run = sys.argv[2]; rck = sys.argv[3] if len(sys.argv) > 3 else None
EV = {r['id']: r for r in csv.DictReader(open(f'{D}/events.tsv'), delimiter='\t')}
def pooled(runs):
    t = defaultdict(lambda: [0, 0])
    for r in runs:
        assert os.path.exists(f'{D}/runs/{r}.done'), r
        for x in csv.DictReader(open(f'{D}/runs/{r}.counts.tsv'), delimiter='\t'):
            t[x['event']][0] += int(x['joined']); t[x['event']][1] += int(x['through'])
    return t
def word(j, th):
    if j >= 2:
        return 'removed' if j / (j + th) >= 0.5 else 'partly removed'
    return 'read through' if th >= 5 else 'not covered'
R, M = pooled(rna), pooled([mg_run])
print('== homologs (RNA pooled over', len(rna), 'runs; metaG', mg_run + ')')
print('id\tL\tjoined\tthrough\tword\tdecoys(j/t)\tmetaG(j/t)')
cuts = [e for e in EV.values() if e['kind'] == 'cut']
for e in sorted(cuts, key=lambda e: e['id']):
    i = e['id']; j, th = R[i]
    dec = ' '.join(f"{d}:{R[d][0]}/{R[d][1]}" for d in sorted(EV) if EV[d]['kind'] == 'decoy' and EV[d]['parent'] == i)
    print(f"{i}\t{e['L']}\t{j}\t{th}\t{word(j, th)}\t{dec}\t{M[i][0]}/{M[i][1]}")
best = max(cuts, key=lambda e: (sum(R[e['id']]), e['id']))['id']
j, th = R[best]; w = word(j, th)
rc_j = rc_x = 0
if rck:
    for x in csv.DictReader(open(rck), delimiter='\t'):
        if x['event'] == best:
            rc_j += int(x['joined_templates']); rc_x += int(x['explained_templates'])
dec_sp = [d for d in EV if EV[d]['kind'] == 'decoy' and EV[d]['parent'] == best and R[d][0] >= 2]
mj, mt = M[best]
dna_ok = mj <= 1 and mj + mt >= 5
print(f'== best-covered homolog {best}: joined {j} through {th} -> {w}; spliced decoys {dec_sp or "none"}; '
      f'metaG {mj}/{mt} -> {"clean" if dna_ok else "FAIL/uncovered"}; read check {rc_x} of {rc_j} joined templates explained'
      f'{"" if rck else " (not run)"}')
if max(sum(R[e['id']]) for e in cuts) < 5:
    v = 'NOT COVERED'
elif w == 'removed':
    v = 'REPLICATES' if (not dec_sp and dna_ok and rck and rc_j and rc_x == 0) else \
        'removed, but: ' + ', '.join(s for s, bad in (('decoy spliced', dec_sp), ('DNA check failed', not dna_ok),
                                                      ('read check pending', not rck), ('joined reads explained', rc_x)) if bad)
elif w == 'read through':
    v = 'DOES NOT REPLICATE'
else:
    v = w
print('   A5 verdict:', v)
# tRNA diagnostic with orientation
ori = defaultdict(lambda: defaultdict(int))
for r in rna:
    for x in csv.DictReader(open(f'{D}/runs/{r}.reads.tsv'), delimiter='\t'):
        e = EV[x['event']]
        if e['kind'] not in ('trna', 'cut'):
            continue
        rev = (int(x['flag']) // 16) % 2
        r1rev = rev if x['read'].endswith('1') else 1 - rev
        fwd = (int(e['strand']) == 1) if e['kind'] == 'trna' else (int(e['strand']) == 0)
        o = 'sense' if r1rev == fwd else 'antisense'
        ori[(e['kind'], e['id'] if e['kind'] == 'cut' else 'tRNA')][(o, x['cls'])] += 1
print('== orientation of spanning reads (mate 1 vs feature; dUTP libraries are mostly antisense)')
for k in sorted(ori):
    print(' ', k, dict(sorted(ori[k].items())))
T = [e for e in EV.values() if e['kind'] == 'trna']
cov = [e for e in T if sum(R[e['id']]) >= 5]
print(f'== tRNA introns: {len(T)} in index, covered (>=5 spanning) {len(cov)}, spliced (joined>=2) {sum(R[e["id"]][0] >= 2 for e in cov)}, '
      f'removed {sum(word(*R[e["id"]]) == "removed" for e in cov)}')
for e in sorted(cov, key=lambda e: e['id']):
    print(f"  {e['id']}\t{e['fam']}\tL {e['L']}\t{R[e['id']][0]}/{R[e['id']][1]}\t{word(*R[e['id']])}")
