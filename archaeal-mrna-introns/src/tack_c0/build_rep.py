"""Replication index (RULES.md A5): the Y674 insertion homolog in the four Methanosuratincola genomes, same-shape decoys
(build_index.py rule: shifted by +-(L+60), kept off other events by 30 nt) and the genomes' tRNA introns.
The insertion is located by its exact sequence (from review/killcheck_c0/homolog_orfs.out) on either strand.
uv run --project /Users/ryan/lab/work/reading-ledger python -I src/build_rep.py  -> out/rep/regions.fa, out/rep/events.tsv"""
import sys, os
import polars as pl
CEN = '/Users/ryan/lab/work/archaeal-intron-census'
sys.path.insert(0, f'{CEN}/src'); sys.path.insert(0, f'{CEN}/src/checker')
from common import load_fasta, rc
from bhbset import slide_set
os.chdir('/Users/ryan/lab/work/tack-rna')
W, DEC_GAP, DEC_SEP = 300, 60, 30
INS = {'GCA_049570475.1': 'ACGGATCTAGAGTCCGGCCGAA', 'GCA_038801655.1': 'ACGGATCTAAAGTCCGGCCGAA',
       'GCA_030522375.1': 'ACGGGTGAGAGAATCCGGCCGAA', 'GCA_038801705.1': 'ACGGATGAGATCCGGCCGAA'}
FLANK = 'CCCAGCTTCACGGTAT'          # 16 nt upstream of every homolog (homolog_orfs/frame_homologs output)
M = pl.read_csv(f'{CEN}/data/genomes_meta.tsv', separator='\t').select('acc', 'phylum', 'class', 'order', 'family', 'genus', 'species')
fa = {}
for f in ('run1.tsv', 'run2.tsv', 'run3.tsv'):
    for l in open(f'{CEN}/data/{f}'):
        a = l.rstrip('\n').split('\t')
        fa.setdefault(a[0], a[1])
P = '/Users/ryan/lab/shared/codon-atlas/products'
ev, regs = [], []
for acc, ins in INS.items():
    seqs = load_fasta(fa[acc])
    hits = []
    for c, s in seqs.items():
        n = len(s)
        for strand, q in ((0, s), (1, rc(s))):
            i = q.find(FLANK + ins)
            while i >= 0:
                x = i + len(FLANK); y = x + len(ins)
                hits.append((c, strand, (x, y) if strand == 0 else (n - y, n - x)))
                i = q.find(FLANK + ins, i + 1)
    assert len(hits) == 1, (acc, hits)
    c, strand, (fx, fy) = hits[0]
    _, k = slide_set(seqs[c], fx, fy)
    cid = f'{acc}:y0'
    E = [('cut', c, fx, fy, k, dict(id=cid, parent=cid, strand=strand, bhb=True, strict=False, fam='Y674', col=None, phase=None))]
    tp = f'{P}/{acc[:8]}/{acc}/trnascan.parquet'
    if os.path.exists(tp):
        t = pl.read_parquet(tp).filter(pl.col('intron_begin').is_not_null())
        j = 0
        for r in t.iter_rows(named=True):
            if r['contig'] not in seqs:
                continue
            for ib, ie in zip(r['intron_begin'], r['intron_end']):
                tx, ty = min(ib, ie) - 1, max(ib, ie)
                _, tk = slide_set(seqs[r['contig']], tx, ty)
                tid = f'{acc}:t{j}'; j += 1
                E.append(('trna', r['contig'], tx, ty, tk, dict(id=tid, parent=tid, strand=r['strand'], bhb=None, strict=None,
                          fam=f"tRNA-{r['type']}", col=None, phase=None)))
    occupied = [(cc, a, b) for _, cc, a, b, _, _ in E]
    kind, cc, a0, b0, _, m = E[0]
    n = len(seqs[cc]); L = b0 - a0
    for sgn, tag in ((-1, 'm'), (1, 'p')):
        d = sgn * (L + DEC_GAP)
        dx, dy = a0 + d, b0 + d
        if dx < 0 or dy > n or any(oc == cc and dx < oy + DEC_SEP and dy > ox - DEC_SEP for oc, ox, oy in occupied):
            continue
        _, dk = slide_set(seqs[cc], dx, dy)
        E.append(('decoy', cc, dx, dy, dk, dict(m, id=m['id'] + tag, bhb=None, strict=None)))
    byc = {}
    for e in E:
        byc.setdefault(e[1], []).append(e)
    for c2, es in byc.items():
        n = len(seqs[c2])
        iv = sorted((max(0, a - W), min(n, b + W)) for _, _, a, b, _, _ in es)
        merged = []
        for a, b in iv:
            if merged and a <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], b)
            else:
                merged.append([a, b])
        for a, b in merged:
            rid = f'{acc}|{c2}|{a}'
            regs.append((rid, seqs[c2][a:b]))
            for kd, c3, x, y, kk, mm in es:
                if a <= x and y <= b:
                    ev.append(dict(region=rid, kind=kd, acc=acc, contig=c3, rx=x - a, ry=y - a, L=y - x, k=kk, **mm))
    print(acc, c, strand, fx, fy, 'k', k)
os.makedirs('out/rep', exist_ok=True)
E = pl.DataFrame(ev).join(M, on='acc', how='left')
E.write_csv('out/rep/events.tsv', separator='\t')
with open('out/rep/regions.fa', 'w') as f:
    for rid, s in regs:
        f.write(f'>{rid}\n{s}\n')
print('regions', len(regs), 'bp', sum(len(s) for _, s in regs))
print(E.group_by('acc', 'kind').agg(n=pl.len()).sort('acc', 'kind'))
