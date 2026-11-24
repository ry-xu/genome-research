"""Window index (RULES.md §1): census cuts (main, placed, kept, k<=3), tRNA introns, same-shape decoys.
Windows (event +-300) are merged per contig into regions in forward orientation.
uv run --project /Users/ryan/lab/work/reading-ledger python src/build_index.py  -> out/regions.fa, out/events.tsv"""
import glob, sys, os
import polars as pl
CEN = '/Users/ryan/lab/work/archaeal-intron-census'
sys.path.insert(0, f'{CEN}/src'); sys.path.insert(0, f'{CEN}/src/checker')
from common import load_fasta, rc
from bhbset import slide_set

W, DEC_GAP, DEC_SEP = 300, 60, 30
TW = [l.strip() for l in open(f'{CEN}/data/twins_excluded.txt')]
C = pl.concat([pl.read_parquet(f) for f in sorted(glob.glob(f'{CEN}/out/cuts1/*.parquet') + glob.glob(f'{CEN}/out/cuts2/*.parquet'))],
              how='diagonal_relaxed').filter(~pl.col('acc').is_in(TW))
C = C.filter(pl.col('placed') & (pl.col('excl') == '') & (pl.col('tier') == 'main') & (pl.col('k') <= 3))
M = pl.read_csv(f'{CEN}/data/genomes_meta.tsv', separator='\t').select('acc', 'phylum', 'class', 'order', 'family', 'genus', 'species')
C = C.join(M, on='acc', how='left')
fa = {}
for f in ('run1.tsv', 'run2.tsv', 'run3.tsv'):
    p = f'{CEN}/data/{f}'
    if os.path.exists(p):
        for l in open(p):
            a = l.rstrip('\n').split('\t')
            fa.setdefault(a[0], a[1])
P = '/Users/ryan/lab/shared/codon-atlas/products'

ev, regs = [], []
for acc in sorted(C['acc'].unique().to_list()):
    seqs = load_fasta(fa[acc])
    cuts = C.filter(pl.col('acc') == acc).to_dicts()
    E = []   # (kind, contig, fx, fy, k, meta)
    for i, r in enumerate(cuts):
        n = len(seqs[r['contig']])
        fx, fy = (r['x'], r['y']) if r['strand'] == 0 else (n - r['y'], n - r['x'])
        cid = f"{acc}:c{i}"
        E.append(('cut', r['contig'], fx, fy, r['k'], dict(id=cid, parent=cid, strand=r['strand'], bhb=r['bhb'], strict=r['strict'],
                  fam=r['fam'], col=r['col'], phase=r['phase'])))
    tp = f'{P}/{acc[:8]}/{acc}/trnascan.parquet'
    if os.path.exists(tp):
        t = pl.read_parquet(tp).filter(pl.col('intron_begin').is_not_null())
        j = 0
        for r in t.iter_rows(named=True):
            if r['contig'] not in seqs:
                continue
            for ib, ie in zip(r['intron_begin'], r['intron_end']):
                fx, fy = min(ib, ie) - 1, max(ib, ie)
                _, k = slide_set(seqs[r['contig']], fx, fy)
                tid = f"{acc}:t{j}"; j += 1
                E.append(('trna', r['contig'], fx, fy, k, dict(id=tid, parent=tid, strand=r['strand'], bhb=None, strict=None,
                          fam=f"tRNA-{r['type']}", col=None, phase=None)))
    occupied = [(c, fx, fy) for _, c, fx, fy, _, _ in E]
    for kind, c, fx, fy, k, m in list(E):
        if kind != 'cut':
            continue
        n = len(seqs[c]); L = fy - fx
        for sgn, tag in ((-1, 'm'), (1, 'p')):
            d = sgn * (L + DEC_GAP)
            dx, dy = fx + d, fy + d
            if dx < 0 or dy > n:
                continue
            if any(oc == c and dx < oy + DEC_SEP and dy > ox - DEC_SEP for oc, ox, oy in occupied):
                continue
            _, dk = slide_set(seqs[c], dx, dy)
            E.append(('decoy', c, dx, dy, dk, dict(m, id=m['id'] + tag, bhb=None, strict=None)))
    # merge windows per contig
    byc = {}
    for e in E:
        byc.setdefault(e[1], []).append(e)
    for c, es in byc.items():
        n = len(seqs[c])
        iv = sorted((max(0, fx - W), min(n, fy + W)) for _, _, fx, fy, _, _ in es)
        merged = []
        for a, b in iv:
            if merged and a <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], b)
            else:
                merged.append([a, b])
        for a, b in merged:
            rid = f"{acc}|{c}|{a}"
            regs.append((rid, seqs[c][a:b]))
            for kind, cc, fx, fy, k, m in es:
                if a <= fx and fy <= b:
                    ev.append(dict(region=rid, kind=kind, acc=acc, contig=c, rx=fx - a, ry=fy - a, L=fy - fx, k=k, **m))
E = pl.DataFrame(ev).join(M, on='acc', how='left')
E.write_csv('out/events.tsv', separator='\t')
with open('out/regions.fa', 'w') as f:
    for rid, s in regs:
        f.write(f'>{rid}\n{s}\n')
print('regions', len(regs), 'bp', sum(len(s) for _, s in regs))
print(E.group_by('kind').agg(n=pl.len(), g=pl.col('acc').n_unique()))
