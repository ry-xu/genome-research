"""C1: the BHB scorer on tRNAscan-SE verified tRNA introns in the census genomes, with R5 decoys.
uv run python src/c1_trna.py  -> out/c1_trna.tsv"""
import os, sys, csv
import polars as pl
sys.path.insert(0, '/Users/ryan/lab/work/archaeal-intron-census/src')
from bhbset import slide_set, label, decoys, ok_window
from checker.common import load_fasta, rc

P = '/Users/ryan/lab/shared/codon-atlas/products'
rows = []
accs = [l.split('\t')[:2] for l in open('data/run1.tsv')]
n_g = 0
for gi, (acc, fa) in enumerate(accs):
    tp = f'{P}/{acc[:8]}/{acc}/trnascan.parquet'
    if not os.path.exists(tp):
        continue
    t = pl.read_parquet(tp).filter(pl.col('intron_begin').is_not_null())
    if not t.height:
        continue
    n_g += 1
    seqs = load_fasta(fa)
    for r in t.iter_rows(named=True):
        if r['contig'] not in seqs:
            continue
        g = seqs[r['contig']]; n = len(g)
        s = g if r['strand'] == 1 else rc(g)
        for ib, ie in zip(r['intron_begin'], r['intron_end']):
            x, y = (ib - 1, ie) if r['strand'] == 1 else (n - ie, n - ib + 1)
            L = y - x
            pairs, k = slide_set(s, x, y)
            if not ok_window(s, pairs):
                continue
            a, b = label(s, pairs)
            dec = decoys(s, pairs, L)
            rows.append(dict(acc=acc, contig=r['contig'], type=r['type'], anticodon=r['anticodon'], note=r['note'] or '',
                             L=L, k=k, bhb=a, strict=b, n_dec=len(dec), dec_bhb=sum(d[1] for d in dec),
                             dec_strict=sum(d[2] for d in dec)))
    if gi % 200 == 0:
        print(gi, acc, len(rows), file=sys.stderr, flush=True)
D = pl.DataFrame(rows)
D.write_csv('out/c1_trna.tsv', separator='\t')
print('genomes with intron tRNAs', n_g, 'introns', D.height)
for nm, f in (('all', D), ('k<=3', D.filter(pl.col('k') <= 3)), ('k>=4', D.filter(pl.col('k') >= 4)),
              ('C1 set: CI-only tRNAs', D.filter(pl.col('note') == 'CI')), ('NCI-only', D.filter(pl.col('note') == 'NCI')),
              ('mixed/other', D.filter(~pl.col('note').is_in(['CI', 'NCI'])))):
    if f.height:
        print(f"{nm}: n={f.height} BHB={f['bhb'].mean():.3f} strict={f['strict'].mean():.3f} "
              f"decoy BHB={f['dec_bhb'].sum()/max(1,f['n_dec'].sum()):.3f} decoy strict={f['dec_strict'].sum()/max(1,f['n_dec'].sum()):.3f}")
