"""K1 (RULES.md, amendment 1): census cuts in the RNA survey's genomes against the survey's removals tables and
through-read coverage. No new mapping. -> out/k1_survey.tsv"""
import csv, glob, os, sys
import numpy as np
import polars as pl
SV = os.environ.get('K1_SV', '/Users/ryan/lab/work/archaeal-introns/survey')   # K1_SV=. for HC1/SC1 (out/units, out/rna)
RUNS = os.environ.get('K1_RUNS', 'out/runs')
sys.path.insert(0, '/Users/ryan/lab/work/archaeal-intron-census/src/checker')
from common import load_fasta
fa = {}
for lst in ('data/run1.tsv', 'data/run2.tsv', 'data/run3.tsv'):
    if os.path.exists(lst):
        fa.update(dict(l.split('\t')[:2] for l in open(lst)))
OUT = sys.argv[-1] if sys.argv[-1].endswith('.tsv') else 'out/k1_survey.tsv'
DIRS = [a for a in sys.argv[1:] if not a.endswith('.tsv')]
D = pl.concat([pl.read_parquet(f) for d in DIRS for f in sorted(glob.glob(f'{d}/*.parquet'))], how='diagonal_relaxed')
units = {os.path.basename(p)[:-13] for p in glob.glob(f'{SV}/out/units/*.removals.tsv')}
D = D.filter(pl.col('acc').is_in(list(units)) & pl.col('placed'))
out = []
for acc in D['acc'].unique().sort():
    rem = list(csv.DictReader(open(f'{SV}/out/units/{acc}.removals.tsv'), delimiter='\t'))
    runs = [l.rstrip('\n').split('\t') for l in open(f'{SV}/out/units/{acc}.runs.tsv')]
    covs = [np.load(f'{SV}/{RUNS}/{os.path.basename(r[0])}.cov.npz') for r in runs if r[3] == '1']
    lens = {c: len(s) for c, s in load_fasta(fa[acc]).items()}
    for r in D.filter(pl.col('acc') == acc).iter_rows(named=True):
        n = lens[r['contig']]
        fs, fe = (r['x'], r['y']) if r['strand'] == 0 else (n - r['y'], n - r['x'])
        hit = [m for m in rem if m['contig'] == r['contig'] and abs(int(m['start']) - fs) <= r['k'] + 3 and abs(int(m['L']) - r['L']) <= 3]
        best = max(hit, key=lambda m: int(m['templates']), default=None)
        assert all(r['contig'] in c.files for c in covs), (acc, r['contig'])
        thr = min(sum(int(c[r['contig']][p]) for c in covs) for p in (fs, fe))
        tm = int(best['templates']) if best else 0
        # amendment 4(e): only 'no row' with >= 10 through-reads is killed; a row with < 2 templates is 'not covered'
        word = 'removed in RNA' if tm >= 2 else ('read through' if (best is None and thr >= 10) else 'not covered')
        out.append(dict(acc=acc, contig=r['contig'], strand=r['strand'], tier=r['tier'], short=r['L'] < 15, fstart=fs, L=r['L'], k=r['k'], fam=r['fam'], col=r['col'], bhb=r['bhb'],
                        strict=r['strict'], excl=r['excl'], dec_bhb=r['dec_bhb'], n_dec=r['n_dec'], templates=tm,
                        fragments=int(best['fragments']) if best else 0, row_start=best['start'] if best else '',
                        row_cls=best['cls'] if best else '', row_bhb=best['bhb'] if best else '', through=thr, word=word))
O = pl.DataFrame(out, infer_schema_length=None)
O.write_csv(OUT, separator='\t')
print(O.group_by('tier', 'bhb', 'word').len().sort('tier', 'bhb', 'word'))
