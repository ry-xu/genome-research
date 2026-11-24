"""Post hoc, descriptive only (not a census rule): R2's competitor models (frameshift by 1-2 nt, stop read-through) exist only
when b_start - a_end <= ~92, because splice_dp.analyse searches x and y within +-45 of the table ends. For placed cuts with a
longer gap, rerun analyse with a window wide enough to reach across the gap and ask whether a competitor beats the census removal.
uv run --project /Users/ryan/lab/work/reading-ledger python src/longcut_check.py OUT.tsv CUTDIR..."""
import glob, sys
import polars as pl
sys.path.insert(0, '/Users/ryan/lab/work/archaeal-intron-census/src/checker')
from common import load_fasta, strand_seq, load_anchor
from splice_dp import analyse

fa = {}
for lst in ('data/run1.tsv', 'data/run2.tsv', 'data/run3.tsv'):
    try:
        for l in open(lst):
            a, p = l.split('\t')[:2]; fa[a] = p
    except FileNotFoundError:
        pass
PROT = load_anchor()
C = pl.concat([pl.read_parquet(f) for d in sys.argv[2:] for f in sorted(glob.glob(f'{d}/*.parquet'))], how='diagonal_relaxed')
C = C.filter(pl.col('placed') & (pl.col('b_start') - pl.col('a_end') > 92)).sort('acc')
rows, cache = [], {}
for r in C.iter_rows(named=True):
    if r['acc'] not in cache:
        cache = {r['acc']: load_fasta(fa[r['acc']])}
    s = strand_seq(cache[r['acc']], r['contig'], r['strand'])
    gap = r['b_start'] - r['a_end']
    win = gap // 2 + 45; ext = win // 3 + 2
    cand = analyse(s, r['a_end'], r['b_start'], PROT[r['sid']][2], r['a_se'], r['b_ss'], EXT=ext, ND=90, WIN=win)
    fs = max([c[5] for c in cand or [] if c[0] == 'strict' and c[3] in (-2, -1, 1, 2)] or [None], key=lambda v: -1e9 if v is None else v)
    rt = max([c[5] for c in cand or [] if c[0] == 'perm' and c[3] == 0] or [None], key=lambda v: -1e9 if v is None else v)
    ex = max([c[5] for c in cand or [] if c[0] == 'strict' and c[3] >= 9] or [None], key=lambda v: -1e9 if v is None else v)
    other = max([v for v in (fs, rt) if v is not None] or [None], key=lambda v: -1e9 if v is None else v)
    rows.append(dict(acc=r['acc'], fam=r['fam'], tier=r['tier'], gap=gap, L=r['L'], k=r['k'], bhb=r['bhb'], excl=r['excl'],
                     census_ex=r['ex_score'], wide_ex=ex, fs=fs, rt=rt, beaten=other is not None and other >= r['ex_score']))
O = pl.DataFrame(rows, infer_schema_length=None)
O.write_csv(sys.argv[1], separator='\t')
print(O.group_by('tier', 'bhb').agg(n=pl.len(), with_competitor=(pl.col('fs').is_not_null() | pl.col('rt').is_not_null()).sum(),
                                     beaten=pl.col('beaten').sum()).sort('tier', 'bhb'))
print(O.filter(pl.col('bhb') | pl.col('beaten')).sort('bhb', 'beaten'))
