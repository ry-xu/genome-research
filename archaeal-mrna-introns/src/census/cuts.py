"""Census cuts of a ledger run (RULES.md R1-R5): one parquet per finished batch, resumable.
uv run --project /Users/ryan/lab/work/reading-ledger python src/cuts.py RUN OUTDIR"""
import os, sys, glob, json
import polars as pl
sys.path.insert(0, '/Users/ryan/lab/work/reading-ledger/snapshots/v0.7')
sys.path.insert(0, '/Users/ryan/lab/work/archaeal-intron-census/src')
import inserts                                   # frozen v0.7 cut rule, run with D_MIN 15 / D_MAX 300 (R1)
inserts.D_MIN, inserts.D_MAX, inserts.DCOL_MIN = 15, 300, -10      # dcol -10..-6 = "wide" tier (amendment 2)
sys.path.insert(0, '/Users/ryan/lab/work/archaeal-intron-census/src/checker')
from common import load_fasta, strand_seq, load_anchor
from splice_dp import analyse
from bhbset import slide_set, label, decoys, ok_window

FUSED = {'PGKT', 'RPOBC', 'ISPDF', 'HDRA', 'SECDF', 'MOACB', 'URE23', 'GSIA'}
run, outdir = sys.argv[1], sys.argv[2]
os.makedirs(outdir, exist_ok=True)
fa = {}
for lst in ('data/run1.tsv', 'data/run2.tsv', 'data/run3.tsv'):     # run3 = HC1/SC1 (amendment 3)
    if os.path.exists(lst):
        for l in open(lst):
            a, p = l.split('\t')[:2]; fa[a] = p
PROT = load_anchor()

def best(lst):
    return max(lst, key=lambda c: c[5]) if lst else None

def one(r, J):
    s = strand_seq(load_fasta(fa[r['acc']]), r['contig'], r['strand'])
    o = dict(r)
    excl = []
    if r['holds']: excl.append('holds')
    if r['has_n']: excl.append('N')
    if J['jtrunc']: excl.append('jtrunc')
    if r['fam'] in FUSED: excl.append('fused')
    fam, king, prot = PROT[r['sid']]
    cand = analyse(s, r['a_end'], r['b_start'], prot, r['a_se'], r['b_ss'])
    o.update(excl=','.join(excl), placed=False, tier='main' if r['dcol'] >= -5 else 'wide')
    if not cand:
        o['status'] = 'edge'; return o
    ex = [c for c in cand if c[0] == 'strict' and c[3] >= 9]
    fs = [c for c in cand if c[0] == 'strict' and c[3] in (-2, -1, 1, 2)]
    rt = [c for c in cand if c[0] == 'perm' and c[3] == 0]
    bex, bfs, brt = best(ex), best(fs), best(rt)
    o.update(ex_score=bex[5] if bex else None, fs_score=bfs[5] if bfs else None, rt_score=brt[5] if brt else None)
    if not bex:
        o['status'] = 'no_removal'; return o
    other = max([c[5] for c in (bfs, brt) if c] or [-1e9])
    x, y, L = bex[1], bex[2], bex[3]
    near = sorted({(c[1], c[2]) for c in ex if c[5] >= bex[5] - 5})
    o.update(x=x, y=y, L=L, phase=bex[4], col=bex[6], g=bex[7], jaa=bex[8], n_near=len(near), near=json.dumps(near),
             stretch=s[x:y], up20=s[x - 20:x], down20=s[y:y + 20])
    if bex[5] <= other:
        o['status'] = 'beaten'; return o
    if not ok_window(s, near):
        o['status'] = 'contig_edge'; return o
    _, k = slide_set(s, x, y)
    a, b = label(s, near)
    dec = decoys(s, near, L)
    o.update(status='placed', placed=True, k=k, bhb=a, strict=b, n_dec=len(dec),
             dec_bhb=sum(d[1] for d in dec), dec_strict=sum(d[2] for d in dec))
    return o

tags = sorted(os.path.basename(p)[:-8] for p in glob.glob(os.path.join(run, 'genomes', '*.parquet')))
for tag in tags:
    dst = os.path.join(outdir, f'{tag}.parquet')
    if os.path.exists(dst):
        continue
    Jt = pl.read_parquet(os.path.join(run, 'joins', f'{tag}.parquet'))
    rows = inserts.batch(run, tag)
    jk = {(j['acc'], j['contig'], j['strand'], j['a_end'], j['b_start']): j for j in
          Jt.select('acc', 'contig', 'strand', 'a_end', 'b_start', 'jtrunc', 'clen').iter_rows(named=True)}
    out = [one(r, jk[(r['acc'], r['contig'], r['strand'], r['a_end'], r['b_start'])]) for r in sorted(rows, key=lambda r: r['acc'])]
    G = pl.read_parquet(os.path.join(run, 'genomes', f'{tag}.parquet'))
    pl.DataFrame(out, infer_schema_length=None).write_parquet(dst + '.tmp') if out else pl.DataFrame({'acc': []}).write_parquet(dst + '.tmp')
    os.replace(dst + '.tmp', dst)
    print(tag, G.height, 'genomes', len(out), 'cuts', sum(o.get('placed', False) for o in out), 'placed', flush=True)
