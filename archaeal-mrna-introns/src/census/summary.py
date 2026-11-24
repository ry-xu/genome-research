"""Census tables (RULES.md R6, R7, C2, verdict; amendment 4). uv run --project ../reading-ledger python src/summary.py OUT CUTDIR...
Writes OUT/genomes_<tier>.tsv, OUT/classes_<tier>.tsv, OUT/sites_<tier>.tsv and prints the verdict lines."""
import glob, os, sys
import numpy as np
import polars as pl
from scipy.stats import fisher_exact
sys.path.insert(0, '/Users/ryan/lab/work/archaeal-intron-census/src/checker')
from common import load_anchor
from splice_dp import SUB, AIDX

NOT_CENSUS = ['AP029002.1', 'AP029247.1']                          # HC1/SC1: RNA step only (amendment 3, 4f)
NOT_CENSUS += [l.strip() for l in open('data/twins_excluded.txt') if l.strip()]   # GCA twins of GCF assemblies (amendment 5)
out = sys.argv[1]; os.makedirs(out, exist_ok=True)
C = pl.concat([pl.read_parquet(f) for d in sys.argv[2:] for f in sorted(glob.glob(f'{d}/*.parquet'))], how='diagonal_relaxed')
C = C.filter(pl.col('fam').is_not_null() & ~pl.col('acc').is_in(NOT_CENSUS))
meta = pl.read_csv('data/genomes_meta.tsv', separator='\t', infer_schema_length=0)
LED = {'out/cuts1': 'out/ledger1', 'out/cuts2': 'out/ledger2'}      # only batches whose cuts exist count
G = pl.concat([pl.read_parquet(f"{LED[d.rstrip('/')]}/genomes/{os.path.basename(f)}", columns=['acc', 'n_loci'])
               for d in sys.argv[2:] for f in sorted(glob.glob(f'{d}/*.parquet'))])
G = G.filter(~pl.col('acc').is_in(NOT_CENSUS)).join(meta.select('acc', 'phylum', 'class', 'order', 'species', 'rep'), on='acc', how='left')
assert G['acc'].n_unique() == G.height, 'a genome appears twice'                 # amendment 7e
assert G['class'].null_count() == 0 and G['phylum'].null_count() == 0, 'genome without class'
assert G['phylum'].null_count() == 0, G.filter(pl.col('phylum').is_null())['acc'].to_list()


def exp_pc(n, dpos, dn):
    """amendment 4(a): n cuts x (decoy BHB + 0.5) / (decoys + 1)"""
    return n * (dpos + 0.5) / (dn + 1.0)


# R7 reference columns (amendment 4c)
PROT = load_anchor()
W, MINSC = 15, 30.0


def enc(s):
    return np.array([AIDX.get(a, 20) for a in s], dtype=np.int64)


def refcol(sid, col, ref):
    """1-based column `col` of anchor sid -> 1-based column of anchor ref by the best ungapped 31-aa window match."""
    if sid == ref:
        return col, None
    q = PROT[sid][2]; r = PROT[ref][2]
    lo, hi = max(0, col - 1 - W), min(len(q), col + W)
    qw, rr = enc(q[lo:hi]), enc(r)
    if len(qw) < 2 * W or len(rr) < len(qw):
        return None, None
    sc = SUB[qw[None, :], np.lib.stride_tricks.sliding_window_view(rr, len(qw))].sum(1)
    bo = int(np.argmax(sc)); best = sc[bo]
    if best < MINSC:
        return None, float(best)
    return bo + (col - 1 - lo) + 1, float(best)


def kmers12(s):
    return {s[i:i + 12] for i in range(len(s) - 11)}


for tier in ('main', 'wide'):
    k = C.filter((pl.col('tier') == tier) & pl.col('placed') & (pl.col('excl') == ''))
    nodec = k.filter(pl.col('n_dec') == 0).height
    k = k.filter(pl.col('n_dec') > 0).with_columns(rate=pl.col('dec_bhb') / pl.col('n_dec').cast(pl.Float64))
    lo = k.filter(pl.col('k') <= 3)
    per = lo.group_by('acc').agg(n_placed=pl.len(), obs=pl.col('bhb').sum(), exp_raw=pl.col('rate').sum(),
                                 dpos=pl.col('dec_bhb').sum(), dn=pl.col('n_dec').sum(),
                                 obs_strict=pl.col('strict').sum(), bhb_fams=pl.col('fam').filter(pl.col('bhb')).unique().sort().str.join(','))
    g = G.join(per, on='acc', how='left').with_columns([pl.col(c).fill_null(0) for c in ('n_placed', 'obs', 'exp_raw', 'dpos', 'dn', 'obs_strict')]
                                                       ).with_columns(pl.col('bhb_fams').fill_null(''))
    g.write_csv(f'{out}/genomes_{tier}.tsv', separator='\t')
    rows = []
    groups = sorted(g.partition_by(['phylum', 'class'], as_dict=True).items(), key=lambda kv: (kv[0][0] or '', kv[0][1] or ''))
    for ci_, ((ph, cl), sub) in enumerate(groups):
        o, n = int(sub['obs'].sum()), int(sub['n_placed'].sum())
        dpos, dn = int(sub['dpos'].sum()), int(sub['dn'].sum())
        e = exp_pc(n, dpos, dn) if n else 0.0
        p = float(fisher_exact([[o, n - o], [dpos, dn - dpos]], alternative='greater')[1]) if n else 1.0
        ob, nb, pb, db = (sub[c].to_numpy() for c in ('obs', 'n_placed', 'dpos', 'dn'))
        if sub.height > 1 and n:
            rng = np.random.default_rng(1 + ci_)
            idx = rng.integers(0, sub.height, size=(2000, sub.height))
            r = ob[idx].sum(1) / exp_pc(nb[idx].sum(1), pb[idx].sum(1), db[idx].sum(1)).clip(min=1e-9)
            r = r[nb[idx].sum(1) > 0]
            ci = (float(np.quantile(r, 0.025)), float(np.quantile(r, 0.975))) if len(r) else (float('nan'),) * 2
        else:
            ci = (float('nan'), float('nan'))
        fams = {f for s in sub['bhb_fams'].to_list() for f in s.split(',') if f}
        rows.append(dict(phylum=ph, cls=cl, genomes=sub.height, loci=int(sub['n_loci'].sum()), placed=n, obs=o,
                         dec_pos=dpos, decoys=dn, exp=round(e, 2), ratio=round(o / e, 1) if e else None,
                         ci_lo=ci[0], ci_hi=ci[1], p=p, bhb_per_1k_loci=round(1000 * o / max(1, int(sub['n_loci'].sum())), 3),
                         genomes_with_bhb=int((sub['obs'] > 0).sum()), fams=len(fams)))
    T = pl.DataFrame(rows).sort(['obs', 'phylum', 'cls'], descending=[True, False, False])
    m = T.filter(pl.col('placed') >= 5).height
    T = T.with_columns(p_bonf=pl.min_horizontal(pl.col('p') * max(m, 1), pl.lit(1.0)), tested=pl.col('placed') >= 5)
    T.write_csv(f'{out}/classes_{tier}.tsv', separator='\t')
    pl.Config.set_tbl_rows(60); pl.Config.set_tbl_width_chars(250)
    k4 = k.filter(pl.col('k') >= 4)
    print(f'== tier {tier}: placed kept cuts {k.height} (+{nodec} with 0 decoys, left out), k<=3 {lo.height}; BHB {int(lo["bhb"].sum())}, '
          f'decoy BHB {int(lo["dec_bhb"].sum())}/{int(lo["n_dec"].sum())}, exp {exp_pc(lo.height, lo["dec_bhb"].sum(), lo["n_dec"].sum()):.1f}; '
          f'k>=4: {k4.height} cuts, BHB {int(k4["bhb"].sum())}, exp {exp_pc(k4.height, k4["dec_bhb"].sum(), k4["n_dec"].sum()):.1f}; Bonferroni m={m}')
    print(T.select('phylum', 'cls', 'genomes', 'placed', 'obs', 'exp', 'ratio', 'ci_lo', 'ci_hi', 'tested', 'p_bonf', 'bhb_per_1k_loci', 'genomes_with_bhb', 'fams').head(40))
    A = T.filter(pl.col('phylum') == 'Asgardarchaeota')
    if A.height:
        ao, an = int(A['obs'].sum()), int(A['placed'].sum())
        ae = exp_pc(an, int(A['dec_pos'].sum()), int(A['decoys'].sum()))
        sig = A.filter(pl.col('tested') & (pl.col('p_bonf') < 0.05)).height
        line = f'Asgard pooled obs {ao} of {an} placed, exp {ae:.2f}, ratio {ao / ae if ae else float("inf"):.1f}; classes passing: {sig}'
        if tier == 'main':
            line += f'; VERDICT BAR: {"MET" if ao >= 20 and ae and ao / ae >= 3 and sig >= 2 else "not met"}'
        else:
            line += ' (wide tier: reported only, no verdict)'
        print(line)
    if tier == 'main':                                   # amendment 6: quality strata, reporting only
        q = g.join(meta.select('acc', 'compl', 'contam', 'contigs'), on='acc', how='left').with_columns(
            pl.col('compl').cast(pl.Float64), pl.col('contam').cast(pl.Float64), pl.col('contigs').cast(pl.Int64),
            asg=pl.col('phylum') == 'Asgardarchaeota')
        hq = (pl.col('compl') >= 90) & (pl.col('contam') <= 5)
        qr = {}
        for nm, f in (('(i) compl>=90 contam<=5', hq), ('(ii) rest', ~hq), ('(iii) <=10 contigs', pl.col('contigs') <= 10)):
            for asg in (True, False):
                x = q.filter(f & (pl.col('asg') == asg))
                o_, n_ = int(x['obs'].sum()), int(x['n_placed'].sum())
                e_ = exp_pc(n_, int(x['dpos'].sum()), int(x['dn'].sum())) if n_ else 0.0
                qr[(nm[:4], asg)] = o_ / e_ if e_ else float('nan')
                print(f"  quality {nm} {'Asgard' if asg else 'other '}: genomes {x.height}, placed {n_}, obs {o_}, exp {e_:.2f}, "
                      f"ratio {o_ / e_ if e_ else float('nan'):.1f}")
        r1, r2 = qr[('(i) ', True)], qr[('(ii)', True)]          # amendment 6 reading, (ii) = Asgard's (amendment 7c)
        print(f"  amendment 6: Asgard (i) {r1:.1f} vs 3 and vs a third of Asgard (ii) {r2 / 3:.1f} -> "
              f"{'excess LEANS on lower-quality MAGs' if not (r1 >= 3 and r1 >= r2 / 3) else 'excess does not lean on lower-quality MAGs'}")
    # R7 recurrent sites on reference columns
    b = lo.filter(pl.col('bhb')).join(meta.select('acc', 'phylum', 'class'), on='acc', how='left')
    sites, unm = [], 0
    for fam in sorted(b['fam'].unique().to_list()):
        s = b.filter(pl.col('fam') == fam)
        ref = s.group_by('sid').len().sort(['len', 'sid'], descending=[True, False])['sid'][0]
        mapped = []
        for r in s.iter_rows(named=True):
            rc_, sc = refcol(r['sid'], r['col'], ref)
            if rc_ is None:
                unm += 1; continue
            mapped.append(dict(r, rcol=rc_))
        if not mapped:
            continue
        M = pl.DataFrame(mapped, infer_schema_length=None)
        for ph_ in sorted(M['phase'].unique().to_list()):
            mp = M.filter(pl.col('phase') == ph_).sort('rcol')
            cols = mp['rcol'].to_list(); seen = set()
            for c0 in sorted(set(cols)):                         # every window centre; identical member sets kept once
                grp = mp.filter((pl.col('rcol') - c0).abs() <= 3)
                key = tuple(sorted(zip(grp['acc'].to_list(), grp['x'].to_list())))
                if key in seen or grp['acc'].n_unique() < 3:
                    continue
                seen.add(key)
                st = grp['stretch'].to_list(); ks = [kmers12(x) for x in st]
                alike = sum(any(i != j and ks[i] & ks[j] for j in range(len(ks))) for i in range(len(ks))) / len(ks)
                sites.append(dict(fam=fam, ref_sid=ref, phase=ph_, rcol=int(np.median(grp['rcol'].to_numpy())), genomes=grp['acc'].n_unique(),
                                  classes=grp['class'].n_unique(), phyla=';'.join(sorted(set(x or '' for x in grp['phylum'].to_list()))),
                                  class_list=';'.join(sorted(set(x or '' for x in grp['class'].to_list()))), L_med=float(grp['L'].median()),
                                  dec_rate=float(grp['rate'].mean()), alike12=round(alike, 2)))
    print(f'R7: {b.height} BHB cuts, {unm} not mapped to a reference column')
    if sites:
        S = pl.DataFrame(sites).sort(['classes', 'genomes'], descending=True)
        S.write_csv(f'{out}/sites_{tier}.tsv', separator='\t')
        print(S.filter(pl.col('classes') >= 2).head(40))
# C2
c2 = {'A. pernix': 'GCF_000011125.1', 'S. solfataricus (rep)': 'GCF_900079115.1', 'S. acidocaldarius': 'GCF_000012285.1', 'M. sedula': 'GCF_000016605.1'}
for nm, a in c2.items():
    s = C.filter((pl.col('acc') == a) & (pl.col('fam') == 'TRUB'))
    done = a in set(G['acc'].to_list())
    print('C2', nm, a, 'not yet run' if not done else
          (s.select('tier', 'd', 'dcol', 'status', 'L', 'k', 'bhb', 'strict').to_dicts() or 'no TRUB cut'))
c2h = C.filter(pl.col('acc').is_in(list(c2.values())) & (pl.col('fam') == 'TRUB') & pl.col('placed') & pl.col('bhb') & (pl.col('excl') == ''))
c2m, c2w = c2h.filter(pl.col('tier') == 'main')['acc'].n_unique(), c2h['acc'].n_unique()
print(f"C2 (amendment 7d): main tier {c2m} of 4 species ({sum(a in set(G['acc'].to_list()) for a in c2.values())} run) -> "
      f"{'PASS' if c2m >= 3 else 'not passed'}; with wide tier {c2w} of 4 (reported only)")
