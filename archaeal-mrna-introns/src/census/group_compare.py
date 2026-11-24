"""Post hoc, descriptive: per-group BHB-cut rates (main tier, k<=3, kept, placed), with TRUB and/or Bathyarchaeia removed,
plus decoy-based obs/exp per group and a genome-bootstrap CI for the Asgard/TACK rate ratio.
Groups: Asgard = Asgardarchaeota; TACK = Thermoproteota + Korarchaeota; DPANN = the DPANN phyla listed below; rest = Euryarchaeota side.
uv run --project /Users/ryan/lab/work/reading-ledger python src/group_compare.py"""
import glob, numpy as np, polars as pl
pl.Config.set_tbl_rows(40); pl.Config.set_tbl_width_chars(220)
DPANN = ['Nanobdellota', 'Nanoarchaeota', 'Micrarchaeota', 'Aenigmatarchaeota', 'Altiarchaeota', 'Undinarchaeota', 'Huberarchaeota',
         'Nanohalarchaeota', 'Iainarchaeota', 'Diapherotrites', 'Pacearchaeota', 'Woesearchaeota', 'Parvarchaeota']
TW = [l.strip() for l in open('data/twins_excluded.txt')]
G = pl.read_csv('out/summary/genomes_main.tsv', separator='\t')
C = pl.concat([pl.read_parquet(f) for f in sorted(glob.glob('out/cuts1/*.parquet') + glob.glob('out/cuts2/*.parquet'))],
              how='diagonal_relaxed').filter(~pl.col('acc').is_in(TW))
K = C.filter(pl.col('placed') & (pl.col('excl') == '') & (pl.col('k') <= 3) & (pl.col('tier') == 'main'))
print('phyla:', sorted(G['phylum'].unique().to_list()))
grp = (pl.when(pl.col('phylum') == 'Asgardarchaeota').then(pl.lit('Asgard'))
       .when(pl.col('phylum').is_in(['Thermoproteota', 'Korarchaeota'])).then(pl.lit('TACK'))
       .when(pl.col('phylum').is_in(DPANN)).then(pl.lit('DPANN')).otherwise(pl.lit('Eury side')))
G = G.with_columns(grp=grp)
print(G.group_by('grp').agg(pl.col('phylum').unique().sort().str.join(',')))
rng = np.random.default_rng(1)
for nm, cf, gf in (('all', pl.lit(True), pl.lit(True)), ('no TRUB', pl.col('fam') != 'TRUB', pl.lit(True)),
                   ('no Bathyarchaeia', pl.lit(True), pl.col('class') != 'Bathyarchaeia'),
                   ('no TRUB, no Bathyarchaeia', pl.col('fam') != 'TRUB', pl.col('class') != 'Bathyarchaeia')):
    k = K.filter(cf).group_by('acc').agg(obs=pl.col('bhb').sum(), n=pl.len(), dpos=pl.col('dec_bhb').sum(), dn=pl.col('n_dec').sum())
    g = G.filter(gf).select('acc', 'grp', 'class', 'n_loci').join(k, on='acc', how='left').fill_null(0)
    T = g.group_by('grp').agg(genomes=pl.len(), loci=pl.col('n_loci').sum(), placed=pl.col('n').sum(), obs=pl.col('obs').sum(),
                              dpos=pl.col('dpos').sum(), dn=pl.col('dn').sum(), with1=(pl.col('obs') > 0).mean())
    T = T.with_columns(per1k=pl.col('obs') / pl.col('loci') * 1000, per_genome=pl.col('obs') / pl.col('genomes'),
                       exp=pl.col('placed') * (pl.col('dpos') + 0.5) / (pl.col('dn') + 1)).with_columns(
        obs_exp=pl.col('obs') / pl.col('exp')).sort('grp')
    print(f'\n== {nm}'); print(T.select('grp', 'genomes', 'loci', 'placed', 'obs', 'exp', 'obs_exp', 'per1k', 'per_genome', 'with1'))
    a = g.filter(pl.col('grp') == 'Asgard').select('obs', 'n_loci').to_numpy(); t = g.filter(pl.col('grp') == 'TACK').select('obs', 'n_loci').to_numpy()
    r = []
    for _ in range(2000):
        ia = rng.integers(0, len(a), len(a)); it = rng.integers(0, len(t), len(t))
        r.append((a[ia, 0].sum() / a[ia, 1].sum()) / max(1e-12, t[it, 0].sum() / t[it, 1].sum()))
    pt = (a[:, 0].sum() / a[:, 1].sum()) / (t[:, 0].sum() / t[:, 1].sum())
    print(f'Asgard/TACK per-locus ratio {pt:.2f} (genome bootstrap 95% {np.percentile(r, 2.5):.2f}-{np.percentile(r, 97.5):.2f})')
nt = K.filter(pl.col('bhb')).join(G.select('acc', 'grp'), on='acc').group_by('grp').agg(cuts=pl.len(), trub=(pl.col('fam') == 'TRUB').sum(), fams=pl.col('fam').n_unique())
print(nt.sort('grp'))
# per-class non-TRUB rate (classes with >= 8 genomes) in Asgard, TACK and DPANN
k = K.filter(pl.col('fam') != 'TRUB').group_by('acc').agg(obs=pl.col('bhb').sum())
g = G.select('acc', 'grp', 'class', 'n_loci').join(k, on='acc', how='left').fill_null(0)
print(g.filter(pl.col('grp') != 'Eury side').group_by('grp', 'class').agg(genomes=pl.len(), loci=pl.col('n_loci').sum(), obs=pl.col('obs').sum(),
      with1=(pl.col('obs') > 0).mean().round(2)).with_columns(per1k=(pl.col('obs') / pl.col('loci') * 1000).round(2))
      .filter(pl.col('genomes') >= 8).sort('per1k', descending=True))
