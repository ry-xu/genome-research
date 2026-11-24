import glob, polars as pl
from scipy.stats import poisson
pl.Config.set_tbl_rows(40); pl.Config.set_tbl_width_chars(200)
TW = [l.strip() for l in open('data/twins_excluded.txt')]
G = pl.read_csv('out/summary/genomes_main.tsv', separator='\t')
C = pl.concat([pl.read_parquet(f) for f in sorted(glob.glob('out/cuts1/*.parquet') + glob.glob('out/cuts2/*.parquet'))], how='diagonal_relaxed').filter(~pl.col('acc').is_in(TW))
base = C.filter(pl.col('placed') & (pl.col('excl') == '') & (pl.col('k') <= 3) & (pl.col('tier') == 'main'))
nt = base.filter(pl.col('bhb') & (pl.col('fam') != 'TRUB')).group_by('acc').agg(nt=pl.len())
pc = base.filter(pl.col('fam') != 'TRUB').group_by('acc').agg(np_=pl.len())          # all placed non-TRUB cuts (any label)
m = pl.read_csv('data/genomes_meta.tsv', separator='\t', infer_schema_length=0).select('acc', 'genus', 'family')
G = (G.join(nt, on='acc', how='left').join(pc, on='acc', how='left').fill_null(0).join(m, on='acc', how='left')
     .with_columns(src=pl.col('acc').str.slice(0, 3)))
G = G.filter(pl.col('phylum').is_in(['Thermoproteota', 'Korarchaeota', 'Asgardarchaeota']))
for lev in ('class', 'order', 'family', 'genus'):
    S = G.group_by(lev, 'src').agg(loci=pl.col('n_loci').sum(), nt=pl.col('nt').sum(), npl=pl.col('np_').sum(), g=pl.len())
    A = S.filter(pl.col('src') == 'GCA').drop('src'); F = S.filter(pl.col('src') == 'GCF').drop('src')
    J = F.join(A, on=lev, suffix='_a').filter(pl.col('loci_a') > 0)
    exp = (J['loci'] * J['nt_a'] / J['loci_a']).sum(); obs = J['nt'].sum()
    # same comparison per placed cut (BHB share), to remove any difference in how many cuts the ledger finds
    J2 = J.filter(pl.col('npl_a') > 0)
    exp2 = (J2['npl'] * J2['nt_a'] / J2['npl_a']).sum(); obs2 = J2['nt'].sum()
    print(f"{lev:7s} strata {J.height:3d}  GCF genomes {J['g'].sum():4d} loci {J['loci'].sum():7d}  per-locus: obs {obs} exp {exp:.1f} P(X<=obs) {poisson.cdf(obs, exp):.2g}"
          f"  | per-placed-cut: GCF placed {J2['npl'].sum()} obs {obs2} exp {exp2:.1f} P {poisson.cdf(obs2, exp2):.2g}")
    if lev == 'genus':
        print(J.filter(pl.col('nt_a') > 0).sort('nt_a', descending=True).head(20))
