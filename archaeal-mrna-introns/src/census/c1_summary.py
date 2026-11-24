"""C1 verdict from out/c1_trna.tsv (amendments 4d, 5). C1 set = tRNAs whose introns are all canonical: a CI token and no
NCI token in tRNAscan-SE's note (pseudo/IPD tokens allowed). This is the lower-share reading of "CI only"; the
note-exactly-CI subset is printed too. GCA twins (data/twins_excluded.txt) are dropped.
uv run --project /Users/ryan/lab/work/reading-ledger python src/c1_summary.py"""
import polars as pl
TW = [l.strip() for l in open('data/twins_excluded.txt') if l.strip()]
D = pl.read_csv('out/c1_trna.tsv', separator='\t').filter(~pl.col('acc').is_in(TW)).with_columns(tok=pl.col('note').fill_null('').str.split(','))
D = D.with_columns(ci=pl.col('tok').list.eval(pl.element() == 'CI').list.any(),
                   nci=pl.col('tok').list.eval(pl.element().str.starts_with('NCI')).list.any(), pseudo=pl.col('tok').list.contains('pseudo'))


def line(nm, f):
    if f.height:
        print(f"{nm}: n={f.height} genomes={f['acc'].n_unique()} BHB={f['bhb'].mean():.3f} strict={f['strict'].mean():.3f} "
              f"decoy BHB={f['dec_bhb'].sum() / max(1, f['n_dec'].sum()):.3f}")


c = D.filter(pl.col('ci') & ~pl.col('nci'))
for nm, f in (('all', D), ('C1 SET: CI and no NCI', c), ('  of which note exactly CI', D.filter(pl.col('note') == 'CI')),
              ('  of which pseudo', c.filter(pl.col('pseudo'))), ('NCI only', D.filter(pl.col('nci') & ~pl.col('ci'))),
              ('mixed CI+NCI', D.filter(pl.col('nci') & pl.col('ci')))):
    line(nm, f)
dr = c['dec_bhb'].sum() / c['n_dec'].sum()
print('C1 verdict (bar BHB >= 0.60, decoys <= 0.08):', 'PASS' if c['bhb'].mean() >= 0.60 and dr <= 0.08 else 'FAIL')
for kk, f in (('k<=3', c.filter(pl.col('k') <= 3)), ('k>=4', c.filter(pl.col('k') >= 4))):
    line('  C1 set ' + kk, f)
