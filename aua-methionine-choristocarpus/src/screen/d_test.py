"""ATA-as-Met test under checks/k1_gb/RULES_ata.md (fixed before this code).
usage: uv run python checks/k1_gb/d_test.py  -> checks/k1_gb/ata_test.parquet (+ printed flags and controls)"""
import os, sys
import polars as pl
from Bio.Data import CodonTable
sys.path.insert(0, 'src')
import read_k1  # noqa: E402

RUNS = [('runs/refseq_v1_mito', 'data'), ('runs/refseq_v1_plastid', 'data'), ('runs/gb_v1_all', 'data_gb')]
ATA_M = sorted(t for t, c in CodonTable.unambiguous_dna_by_id.items() if c.forward_table.get('ATA') == 'M')
COD = {32: 'ATT', 33: 'ATC', 34: 'ATA', 35: 'ATG'}


def declared(d):
    read_k1.REC = d
    return read_k1.declared()


anc = pl.read_csv('data/anchors_tax.tsv', separator='\t', columns=['sid', 'acc'], schema_overrides={'sid': pl.Int32})
anc = anc.join(declared('data').select('acc', 'tt'), on='acc', how='left')
bad_sid = anc.filter(pl.col('tt').is_in(ATA_M))['sid']          # anchors whose table reads ATA as Met; others (incl. prok) kept
out = []
for run, d in RUNS:
    for f in sorted(os.listdir(os.path.join(run, 'cols'))):
        c = pl.scan_parquet(os.path.join(run, 'cols', f)).filter(pl.col('codon').is_in(list(COD)) & pl.col('res').is_in([7, 10])
                                                                  & (pl.col('spos') > 30) & ~pl.col('sid').is_in(bad_sid.implode()))
        L = pl.scan_parquet(os.path.join(run, 'loci', f)).select('acc', 'locus', 'pid')
        x = c.join(L, on=['acc', 'locus'], how='left').filter(pl.col('pid') < 95).unique(['acc', 'sid', 'spos', 'codon', 'res'])
        out.append(x.group_by('acc').agg(*[((pl.col('codon') == k) & (pl.col('res') == r)).sum().alias(f'{n}_{a}')
                                           for k, n in COD.items() for r, a in [(10, 'M'), (7, 'I')]]).collect()
                   .with_columns(pl.lit(d).alias('d')))
t = pl.concat(out).unique('acc', keep='first')
D = pl.concat([declared(d).with_columns(pl.lit(d).alias('d')) for d in ('data', 'data_gb')])
t = t.join(D, on=['acc', 'd'], how='left').with_columns(pl.col('tt').fill_null(1), pl.col('tt_src').fill_null('default'))
R = pl.concat([pl.read_parquet(f'{d}/records_ranked.parquet', columns=['acc', 'kind', 'organism', 'order', 'class', 'phylum'])
               for d in ('data', 'data_gb')]).unique('acc')
t = t.join(R, on='acc', how='left').with_columns(pl.col('tt').is_in(ATA_M).alias('decl_ATA_M'))


def flag(cod, ctl='ATT'):
    m, i = pl.col(f'{cod}_M'), pl.col(f'{cod}_I')
    ok = (m >= 10) & (m >= 2 * i)
    if ctl:
        ok = ok & (pl.col(f'{ctl}_M') <= 0.05 * (pl.col(f'{ctl}_M') + pl.col(f'{ctl}_I')))
    return ok


t = t.with_columns(flag('ATA').alias('flag_ATA'), flag('ATT', None).alias('flag_ATT'), flag('ATC', 'ATT').alias('flag_ATC'))
t.write_parquet('checks/k1_gb/ata_test.parquet')
pl.Config.set_tbl_rows(100); pl.Config.set_tbl_width_chars(220); pl.Config.set_fmt_str_lengths(30)
pop = t.filter(~pl.col('decl_ATA_M'))
print('ATA=Met tables:', ATA_M, '| anchors dropped as echo risk:', len(bad_sid))
print('population (declared ATA=Ile):', pop.height, '| flag_ATA', pop['flag_ATA'].sum(), '| flag_ATT', pop['flag_ATT'].sum(), '| flag_ATC', pop['flag_ATC'].sum())
pos = t.filter(pl.col('decl_ATA_M'))
ev = pos.filter((pl.col('ATA_M') + pl.col('ATA_I')) >= 10)
print('recall: declared ATA=Met genomes', pos.height, '| with >=10 ATA at M/I', ev.height, '| flagged', ev['flag_ATA'].sum())
F = pop.filter('flag_ATA').select('acc', 'kind', 'tt', 'tt_src', 'organism', 'order', 'class', 'ATA_M', 'ATA_I', 'ATT_M', 'ATT_I', 'ATG_M', 'ATG_I')
print(F.sort('kind', 'class', 'order', 'organism'))
print(F.group_by('kind', 'class', 'order').agg(pl.len().alias('records'), pl.col('organism').n_unique().alias('organisms')).sort('kind', 'class'))
print('flagged ATT/ATC genomes:'); print(pop.filter(pl.col('flag_ATT') | pl.col('flag_ATC')).select('acc', 'organism', 'order', 'ATT_M', 'ATT_I', 'ATC_M', 'ATC_I'))
