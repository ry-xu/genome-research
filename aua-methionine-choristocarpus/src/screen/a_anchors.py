"""GenBank K1 remainder (ATA read as Met): who supplies the M? Per genome, ATA observations (pid<95, K1-2 collapsed)
by anchor kingdom/phylum and by whether the anchor genome's declared table reads ATA as Met (echo risk).
usage: ORG_DATA=data_gb uv run python checks/k1_gb/a_anchors.py"""
import polars as pl
from Bio.Data import CodonTable
pl.Config.set_tbl_rows(80); pl.Config.set_tbl_width_chars(200)
A = ['PV826132.1', 'OR830334.1', 'CM072268.1', 'LC716139.1', 'LC716140.1']
RES = 'ACDEFGHIKLMNPQRSTVWY'
x = pl.scan_parquet('runs/gb_v1_all/cols/*.parquet').filter(pl.col('acc').is_in(A) & (pl.col('codon') == 34)).collect()
L = pl.read_parquet('runs/gb_v1_all/loci/*.parquet', columns=['acc', 'locus', 'pid'])
x = x.join(L, on=['acc', 'locus'], how='left').filter(pl.col('pid') < 95).unique(['acc', 'sid', 'spos', 'codon', 'res'])
x = x.with_columns(pl.col('res').map_elements(lambda r: RES[r], return_dtype=pl.String).alias('aa'))
T = pl.read_csv('data/anchors_tax.tsv', separator='\t', schema_overrides={'sid': pl.Int32}).select('sid', pl.col('acc').alias('aacc'), 'kind', 'phylum', 'kingdom')
c = pl.concat([pl.read_parquet(f'data/cds_{k}.parquet', columns=['acc', 'transl_table']) for k in ('mito', 'plastid')]).filter(pl.col('transl_table') != '')
tt = c.group_by('acc', 'transl_table').len().sort('len', descending=True).group_by('acc').first().select(pl.col('acc').alias('aacc'), pl.col('transl_table').cast(pl.Int64).alias('att'))
ataM = {t for t in CodonTable.unambiguous_dna_by_id if CodonTable.unambiguous_dna_by_id[t].forward_table.get('ATA') == 'M'}
x = x.join(T, on='sid', how='left').join(tt, on='aacc', how='left').with_columns(
    pl.when(pl.col('aacc').is_null()).then(pl.lit('prok')).when(pl.col('att').is_in(list(ataM))).then(pl.lit('ATA=M table')).otherwise(pl.lit('ATA=I table')).alias('anchor_code'))
print(x.group_by('acc', 'anchor_code', 'aa').len().filter(pl.col('aa').is_in(['M', 'I', 'L', 'V'])).sort('acc', 'anchor_code', 'aa'))
print(x.filter(pl.col('aa') == 'M').group_by('acc', 'kind', 'kingdom', 'phylum').len().sort('acc', 'len', descending=[False, True]))
