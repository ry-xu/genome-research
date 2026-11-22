"""For the ATA->Met genomes: residue profile of ATT/ATC/ATA/ATG observations (pid<95, collapsed) and codon usage
in observed positions. If ATA were Ile, it should sit at I columns like ATT/ATC do.
usage: ORG_DATA=data_gb uv run python checks/k1_gb/b_codons.py"""
import polars as pl
pl.Config.set_tbl_rows(80); pl.Config.set_tbl_width_chars(200)
A = ['PV826132.1', 'OR830334.1', 'CM072268.1', 'LC716139.1', 'LC716140.1']
C = {'ATT': 32, 'ATC': 33, 'ATA': 34, 'ATG': 35}
RES = 'ACDEFGHIKLMNPQRSTVWY'
x = pl.scan_parquet('runs/gb_v1_all/cols/*.parquet').filter(pl.col('acc').is_in(A) & pl.col('codon').is_in(list(C.values()))).collect()
L = pl.read_parquet('runs/gb_v1_all/loci/*.parquet', columns=['acc', 'locus', 'pid'])
x = x.join(L, on=['acc', 'locus'], how='left').filter(pl.col('pid') < 95).unique(['acc', 'sid', 'spos', 'codon', 'res'])
inv = {v: k for k, v in C.items()}
x = x.with_columns(pl.col('codon').map_elements(lambda c: inv[c], return_dtype=pl.String).alias('cod'),
                   pl.col('res').map_elements(lambda r: RES[r], return_dtype=pl.String).alias('aa'))
t = x.group_by('acc', 'cod').agg(pl.len().alias('N'), (pl.col('aa') == 'M').sum().alias('M'), (pl.col('aa') == 'I').sum().alias('I'),
                                  pl.col('aa').is_in(['L', 'V']).sum().alias('LV')).with_columns(
    (pl.col('M') / pl.col('N')).round(2).alias('Mshare'), (pl.col('I') / pl.col('N')).round(2).alias('Ishare'))
print(t.sort('acc', 'cod'))
