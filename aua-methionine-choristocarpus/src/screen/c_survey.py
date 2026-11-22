"""All genomes, three runs: ATA and ATT observations at anchor M vs I (pid<95, K1-2 collapsed), with group/order and
declared table. Finds ATA-behaves-as-Met genomes below K1's N>=10 too. Per batch (memory).
usage: uv run python checks/k1_gb/c_survey.py  -> checks/k1_gb/ata_survey.parquet"""
import os, sys
import polars as pl
sys.path.insert(0, 'src')
RUNS = [('runs/refseq_v1_mito', 'data', 'refseq'), ('runs/refseq_v1_plastid', 'data', 'refseq'), ('runs/gb_v1_all', 'data_gb', 'gb')]
out = []
for run, rec, src in RUNS:
    for f in sorted(os.listdir(os.path.join(run, 'cols'))):
        c = pl.scan_parquet(os.path.join(run, 'cols', f)).filter(pl.col('codon').is_in([32, 34, 35]) & pl.col('res').is_in([7, 10]))
        L = pl.scan_parquet(os.path.join(run, 'loci', f)).select('acc', 'locus', 'pid')
        x = c.join(L, on=['acc', 'locus'], how='left').filter(pl.col('pid') < 95).unique(['acc', 'sid', 'spos', 'codon', 'res'])
        out.append(x.group_by('acc').agg(*[((pl.col('codon') == k) & (pl.col('res') == r)).sum().alias(f'{n}_{a}')
                                           for k, n in [(32, 'ATT'), (34, 'ATA'), (35, 'ATG')] for r, a in [(10, 'M'), (7, 'I')]]).collect()
                   .with_columns(pl.lit(src).alias('src'), pl.lit(os.path.basename(run)).alias('run')))
t = pl.concat(out).unique('acc', keep='first')
R = pl.concat([pl.read_parquet(f'{d}/records_ranked.parquet', columns=['acc', 'kind', 'organism', 'order', 'class', 'phylum', 'group', 'length'])
               for d in ('data', 'data_gb')]).unique('acc')
t = t.join(R, on='acc', how='left')
t.write_parquet('checks/k1_gb/ata_survey.parquet')
print(t.shape)
