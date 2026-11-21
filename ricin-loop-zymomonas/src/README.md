# Code

Run from the folder root with `uv run python -I <script>`. Usage is in each file's first lines.

**Runs on what is here:**
- `fig_srl_rows.py OUT.png [--pfreud]`: the README figure (E. coli and *Z. mobilis*; `--pfreud` adds a
  *P. freudenreichii* row). The sequences are written into the script; where they come from is in its first lines.
- `fig_srl.py data/census OUT.png`: the census figure (`../figures/srl_census.png`).

**Runs on public downloads** (not included):
- `srl_genomes.py ZIP [ZIP...]`: exact SRL contexts in NCBI Datasets genome zips (`../data/propionibacterium/srl_counts.tsv`).
- `srl_zymo_recount.py DIR`: the second agent's *Z. mobilis* recount over a folder of `.fna.gz` assemblies
  (`../data/zymomonas/altair_recount.out`).
- `srl_silva_a1.py SILVA_138.2_LSURef_tax_silva.fasta.gz OUT.tsv`, then `srl_silva_summary.py OUT.tsv`
  (`../data/silva/`).
- `checks_pfreud/`: the independent *P. freudenreichii* check (BLAST+ `blastn` against the *E. coli* 23S, planted
  variants, short reads, missed-copy sweep). Its scripts expect the genomes in `data/dl/<accession>/` and write
  `out_pf/` in the working folder; outputs are in `../data/propionibacterium/check/`.
- `srl_readcheck.sh ACC RUN...` (SRA download, fastq-dump, minimap2) and `srl_readcheck.py`: the read check of
  metagenome-assembled carriers.

**Needs the lab's census working tree, which is not included** (the 23S of 15,739 GTDB genomes, GTDB metadata and
marker alignments; its location is set with `FROZEN_CORE`):
- `lsu_candidates.py`: the census candidate list.
- `fig_srl_numbers.py`: wrote `../data/census/loops.tsv` and `phyla.tsv`.
- `srl_covar.py`: the covariation test of marker residues (`../data/census/gaaa85.*`).
- `srl_readcheck_summary.py` (`../data/reads/summary.tsv`); `srl_readcheck.py` and `.sh` also read the census 23S.
