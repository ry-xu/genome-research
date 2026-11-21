# Code

Run everything from the folder root with `uv run python -I <script>`.

**Runs on what is here:**
- `check_numbers.py`: recounts the README's numbers (286 of 290, the per-genome table, the other 17 *Wickerhamiella*,
  the Fisher tests) from `../data/census/`.
- `fig_proline_rows.py OUT.png`: the figure in the README. The 27-nt stretches are written into the script; each
  occurs once in `../data/genomes/all.fa` at the first base given.
- `followup/trnas.py OUT.tsv data/genomes/all.fa`: the tRNA calls of the follow-up (reproduces
  `../data/followup/trna/trnas.tsv`); `followup/origin.py`, `stems.py`, `branch_changes.py` and `t0b.py` read that
  table (usage in each file).

**Needs the census, which is not included** (about 22,700 genomes; three files, about 215 MB, built by the lab's
genome-fetch and profile-building steps, which are not included either):
- `codes_census/census_org.py` (with `../scan_org.py` and `../codes.py`): the census: translate each genome with its
  declared code, align each protein to the profile models, record the codon at every match column.
- `codes_census/sitescan_org.py`: the census loader and the column consensus (strict = at least 90% of the
  family-weighted residues agree).
- `cun_pro/cun_by_res.py`: the 286-of-290 count. `cun_pro/export_sites.py`: wrote `../data/census/cun_sites.tsv` and
  `strict_pro_codons.tsv` for this folder.
- `codes_census/pro_leu.py` (`pro_leu_wss.tsv`), `codon_contrast.py` (`cun_CTA.tsv`, `cun_CTT.tsv`),
  `reads_codons.py` (`../data/reads/`; also needs minimap2, samtools and the read runs), `pro_grid.py`
  (`../figures/pro_grid_tree.png`), `org_trnas.py` (tRNA calls).
- `tree/`: `concat.py` (the concatenated protein alignment, CUN masked), then IQ-TREE; `support.py` and `show.py` read
  the tree.
- `followup/census_genes.py`, `map_run.sh`, `r2_batch.sh`, `rdd.py`, `species_check.py`: the RNA-editing arm (gene
  tables from the census inputs, which it expects under `$NEW_STOPS`; reads streamed from SRA with fastq-dump into
  minimap2; the RNA-DNA difference caller). `followup/audit_rna/`: the audit's scripts for that arm.

**Kept as run, by the checking agents** (paths refer to the lab's working tree):
- `audit/`: the code audit's scripts (`indep.py` is its reproduction without profile models; `census_set_reads.py` and
  `iso_census.py` made `../data/reads/*.out`).
- `blind_recheck/`: the blinded re-check's own pipeline (gene finding, reference alignments, codon tallies, report).
