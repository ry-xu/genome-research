# Code

Python, run with `uv run python -I <script>` from the folder root. The scripts find each other through relative
paths (`codes.py` here, the census modules in `codes_census/`).

**Runs on what is here:**
- `check_numbers.py`: recounts the README's numbers (23 and 0 positions, 179 of 184, 564 of 576, 7,109 of 7,298
  reads, 3 CGA, 50 of 50) from `../data/`.
- `fig_uga_arg_rows.py OUT.png`: the figure in the README, from `../data/genomes/fig_genomes.fa`.
- `tga_arg/fig_pneumo.py ../data/sites/pneumo_cols.tsv OUT.png`: the figure in `../notes.md`.

**Needs the census, which is not included:**
- `codes.py`: genetic code tables and the translation used throughout.
- `codes_census/`: the mitochondrial census.
  - `scan_org.py`, `census_org.py`: translate each genome's protein genes with its declared code and align them to
    profile HMMs (pyhmmer), storing every placed codon with its profile column (`*.census.npz`). The profiles and
    the gene tables (from RefSeq and GenBank organelle records) are not included.
  - `sitescan_org.py`: per-column consensus, family-weighted, and the strict (≥ 90%) columns.
  - `profile_scan.py`: the all-codons scan that flagged UGA in *P. jirovecii*.
  - `reads_codons.py`: map RNA reads (minimap2) to a genome's genes and tally the codon at each placed site.
- `tga_arg/`: the *Pneumocystis* analysis (`pneumo_cols.py`, `pneumo_same_cols.py`, `isolates.py`, `tally_rna.py`,
  `fig_pneumo.py`, plus `clade_tga.py` and `tga_by_gene.py`, used while exploring). Each script's docstring gives its
  usage.
- `pneumo_cun.py`: the CUN = Thr check.

Memory: `pneumo_cols.py` on the full census needs more than 6 GB.
