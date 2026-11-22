# Code

Run everything from the folder root with `uv run python -I <script>`.

**Runs on what is here:**
- `check_numbers.py`: recounts the README's numbers (32 genomes, 124 sites, 28 of 32, 0 of 4,322, the reads and the
  protein tallies) from `../data/`, and prints the kept outputs of the census-based tests.
- `aga_ala/fig_aga.py data/census/trna_vs_usage_fungi.tsv OUT.png`: the figure in the README (reproduces it exactly).
- `aga_ala/trna_vs_usage.py contrast.tsv aragorn.tsv OUT.tsv`: the tRNA join; with `../data/census/aga_ala_Fungi.tsv`
  and `../data/trna/fungi_ucu.aragorn.tsv` it reproduces `trna_vs_usage_fungi.tsv`.
- `aga_ala/trna_origin.py data/trna/basidio.aragorn.tsv`: each tRNA(UCU)'s identity to the genome's other tRNAs.

**Needs the census, which is not included** (about 22,700 genomes; built by the lab's genome-fetch and
profile-building steps, which are not included either):
- `codes_census/census_org.py` (with `../scan_org.py` and `../codes.py`): the census: translate each genome with its
  declared code, align each protein to the profile models, record the codon at every match column.
- `codes_census/sitescan_org.py`: the census loader and the column consensus (strict = at least 90% of the
  family-weighted residues agree).
- `codes_census/codon_contrast.py` (`../data/census/aga_ala_*.tsv`), `codes_census/org_trnas.py` (tRNA calls),
  `codes_census/reads_codons.py` (`../data/reads/`; also needs minimap2, samtools and the read runs).
- `aga_ala/lineage_nosel.py` (`lineage_nosel.txt`), `aga_ala/all_columns.py` (`all_columns.txt`).

**Protein search** (`prot/`; needs Sage v0.14.7, ThermoRawFileParser and the .raw files from iProX):
`build_db.py` (the three-version database), `fetch_convert.sh`, `one_tmt.sh`, `stream_tmt2.sh` (download and convert),
`run_sage.sh`, `sage_cfg.py` (the configuration), `tally_aga.py` (the tallies), `reporters.py` (TMT reporters),
`pep_in_genome.py` (the peptide in six-frame genome translations).

**Kept as run, by the checking agents** (paths refer to the lab's working tree or the agents' temporary folders):
- `audit/`: the code audit's scripts (`repro_star.py` is its reproduction without the census; `distinct_cols.py` the
  pseudoreplication check; `shared_cols.py` the shared-column test, P ≈ 7e-4).
- `blind_recheck/`: the blinded re-check's pipeline (Codetta runs, relatives, ortholog-column tallies).
- `prepost_check/`: the pre-post checker's scripts (`analyze.py` and `extra*.py` give the 78% / 97.5% and the
  non-strict 23% / 1%).
