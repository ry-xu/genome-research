# Data

Small summary tables and audit reports only. The genomes (GTDB release 232 / NCBI), raw sequencing runs (SRA), and
large BAM/SAM alignments are public or too large and are not included. Reports under `checks/` and `asgard_b35/` were
copied from the lab record; local paths inside them refer to the original working folders.

## `census/` (3,666-genome archaeal DNA scan)
- `genomes_main.tsv`, `genomes_wide.tsv`: per-genome placed cuts, BHB-motif calls, and decoy-expected counts across
  3,666 archaeal genomes.
- `classes_main.tsv`, `classes_wide.tsv`: class-level enrichment tests (Bonferroni-corrected).
- `sites_main.tsv`, `sites_wide.tsv`: reference-protein column clusters across archaeal classes.
- `summary.txt`, `group_compare.txt`: supergroup summaries (Asgard, TACK, DPANN, Euryarchaeota side) with and without
  TRUB (*cbf5*) and Bathyarchaeia.
- `RULES.md`: frozen rules and amendments written before scoring.

## `asgard_b35/` (*Ca.* Lokiarchaeum ossiferum B35, `L02` and `S859`)
- `rna_mol_templates.tsv`: deduplicated RNA molecules spanning `L02` (385 lacking the 31 bases, 26 read-through) and
  control sites from SRR27313660 and SRR32148620 (input to `../src/fig_l02_public.py`).
- `rna_mol_gaps.tsv`, `rna_mol_sites.tsv`: per-gap and per-site molecule summaries (`383` exact 31-base matches).
- `fig_l02_dna_templates.tsv`, `dna_sites.tsv`: DNA fragments spanning `L02` from SRR22273278 (0 lacking, 46 through;
  input to `../src/fig_l02_public.py`).
- `c1.out`: PCR-skipping duplicate-family test (0 mixed families vs 13.47 expected, Poisson P = 1.41e-6).
- `c2.out`: DNA-hairpin skipping scan across 3,735 loci folding at least as tightly as `L02` (0 lose >=10%).
- `INDEPENDENT_agent.md`: independent kill-check agent report (aligner-free k-mer counts, strand balance, B-36 check).

## `tack_c0/` (*Ca.* Methanosuratincola verstraetei LCB70, `c0`)
- `f1.events.tsv`, `f1.txt`: RNA template counts across studies S1–S3 (`GCA_049570475.1:c0`: 45 joined, 6 through).
- `c0_reads.tsv`: per-read details at `c0` across the 9 PRJNA916083 metatranscriptomes.
- `metag_words.tsv`: metagenome DNA check (SRR23044201: 0 joined, 58 through).
- `gap_census_s2.txt`, `hairpin_s2_W22.txt`, `hairpin_s2_W22_loci.tsv`: background gap and 22-nt RNA hairpin scan
  (84 loci with dG <= -5.6 kcal/mol, 21,558 spanning reads, 0 removal reads).
- `SURVEY62.md`: survey of all 62 Methanomethylicia species representatives (25 copies in 9 Methanomethylicaceae
  genera).
- `c0_homologs2.out`, `covar.out`, `tree.out`, `decoy_split.out`: BHB fold check across all 25 copies (25/25 paired
  central helix at exact ends), compensatory inner-stem covariation, parsimony tree count (>=4 independent double
  changes), and 112-pair decoy check (0/112).

## `checks/`
Independent referee, kill-check, figure fact-check, and prior-art reports (`sabik_asgard_CHECK.md`,
`FACTCHECK_fig_l02_agent.md`, `killcheck_c0_REPORT.md`, `killcheck_a5_REPORT.md`, `PRIOR_ART_asgard.md`,
`PRIOR_ART_c0.md`).
