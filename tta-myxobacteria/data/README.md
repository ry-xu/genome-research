# Data

Small tables only. Genomes are not included: NCBI accessions download on demand (`src/fetch.py`), and catalogue
genomes come from the catalogues listed in `outside_gtdb/searched.md`. Reference data: `scripts/setup_reference_data.sh`.
Everything here was produced by the code in `src/` (the post hoc checks by `src/checks/`). Genome accessions with a
`GB_`/`RS_` prefix are GTDB R232 identifiers; drop the prefix for NCBI.

## The genomes
- `genomes.tsv`: the 10 GTDB R232 genomes of the clade (8 species representatives and the only 2 other genomes of
  their genera) and the 7 genomes outside GTDB that pass the same test. Columns: source, place, completeness, Arm B
  result, the four frozen checks, TTA ends against sister orthologs, RF2 motif, tRNA-Leu(UAA), status. A dash means
  not run or not applicable.
- `relatives.tsv`: the 8 nearest relatives used as controls throughout (same columns where they apply).

## The screen
- `screen/arm_a_flags.tsv.gz`: Arm A, every flagged (genome, codon) among the 199,923 GTDB R232 representatives
  (`src/armA_run.py`).
- `screen/arm_b_candidates.tsv`: the 31 (genome, codon) pairs that passed the frozen Arm B rule, with GTDB/NCBI
  metadata (`src/candidates.py`). It includes candidates for other codons and lineages (archaeal CGA,
  Streptosporangiaceae TTA and others) that are not claimed here; see `notes.md`, "What could still be wrong".
- `prereg.md`: the logged plan. Every pass rule and check, with the time it was written, including rules added later
  and why. It also covers the other candidates.

## The four frozen checks
- `frozen_checks/check1_verdict.txt`: check 1 (another assembly of the genus passes).
- `frozen_checks/check4_verdict.txt`: check 4 (the blinded tblastn re-check on raw contigs).
- `reads/`: check 3. `<genome>.txt` is the summary (sites, depth, support, with a control codon); `<genome>.tsv` the
  per-site codon tallies. Each read pair counts once. `compare.tsv` compares this with the earlier count of both mates
  (no verdict changed; it also lists genomes not claimed here).

## Post hoc tests
- `ortholog_ends/`: the adversarial review on close-relative orthologs. `end_codon_summary.tsv` gives, per genome and
  reference (clade or sister), the codon just past the ortholog's end (`TTA_p1` of `runon`). `runon_events.tsv.gz`
  lists every such event; `interior.tsv` and `interior_tta_sites.tsv.gz` the TTAs inside genes;
  `adversarial_review.md` the reviewer's report; `stop_background_tests.out` and `context.out` the tests behind it.
- `shared_ends/`: TTA against the canonical stops at the same gene end across genera (`SUMMARY.md`; `pairs.tsv.gz`
  every resolved pair).
- `blinded_panel/`: the 16-genome panel. `blind_panel_key.tsv` maps panel ids to genomes; `blind_panel_unblind.tsv`
  and `.txt` are the blinded agent's scores joined to the key.
- `machinery/`: tRNA calls (`trnascan_panel.tsv`, tRNAscan-SE 2 with *E. coli* K-12 as control; `trna_panel.tsv`,
  ARAGORN), release factors on raw contigs (`rf_panel.tsv`), RF2 site table (`rf2_sites.out`), SPY/F201 prevalence
  across GTDB (`spy_prevalence.out`), and fixed differences across 205 shared genes (`rf_fixed_diffs.out`,
  `fixed_diff_polarised.out`, `fixed_diff_background.out`).

## Members beyond the GTDB representatives
- `gtdb_nonrep_members/`: Arm B on the two non-representative GTDB genomes of the clade's genera.
- `outside_gtdb/`: the search of NCBI and genome catalogues. `searched.md` lists every source and what was screened;
  `catalogue_members_screen.tsv` the marker-test result per candidate; `cat_members.tsv` (with `cat_members_NOTE.md`)
  the members; `armB_cat*.tsv` the frozen Arm B rule on them; `r232_polyangiaceae_genera.txt` the genera used as
  references.
