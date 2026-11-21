# Data

Small tables only. The census files (about 22,700 mitochondrial genomes, about 215 MB), the profile models and the read
runs are not included. Accessions are GenBank (BK07xxxx are third-party assemblies from Sarhan et al. 2025). Times
are US Eastern, November 2026. The four genomes: *W. infanticola* BK074532.1, *W. siamensis* BK074538.1, *W.
verensis* BK074542.1, *S. pinicola* BK074262.1.

## `census/`
- `cun_sites.tsv`: every aligned CUU/CUA codon of the four, with its profile column, the column's consensus residue,
  its weighted share and whether the column is strict (share at least 0.9). 286 of the 290 strict rows are proline
  (`src/cun_pro/export_sites.py`, same loader as `src/cun_pro/cun_by_res.py`).
- `strict_pro_codons.tsv`: per genome, all codons at strict proline columns (91–93): CUU/CUA, CCN, other.
- `pro_leu_wss.tsv`: per genome of the order (54 genomes): codons at strict proline columns (`P:CCN`, `P:CTN`) and
  strict leucine columns (`L:CTN`, `L:TTR` = UUA/UUG), and all aligned CCN and CTN (`src/codes_census/pro_leu.py`).
  The `genus_species` column holds the genus.
- `cun_CTA.tsv`, `cun_CTT.tsv`: the codon contrast over 151 Dipodascomycetes (`src/codes_census/codon_contrast.py`):
  the codon at strict proline and leucine columns, its leucine synonyms at the same, and the one-sided Fisher p. The
  last row pools all genomes.
- `pro_grid_tree.tsv`: per genome in `../figures/pro_grid_tree.png`, the codons drawn (CCN, CUN, other, no site).

## `reads/`
Codon tallies from reads at each CUU/CUA and CCN site (`src/codes_census/reads_codons.py`; columns: codon, profile
column, consensus, reads, reads matching the assembly, other triplets).
- `BK0*.tsv`: the source run of each assembly (mapping quality 0 allowed). These tables were labelled with a consensus
  from one of the three census files only, which gives 77/81/81/43 CUU/CUA-at-proline sites instead of 78/82/82/44;
  `census_set_reads.out` is the recount on the census site set (the numbers in the README).
- `mapq20/`: the same at mapping quality 20 or more.
- `other_sequencings/`: *W. siamensis* PYCC 7069 (ERR11532760, another lab, 2024), recounted on the census set in
  `iso_census.out`; DS02 PacBio subreads (SRR4787334, first 26,351 subreads) against BK074532.

## `tree/`
IQ-TREE runs on 49 Dipodascales mitochondrial genomes (`accs.txt`, names in `names.fa`; concatenated proteins, CUN
masked; `src/tree/concat.py`): `wt.*` (full search, 1,000 ultrafast bootstraps) and `wtfast.*` (fast search, SH-aLRT).
`support_*.txt` list the support of the nodes joining the four to their relatives (`src/tree/support.py`).

## `genomes/`
- `all.fa`: 57 mitochondrial genomes from GenBank (the 49 of the tree and 8 controls, among them *S. cerevisiae* and *Ashbya gossypii*), used by the follow-up tRNA arm
  and to check the figure's sequences.

## `checks/`
- `audit_report.md`: the code audit (fresh agent), with an independent reproduction without profile models.
- `blind_recheck/report.md`: the blinded re-check of eight genomes.
- `factcheck.md`, `prepost_check.md`, `postcheck_additions.md`: fact-checks of the write-up by fresh agents, in order.
- `tweet_check.md`: line-by-line check of the thread and its figure (21 Nov).
- `ds02_identity/`: DS02's mitochondrial contigs against the *Wickerhamiella* genomes (minimap2 PAF), and DS02's NCBI
  assembly report.
- `OR703884_vs_BK074532.paf`: Wolters et al.'s *W. infanticola* assembly against BK074532 (identical after rotation).

These reports were written in the lab's working tree; paths in them (`results/cun_pro/`, `work/...`) refer to it.
`results/cun_pro/` there holds what `census/`, `reads/` and `checks/` hold here.

## `followup/` (20 Nov; preliminary)
- `RULES.md`: the rules, frozen before any target analysis (with two amendments). `RESULT.md`: the results.
- `trna/`: ARAGORN tRNA calls for all genomes (`trnas.tsv`), the origin statistics for the targets and the controls
  (`target_T.txt`, `controls_T.txt`), acceptor-stem bases (`stems.txt`), changes on the branch (`branch_changes.txt`),
  and the audit (`audit_REPORT.md`).
- `rna_editing/`: the search for RNA and protein data (`R0_summary.md`, `sra_strain_search.txt`, `pride.txt`,
  `rna_runs_strain.tsv`, `bk_sources.txt`), the per-run results (`r2_*.txt`), the RNA-DNA differences in
  *W. versatilis* SRR30194995 (`SRR30194995.rdd.tsv`) and the *Arabidopsis* control (`ath_ERR15942850.rdd.tsv`,
  `pcr_ath.txt`, `pcr_known_edits.txt`), coverage (`coverage_*.txt`), and the same-strain Nanopore DNA filter
  (`rdd_with_ont_filter.txt`, `ont_check.txt`).
