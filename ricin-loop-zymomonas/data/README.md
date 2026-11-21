# Data

Small tables only. The genomes (GTDB release 232, NCBI Datasets zips), the SILVA 138.2 LSU Ref FASTA and the read runs
are public and not included. Times are US Eastern, November 2026.

## `census/` (GTDB 23S census, 15,739 bacterial genomes)
- `recount.json`, `recount.out`: a fresh check agent's re-read of every genome the pipeline listed for a 2661 change,
  with its own code (the lab calls it the "kill" agent). 116 listed; 10 dropped (a long insert inside the loop); 102
  keep G2661A in every full-length copy, and 85 carry the full GAAA loop. The final carrier list is `recount.json`.
- `loops.tsv`, `phyla.tsv`: the numbers behind `../figures/srl_census.png` (`src/fig_srl_numbers.py`, drawn by
  `src/fig_srl.py`: `uv run python -I src/fig_srl.py data/census OUT.png`).
- `gcaa17.tsv`: the 17 genomes that read GCAA in every full-length copy (in both the `G2661A` and `A2660C` lists of
  `recount.json`), with their GTDB taxonomy: 6 phyla.
- `gaaa85.tsv`, `gaaa85.out`: the covariation test of GTDB marker residues (IF2, EF4, BipA and the rest of the marker
  set) against GAAA origins (`src/srl_covar.py`). 0 of 149 columns pass.

## `zymomonas/`
- `altair_recount.out`: every complete *Z. mobilis* assembly in GenBank, recounted by a second agent with its own
  code (`src/srl_zymo_recount.py`; the `loose` count, AGTAC NNNN GGACCG on both strands). 21 rows, because three
  assemblies appear in two versions: 18 assemblies, 2 or 3 copies each, all GAAA. The `strict` pattern adds the
  *E. coli* bases upstream, which *Z. mobilis* does not share, so it finds nothing.

## `propionibacterium/`
- `srl_counts.tsv` (no header: assembly, matches, loops): per assembly, the number of exact AGTAC NNNN GGACC matches
  and the loop in each (`src/srl_genomes.py`): the 55 complete *P. freudenreichii* assemblies (two are a phage and a
  fragment, with no 23S; the other 53 are all GCAA, 120 copies) and the *E. coli* K-12 (GCF_000005845.2, 7 GAGA) and
  *Cutibacterium acnes* (GCF_000008345.1, 3 GAGA) controls. Some rows have extra matches outside the 23S copies
  (ATGG in the type strain GCA_000940845.1, CTCA in seven others, and three non-GAGA matches in each control); the
  independent BLAST-based check counts only the 23S copies and agrees (`check/tally_out.txt`). The sister species are
  in the independent check only (`../checks/pfreudenreichii_independent.md`).
- `check/`: printed tallies of the independent check (`src/checks_pfreud/`): `tally_out.txt` (copies and loops),
  `strains_out.txt` (duplicates and strains), `meta_out.txt` (submitters, platforms, type material),
  `sra_loop_counts.txt` (short reads).

## `reads/`
- `summary.tsv`: the read-level check of 3 metagenome-assembled GAAA carriers (and the closing-pair carriers of a
  separate result), one row per genome (`src/srl_readcheck.sh`, `srl_readcheck.py`, `srl_readcheck_summary.py`).

## `silva/`
- `summary_Ref.txt`: SILVA 138.2 LSU Ref read with an exact AGTAC NNNN GGACC match in the last 30% of each sequence
  (`src/srl_silva_a1.py`, `src/srl_silva_summary.py`): controls, bacterial loop counts (organelles excluded), non-GAGA
  genera, organelles. Archaea read 0 because their stem differs (GGAAC); the rule is bacterial.
- `non_gaga_bacteria_LSURef.tsv`: every bacterial entry not reading GAGA (or matching more than once).

## `checks/`
- `priorart_crw_silva.md`: the second prior-art search (CRW diagrams, frequency file and 2001 alignment, RiboVision,
  Blanchet 2014 Table 1, Kudrin 2018).
- `pfreudenreichii_independent.md`: the independent *P. freudenreichii* check. Verdict: NARROWED (real and
  species-wide; the sequence is an unnamed row in Blanchet 2014 and the change is Kudrin's mutant).
