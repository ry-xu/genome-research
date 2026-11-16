# Data

Small tables only. Genomes, reads, Swiss-Prot and the DIAMOND database are not included; `scripts/setup_reference_data.sh`
says how to get them. KACI's and Codetta's tables are cited, not copied ([`../notes.md`](../notes.md), References).
Everything here was produced by the code in `../src/`. GTDB accessions with a `GB_`/`RS_` prefix are GTDB R232
identifiers; drop the prefix for NCBI. Times in the rule files are US Eastern, November 2026.

- `prereg.md`: the logged plan. The frameshift scan's gate, the CGG test fixed before the TANB77 scan, and the log of
  11 Nov 15:39 to 16:17 (the discovery, the test's results, the tRNA census, the turnover and CAG-698 tests, the
  blinded re-check), with the time each entry was written. It points to the other rule files below.

## Genome sets (`genome_sets/`)
GTDB R232 metadata rows (taxonomy, CheckM2 completeness and contamination, contigs, GC, translation table) for each set
scanned. `local` is the file in the genome cache used for the run (see `../src/paths.py`).
- `tanb77_set.tsv`: test set T, one genome per TANB77 genus (64), plus GCA_944393025.1 (a second genome of the
  clearest species, used for the read check).
- `clos_lowgc_set.tsv`: the 50 low-GC Clostridia family representatives outside TANB77 (the 6 frozen controls and the
  extension).
- `cag698_set.tsv`: the 22 usable CAG-698 species representatives (second-lineage test).
- `moll_set.tsv`: the Mollicutes family representatives.
- `tanb77_reader_set.tsv`: the two CAG-508 genomes that keep a CGG tRNA (CAG-793, CAG-452), scanned afterwards.

## The frozen test (`frozen_test/`)
- `cgg_test_T.tsv`: per genome of set T (`src/analysis/cgg_test.py`): aligned loci, clean +1 shifts (`p1`), those at
  CGG (`cgg`), their share, the call (1 = skip), ARAGORN tRNA-CCG count (`ccg`) and tRNA count.
- `cgg_test_C.tsv` (the 6 frozen controls), `cgg_test_C50.tsv` (all 50), `cgg_test_cag698.tsv` and
  `cgg_test_cag698r.tsv` (CAG-698; the second with tRNA status from the screen), `cgg_test_moll.tsv` (Mollicutes): the
  same columns.
- `inframe_{tanb77,cag698,clos,moll}.tsv`: in-frame use of CGG and the other arginine codons per 1,000 codons, in hits
  that align without a shift (`src/analysis/inframe_cgg.py`).
- `shift_T_pooled.tsv`: the related protein's residue opposite each skipped CGG, pooled over set T (1,458 rows; Gly at
  1,237). `shift_residue_445.tsv`: the same for the discovery genome, with the other breaks
  (`src/analysis/shift_residue.py`).
- `figure_data.tsv`: the numbers behind `../figures/skip_vs_reader.png` (`src/analysis/cgg_figure.py`).

## Per-site frameshift calls (`scans/`)
Output of `src/fs_scan.py`. `<genome>.fs.tsv.gz` has one row per frame break (position, type, flank identity,
homopolymer run, distance to the old frame's stop, sequence context, homolog residue); `src/summarize.py` `flag()`
classes them (clean, paired, programmed). `<genome>.hits.tsv.gz` has one row per aligned locus.
- `tanb77/`: set T (breaks and hits); `tanb77_reader/`: the two reader genomes; `discovery/`: GCA_963931445.1 from the
  family scan.
- `controls_clos_lowgc/`, `cag698/`, `mollicutes/`: breaks only, for the control sets.
- `pig_check/`: the pig-gut genome MGYG000302328 (RNA test, Set B).

## tRNA (`trna/`)
- `cgg_reader_order.tsv`, `cgg_reader_family.tsv`: per GTDB order and family, usable species representatives and how
  many have no CGG-reading tRNA (CCG or TCG anticodon), from an ARAGORN screen of the GTDB R232 representatives
  (`src/analysis/cgg_reader_census.py`; usability rule in its docstring). `cgg_reader_genomes.tsv`: the usable genomes
  with no reader.
- `pig_trnascan/`: tRNAscan-SE 2 calls (`.out`) and secondary structures (`.ss`) for the 24 pig-gut TANB77 genomes,
  and the planted 8-nt-loop tRNA used as a positive control (`control_8loop.*`; `src/trna_pig/`).
- The blinded re-check's own tRNA search is in `blinded_recheck/trna/`.

## Raw reads of the discovery genome (`reads/`)
- `allele_cgg.tsv`: Illumina run ERR11860625, reads carrying the with-C and without-C 25-mers at each of the 21 CGG
  skip sites (140 with, 0 without; site 18 has almost no reads). `allele_other.tsv`: the same for the other breaks.
- `cgg.tsv`, `other.tsv`, `breaks.tsv`: the sites; `kmers_*.tsv`: the 25-mers (`src/analysis/allele_kmers.py`);
  `read_check.tsv`: the stage-3 read check per site.

## Turnover (`turnover/`)
- `site_turnover.tsv.gz`: one row per (window, other genome): skip sites, GGN starts and other codon starts of each set
  T genome aligned (blastn dc-megablast) to every other genome, with relation (species, genus, family, other), the
  base opposite the C and net length change (`src/analysis/site_turnover.py`).
- `site_turnover_prfb.tsv`: the same for the *prfB* shift site (positive control for a conserved site).

## The blinded re-check (`blinded_recheck/`)
- `VERDICT.md`: the agent's verdict, written before unblinding (genomes G1-G9; the key is the table in
  `prereg.md`, entry 16:17).
- `out/`: per genome, miniprot CGG loci (`cgg_loci`), frame events (`events`), shift sites at CGG with the residue
  opposite (`fs_cgg_sites`), and summaries; `plant_control.txt` and `G6planted.*`/`G6plantR.*` are the planted-site
  controls; `rf2_control.tsv` the *prfB* control. Alignments (`.paf`) are not included. Code: `src/recheck/`.

## RNA test (`rna/`)
- `RULES_AND_LOG.md`: the rule (12 Nov 07:00, before any RNA read was mapped), Amendments 1-5 with times, and the log.
  Amendment 5 (the skip / read-through split) was made after the count.
- `sites.tsv.gz`: every site of the test: Set A (the C of each in-frame CGG met by a skip walk), Set C (CGA/CGC/CGT
  controls), A0 and B, with contig position, strand and CDS (`src/rna_cgg.py sites`). `pig_genomes.tsv`: the 24
  genomes, their tRNA-Arg anticodons and which are primary.
- `skiplike/`: Amendment 5. `site_class.tsv.gz` classes every Set A site (skip, inframe, unclear) with the codons each
  frame runs; `summary.tsv` and `decide.txt` give read pairs, lacks and sites per class for both pipelines;
  `S*_mine_*.obs.tsv` are per-site keeps/lacks (minimap2 pipeline) and `S*_recheck_skip.obs.tsv.gz` the aligner-free
  recount at skip sites. S8 = the 8 stranded libraries, S9 = with SRR37726778.
- `followS/`: the stranded follow-up: per-run and pooled per-site counts (`.obs.tsv`), totals and the decision.
- `main/`, `dna/`: the first 10 RNA runs (9 not clearly RNA; see the log) and the paired DNA runs: totals per run and
  pooled per-site counts.
- `strand_probe.tsv`: the strandedness probe per run; `screen_*.tsv` and `*_runs.txt`: run screening and selection.
- `recheck/`, `recheck_S/`: the independent recount without an aligner (reports, summaries, per-site counts).

## CUU context (`cuu_context/`)
- `README.md`: the rule (15 Nov 02:55, before any upstream codon was read) and log, including the corrected sizes.
- `result.txt`, `codon_ratios.tsv`: observed/expected for every upstream codon at skip-like and in-frame-like sites
  (`src/skip_context_chicken.py`). The 7.0-fold there uses the GGN background, later found too generous.
- `recheck/`: the independent re-derivation (`RECHECK.md`; the CGA/CGC/CGT background gives 3.9, and 2.8 with related
  sites counted once).
- `chicken_genomes.tsv`, `chicken_sites_setA.tsv.gz`, `chicken_site_class.tsv.gz`: the 134 chicken-gut TANB77 genomes
  of the site build (100 primary, which give the Set A sites; 94 used after excluding the original count's genomes)
  and their CGG sites with sequence classes.

## Phages (`phages/`)
- `RULES.md`: the phage analysis rules, frozen 15 Nov 17:50 before any target statistic.
- `main/`: phage groups (`T.tsv` the TANB77-host phages with GC and host; `run_ids.tsv.gz` all groups; `recoded_clos.tsv`
  the host genera), per-phage codon statistics (`target.tsv`, `ctrl11.tsv.gz`), and the scores (`score_T1_R.txt`:
  skip ratio and CGG use against other Clostridia phages). T1 = phages with a skipping TANB77 genus as CRISPR host;
  R = phages of other Clostridia orders.
- `posthoc/`: descriptive follow-ups (distinct conserved sites, `sites_t1.txt`; viral families, `vfam.txt`; high-GC
  phages, `highgc*.txt`; host evidence, `evidence.txt`).
- `controls/`: null and planted-recall controls.
- `audit/`: an independent audit with its own code (`REPORT.md`, logs; the 31-of-40 read-through check is task 3).
- `out_of_sample/`: the referee check (rules, report, and outputs; `out_h4.txt` is the k-mer-host test, 20/50 against
  1/146).
- `pfam_junction/`: the Pfam check of proteins read through the skip (`score.txt`: 239/600 against 18/600 and 3/600).
- Phage genomes and host predictions come from UHGV (not included).
