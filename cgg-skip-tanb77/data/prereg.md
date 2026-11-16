# Logged plan: rules, predictions and changes, with the time each was written

This is copied from the project's working log (times are US Eastern, November 2026). Rules were written before the data
they govern were scored; later changes say when and why. Lab-internal agent names are replaced by roles, and file paths
are those of the working project: `scratch/<name>.tsv` and `results/...` tables that support the result are in this
folder's `data/` (see [`README.md`](README.md) for where each one went). Code paths `src/...` and `scratch/<name>.py` are
in `../src/` and `../src/analysis/`.

Other rule files, each with its own timestamps:
- RNA test (12 Nov 07:00, Amendments 1-5): [`rna/RULES_AND_LOG.md`](rna/RULES_AND_LOG.md)
- CUU-context replication (15 Nov 02:55): [`cuu_context/README.md`](cuu_context/README.md)
- Phage analysis (15 Nov 17:50) and its referee check (15 Nov 18:45): [`phages/RULES.md`](phages/RULES.md),
  [`phages/out_of_sample/RULES.md`](phages/out_of_sample/RULES.md); Pfam junction check (15 Nov 18:43):
  [`phages/pfam_junction/RULES.md`](phages/pfam_junction/RULES.md)
- Blinded re-check (11 Nov, verdict written 16:12 before unblinding): [`blinded_recheck/VERDICT.md`](blinded_recheck/VERDICT.md)

The family-wide frameshift scan (one genome per GTDB family) used the gate below; its other flags are not part of this
result and their tables are not included.

## 1. The frameshift scan gate and the CGG test (11 Nov)

## Gate (written 2026-11-11 12:25 EST, amended 12:42 after the independent audit, both before any family or RF-set
result was read; calibrated on controls and injections only)
**Production path** (`src/fs_scan.py`, `src/summarize.py`, `src/run_set.py`): DIAMOND 2.2.8 `blastx -F 15
--range-culling -k 10 -e 1e-10 -c 1`, culled to one hit per range in `fs_scan` (was `-k 1` until 13:15, see log),
against `data/sprot_p50.dmnd` (Swiss-Prot clustered at 50%, the 49,607 clusters
with a bacterial, archaeal or viral member). Each genome is run with its GTDB translation table (4, 25, else 11).
**Break classes:**
- **pair:** opposite-type breaks within 150 nt in one alignment (a detour around an in-frame stop, double indels;
  30 nt until the audit showed detour spacing is flat out to ~150 nt). Excluded everywhere.
- **clean:** not a pair; alignment ≥ 50 aa; both 15-column flanks fully aligned with identity ≥ 0.20; longest run
  touching the tie window ≤ 4; no run ≥ 5 reaching within 6 nt of the window (measured whole, not cut at the 6-nt
  edge).
- **prog:** the clean tests except the two homopolymer ones (prfB's site sits next to its GGGGG Shine–Dalgarno run),
  plus a stop in the abandoned frame within 3 codons of the tie window.
**A genome is flagged** when either route passes:
- **Route S (stop-anchored motif):** the commonest (type, P-site codon, stop) motif among prog breaks occurs in ≥ 6
  distinct subject proteins and in ≥ 0.25 of all prog proteins.
- **Route K (motif-agnostic):** among clean breaks (one per subject protein), the best non-homopolymer 4-mer in the
  30-nt break windows has binomial −log10 p ≥ 6 against random windows from the same genome's aligned loci, is in
  ≥ 8 windows, and its window fraction exceeds the background fraction by ≥ 0.20.
**Calibration** (`results/controls/summary.tsv`, rerun with the 13:15 search; every flag and route call is the
same as under `-k 1`, whose tables are in `results/old_k1/controls/`):
- **No flags among 23 nulls:**
  - E. coli K-12, B. subtilis 168, M. tuberculosis H37Rv;
  - pseudogene-rich M. leprae and *Sodalis*;
  - a CPR (GCA_047456305.1), an Asgard (*Heimdallarchaeum*, GCA_021513715.1) and a Dehalococcoidia
    (*Dehalogenimonas*, GCA_050859015.1);
  - homopolymer indels at 0.1/kb in 4 genomes;
  - random +1 and −1 indels in 1–10% of genes.
  - Null maxima, route S: 6 genes at share 0.05 (*Sodalis*, and random +1 plants at 3%), or 2 genes at a higher
    share.
  - Null maxima, route K: lp 4.8 without random plants. Random +1 plants reach lp 6.5 (TAAA in 34 windows,
    effect 0.15) at 3% and lp 16.3 (CTAG, created by the plants themselves, effect 0.05) at 10%. Only the 0.20
    effect floor keeps these out.
- **Positive:** *Euplotes focardii* is flagged by both routes (S: +1 AAA-TAA in 32 proteins, share 0.33; K: AATA,
  lp 16.6). This is despite the prokaryote-only database; against all of Swiss-Prot (`-k 1` search) it gives 110
  proteins, share 0.47, lp 31.
- **E. coli prfB:** prog, +1, a_stop 2, motif CTT-TGA.
- **Smallest recurring code flagged** (planted motifs: +1 AAA-T or CTG-T, quadruplet GAT+C or CCG+A, −1 after AAG):
  - E. coli (4,000 aligned loci): +1 and quadruplet at 1% of genes; −1 at 3%.
  - Asgard and Dehalococcoidia (780–900 loci): +1 and quadruplet at 3%; −1 not at 3%.
  - CPR (277 loci): +1 and quadruplet at 10%; −1 not at 10%.
- **Per-event recall** (planted event found within ±9 nt, same type and strand):
  - In E. coli, 95–100% of events fall in an aligned locus (71–79% with `-k 1`), and 84–97% of those are found as
    clean.
  - In the divergent genomes, fewer fall in a locus: Asgard 25–39%, Dehalococcoidia 43–58%, CPR 29–71% (few
    events).
  - Tables: `scratch/k10/recall.txt`.
- **Known limits:**
  - −1 codes are weakly detectable.
  - Route K drops only the four pure 4-mers; clean breaks can still sit next to runs of 3–4, so 4-mers such as AAAC
    are enriched by homopolymer indels. At 1 indel/kb (the reviewer's injections) E. coli reaches AAAC lp 24.0 in
    279 windows, but effect 0.12 (< 0.20); the effect threshold is what stops this class, and stage 3 (raw reads)
    is the backstop.
  - Genomes far from Swiss-Prot have low locus coverage, so per-genome loci are reported.
**A flagged genome becomes a candidate** only if, in order:
1. **Shifted windows:** the flagged motif or 4-mer is at most half as frequent in windows 60 nt up- and downstream
   of the same breaks.
2. **Not mobile elements:** fewer than half of the motif proteins are transposases, IS, phage or RF2/dnaX.
3. **Raw reads carry the breaks:** the 25-mer ratio is ≥ 0.5 at ≥ 80% of motif sites with flank depth ≥ 5
   (another agent's `reads_check.py`).
4. **Independent assembly:** the same motif breaks appear at orthologous positions in another assembly of the genus
   or family, and not in a sister clade.
5. **Different aligner:** a subagent re-checks the motif genes independently (miniprot or exonerate).
6. **Prior art:** the motif and lineage are checked against the literature.

**Pooled test** (`src/pool.py`; written 2026-11-11 14:50 EST, amended 15:20 after an independent audit, both before
any pool of family reps was read except Mycoplasmatales, 5 members, not flagged). Small and
divergent genomes have few aligned loci, so a code used in a few percent of genes gives each genome only 2–5 motif
proteins, below the route S floor. A clade that shares the code pools them.
- **Lab prior work:** pooling comes from earlier lab work (another agent's idea, and genus pooling of Arm A in the stop-codon screen)
  (recall for rare stops reached ≥ 0.96 only at ≥ 5 members). Two other agents found that pooling one gene model
  across many genomes counts a conserved orthologous feature once per genome, the pooled conserved-tail confound
 . Here the unit is a distinct Swiss-Prot subject, so a break conserved at
  one orthologous site across a clade counts once.
- **Pools:** family reps (`data/family_set.tsv`) grouped by GTDB order. A pool needs ≥ 2 members used. Members with
  ≥ 100 clean breaks per 1,000 loci, or with no background table, are left out before grouping (the ONT-type indel
  load, log 12:56). Singleton orders are covered only by the per-genome test.
- **Route S:** each subject protein with a prog break anywhere in the pool counts once, under the motif read by the
  most members (ties go to the first member). Thresholds are the same as for one genome (≥ 6 proteins, share
  ≥ 0.25). Before the audit, a subject broken with different motifs in several members counted once per motif, which
  inflated the top motif about 2x in large random pools.
- **Route K:** one clean window per subject protein across the pool. A window from genome g holds 4-mer k with g's
  own background fraction q_g(k). The expected count is Σ n_g q_g(k) and the tail is exact Poisson-binomial, so a
  mix of genomes cannot fake enrichment through composition. Every 4-mer gets the exact tail, and the pool passes
  if any 4-mer meets all three thresholds (lp ≥ 6, ≥ 8 windows, effect ≥ 0.20). The per-genome gate tests only its
  best-lp 4-mer and stays frozen. Large pools give lp ≥ 6 to small effects (ATGA/GAGG gene-junction 4-mers at
  effect 0.03–0.05), so the effect threshold carries the test there. `kmer_junction` reports the share of the 4-mer's
  windows that are −1 junctions of two adjacent genes aligned to one fused subject.
- **Checks:** over the 52 control tables, one-genome pools reproduce `summarize.py`, except that a subject with two
  motifs in one genome now counts once. That lowers *Euplotes* from 32 to 29 proteins (still flagged) and *M. leprae*
  from 5 to 4.
  `scratch/pool_controls.tsv` lists the control pools and `scratch/pool_controls.out.tsv` holds their results.
  - Null pools (all eight null genomes; the three divergent ones; 0.1/kb homopolymer indels in four genomes): no flags.
    The maxima are route S 8 proteins at share 0.07 and route K lp 3.5.
  - At 1/kb homopolymer indels, every member exceeds the break cap. With the cap off, the pool reaches route S 42
    proteins at share 0.06 and GAAA lp 21.6 at effect 0.10, so it is not flagged.
  - Random +1 and −1 plants at 3%: not flagged (lp 4.0 and 2.9).
- **Planted codes in the three divergent genomes pooled** (each of which alone needs 3–10%):
  - +1 AAA-T: flagged at 1% (S 9 proteins at share 0.31; K AAAT lp 8.5, effect 0.40).
  - GAT+C quadruplet: flagged at 1% (K GATC lp 8.0, effect 0.33).
  - −1 after AAG: flagged at 3% (K AAGT lp 7.2, effect 0.23), but not at 1%.
  - These are unrelated genomes. In a real order, related members share subjects, and the first member's window
    wins, so each genome contributes fewer windows (28.5 against 33.5 per genome in the audit's pools). A code in
    10 members of a 148-member order is diluted to about 0.07 of its effect.
- **Empirical null:** the same used family reps are shuffled into random pools with the order pools' sizes (seeds
  1–3, in the same run on the same snapshot of scanned genomes). A clade-specific code is diluted there, so flags in
  shuffled pools estimate the pooled false-flag rate. Order pools share more subjects than shuffled ones (10.6%
  against 4.4% in the audit's RF-set pools), and the one-motif-per-subject rule keeps that from inflating route S.
- **A flagged pool becomes a candidate** only if, in order:
  1. Stage 1 (shifted windows) passes, counted over the pool.
  2. Stage 2 (mobile elements) passes.
  3. Each member holding ≥ 2 motif proteins or excess windows goes through stages 3–6 on its own.
  - Per-member counts are reported. A pool driven by one member is that member's case.
  - For a route K flag, the `kmer_junction` share is reported. If most of the 4-mer's windows are gene junctions,
    the flag is that artifact.

**Side readout, release factors** (`src/prfb_sensor.py`, from the same scan tables). For each genome it reports
whether RF1 and RF2 have a locus, and whether RF2 carries the +1 CUU-UGA autoregulatory shift. A genome that reads UGA
as Trp should have lost RF2 or its shift, and a genome that reads UAG as an amino acid should have lost RF1. This
gives a second view of stop reassignment that does not use the code-inference tools.
  - DIAMOND often aligns straight past the shift (the ~25 aa before it do not pay the frameshift penalty), so the
    shift is tested on the genome sequence: an in-frame CTT-TGA near the RF2 hit's start whose next codon is in the
    frame of the hit's 3' end.
  - Controls: the shift is found in E. coli and *Sodalis* (codon 26), *B. subtilis* 168, CPR, Dehalococcoidia and the
    *Euplotes* assembly's bacterial RF2, and not in *M. tuberculosis* or *M. leprae*. Heimdallarchaeota have neither
    RF1 nor RF2.
  - The DIAMOND-break version missed *B. subtilis* and 133 of 158 family reps whose hit starts near RF2's N-terminus.

**CGG skip test, order TANB77** (written 2026-11-11 15:35 EST, before any TANB77 genome other than the discovery
genome GCA_963931445.1 was scanned; the discovery genome is not counted). Hypothesis: in-frame CGG is not read as a
triplet; the C is skipped and the GGN after it is read as Gly (a +1 shift onto GGN).
- **Test set T:** one genome per GTDB genus in order TANB77 (64 genera; best completeness − 5 × contamination, then
  fewest contigs; for *CALXVA01*, the best genome other than the discovery genome). Scanned like the family reps.
- **Control set C:** the family reps already scanned in class Clostridia outside TANB77 with GC ≤ 0.35. No new
  compute.
- **Per genome:**
  - `cgg`: clean +1 breaks whose old-frame break codon is CGG.
  - `p1`: all clean +1 breaks.
  - `ccg`: an ARAGORN tRNA with anticodon CCG (the CGG reader), from the atlas trna product.
  - Call "CGG skip" if `cgg` ≥ 3 and `cgg / p1` ≥ 0.25. CGG is under 0.2% of codons in these genomes, so random
    indels almost never land on it.
- **Predictions:**
  1. In T, CAG-508 genomes with ≥ 300 loci: a majority call. In C: ≤ 5% call.
  2. In T, a call goes with the absence of `ccg` (one-sided Fisher).
  3. Pooled over T's CGG breaks, the homolog residue opposite the skipped codon (`scratch/shift_residue.py`) is Gly in
     ≥ 50%. Gly is about 7% of residues.
  4. The discovery genome's CGG sites are in its raw HiFi reads (stage 3, `reads_check.py`). The same-species
     assembly GCA_944393025.1 (a chicken-caecum MAG) carries the same CGG + GGN sites wherever the locus is
     assembled (stage 4).
- **Prior lab signals, not yet explained:**
  - The lab's earlier stop-codon screen (Arm A, `../tta-myxobacteria`) flagged CGG at gene ends in 826 TANB77 genomes (919 CGG flags in all of GTDB). Arm B did not
    confirm a stop. A +1 skip ends Prodigal's old-frame gene a few codons after the CGG.
  - codon-hunt's census called CGG = Gly in family CAG-508 (892 genomes, margin 935;
    `work/codon-hunt/results/candidates_family.tsv`). An in-frame alignment puts CGG in the Gly column of the GGN
    that follows.
  - Prior art read so far: CGG reassigned to Gln/Trp in other low-GC Clostridia (Shulgina & Eddy 2021), to Ala in
    *RGIG3102* (KACI) and *RUG521* (codon-hunt), and unassigned CGG in *Mycoplasma capricolum* (Oba et al. 1991).
    In *E. coli*, CGG-rich genes without the cognate tRNA show +1 shifts and CAG misreading. No CGG skip in TANB77
    found yet.

## 2. Log, 11 Nov 15:39 to 16:17: the discovery, the frozen test, the census, turnover, CAG-698, the blinded re-check

- 15:39 **A CGG skip in GCA_963931445.1** (*CALXVA01* sp944393025, family CAG-508, order TANB77, Clostridia; a
  one-contig HiFi MAG, GC 25.9%). It came out of the family scan as a route-K CGGA flag (lp 28.5).
  - 21 of its 32 clean breaks are +1 at an old-frame CGG. The 4 nt from the break are always C + GGN, and CGGA, CGGT
    and CGGG all occur, so the flag's 4-mer is only the commonest of them.
  - Homolog residue at the skipped position (`scratch/shift_residue.py`, output `scratch/shift_residue_445.tsv`):
    Gly at 16 of 20 aligned CGG breaks (Ser 2, Ala 1, Phe 1). At the other 10 breaks it is never Gly, except two
    ACA breaks placed one codon before a CGG.
  - ARAGORN (atlas trna product) finds no tRNA with anticodon CCG. The Arg tRNAs are ACG, CCT and TCT. There are two
    Gly-TCC.
  - In aligned genes, in-frame CGG is 16 of 216,088 codons (CGT 303, CGA 132, CGC 45).
  - Reading: CGG has no reader, and the ribosome skips the C and reads the GGN as Gly. This would explain the stop-codon screen's
    Arm A CGG flags (90% of them in TANB77, 86% in CAG-508) and codon-hunt's CGG = Gly call for CAG-508. Neither
    had been followed up.
  - Frozen test for the rest of the order: Gate, "CGG skip test, order TANB77" (15:35).
  - One change to prediction 4: the raw-read check uses the sample's Illumina run (ERR11860625, NovaSeq, 13.4 Gb)
    in place of the HiFi reads. The HiFi reads built the assembly, so the Illumina reads are the independent check.
- 15:51 **CGG skip test, order TANB77: all four frozen predictions pass** (`scratch/cgg_test.py`; per-genome table
  `scratch/cgg_test_T.tsv`). T = 64 genomes, one per genus; the same-species MAG GCA_944393025.1 is used only for
  stage 4.
  - **1. PASS.** CAG-508 genomes with ≥ 300 loci: 46 of 50 call. The 4 that don't are CAG-492, CAMPYH01, JABCPO02
    and UMGS1901. Other TANB77 families: SIG350 2/2, UBA1234 4/4, CAG-465 2/6, RGIG9151 0/2.
    - Frozen C (6 scanned genomes at 15:35): 0 call (`scratch/cgg_test_C.tsv`). Six genomes have no power, so I
      also scanned all 50 low-GC Clostridia family reps outside TANB77 (`data/clos_lowgc_set.tsv`). That is an
      extension, not the frozen set. Result: 0 of 50 call, and only 2 CGG breaks among 575 clean +1 breaks
      (`scratch/cgg_test_C50.tsv`).
  - **2. PASS.** Calls against `ccg`:
    - call and no CCG tRNA: 54
    - call and CCG tRNA: 0
    - no call and no CCG tRNA: 8
    - no call and CCG tRNA: 2 (both RGIG9151)
    - One-sided Fisher p = C(10,2)/C(64,2) = 0.022.
    - No TANB77 genome has a TCG tRNA, and only RGIG9151 has CCG. 44 of the 50 C genomes have CCG or TCG.
  - **3. PASS.** Homolog residue opposite the skipped codon, pooled over T: Gly at 1,237 of 1,458 CGG breaks (84.8%;
    then Ala 70, Ser 59). Unfrozen controls on the same genomes: other clean +1 breaks give Gly at 40 of 663 (6.0%),
    and breaks placed one codon before a CGG give 20 of 236 (8.5%).
  - **4. PASS, stage 3 (changed to Illumina, logged at 15:39).** ERR11860625: 374,970 of 89.1M reads map to the
    MAG (minimap2 sr). With-C vs without-C 25-mers (`scratch/allele_kmers.py`) at the 21 CGG break sites:
    - 140 reads carry the with-C allele, 3–14 per site, about the flank depth.
    - 0 reads carry the without-C allele.
    - Site 18 has almost no reads at all.
    - `reads_check.py`: at the 14 CGG sites with flank ≥ 5, the median ratio is 1.00; 86% supported, 0% unsupported.
  - **4. PASS, stage 4.** The 21 sites in GCA_944393025.1 (another sample):
    - 17 match the with-C 25-mer exactly.
    - 3 align with the C present and no gap. At sites 15 and 17 the SNP is at the N of CGGN (A↔T): synonymous Gly if
      the C is skipped, a change to the next codon if CGG is read as a triplet.
    - 1 locus is not assembled.
    - 0 match the without-C 25-mer.
  - **Context** (`scratch/skip_context.py`, T without the extra genome, 1,477 sites; not part of the frozen test):
    - In frame-continuous hits, in-frame CGG is 364 of 10.27M codons (0.035 per 1,000).
    - The N of CGGN follows ordinary Gly codon usage: A 67 / C 5 / G 6 / T 21%, against GGN usage of
      GGA 57 / GGC 7 / GGG 6 / GGT 29%.
    - The codon before the CGG is skewed against codons that precede in-frame GGN (χ² 2,343). Enriched: CTT 153 vs
      25.7 expected, TTG 55 vs 11.1, TCA 92 vs 26.6, CCA 84 vs 28.9, GAC 34 vs 8.9. Depleted: TTT 0 vs 55, AAA 21
      vs 99, GAA 23 vs 82, TAT 13 vs 56.
    - CTT before the shift site matches the prfB CUU-UGA slip, where peptidyl-tRNA-Leu re-pairs at UUU. One caveat:
      some of the depletion of homopolymer codons may come from the clean-break filter.
  - **Census** (`scratch/cgg_test_fam.tsv`): of 1,289 family reps scanned so far, only the discovery genome calls.
    High-GC genomes reach cgg ≥ 3 but share ≤ 0.17. In codon-hunt's candidates, CGG = Gly occurs only for CAG-508.
  - **Still to come:** a blinded re-check with a different method (miniprot on 9 coded genomes, including CCG-tRNA
    and non-call controls) and the prior-art search.
- 15:58 **Census: where CGG has no reader** (`scratch/cgg_reader_census.py`; rules fixed in its docstring before
  the run; outputs in `results/cgg_reader/`). Source: codon-hunt's ARAGORN screen of 191,958 R232 species reps
  (archive, read-only: `project-13/results/trna_screen/genomes.tsv`). Usable genomes: completeness ≥ 90,
  contamination ≤ 5, ≥ 30 tRNAs, and at least one other Arg tRNA (ACG, CCT or TCT). 108,450 are usable, and 4,536
  have no CCG or TCG tRNA.
  - TANB77: 347 of 352 have no reader. By family: CAG-508 287 of 289, CAG-465 40/40, UBA1234 17/17, SIG350 3/3.
    The rest (RGIG9151) keep CCG.
  - The only other orders with a majority: Mycoplasmatales (173/235) and RF39 (265/495). Acholeplasmatales is
    38/108. Most orders are at 1–7% (likely missed tRNAs in MAGs).
  - **Losing the reader is not enough for a skip.** 7 of the scanned Mollicutes family reps have no CGG reader in the
    screen. Their clean +1 breaks: 0 of 28 at CGG. In T's genomes with no CCG tRNA: 1,500 of 2,477 (61%). At
    TANB77's rate, 0 of 28 has probability about 1e-11. Mollicutes drop CGG; TANB77 keep CGGN where Gly is coded.
  - In T, the 8 no-call genomes with no CCG tRNA have few breaks (1–10 each). Pooled, 9 of 36 are at CGG (25%),
    against 0 of 8 in the two RGIG9151 genomes.
- 15:58 **Turnover test, predictions fixed before the run** (`scratch/site_turnover.py`).
  - **Question:** if CGGN reads as Gly, inserting or deleting the C costs nothing. Are skip sites then gained and lost
    between relatives like synonymous changes?
  - **Method:** a 121-nt window (gene strand) is taken around:
    - each CGG skip site in each T genome (plus the discovery genome);
    - an equal number of in-frame GGN starts not preceded by C;
    - an equal number of other in-frame codon starts.
    - The GGN and other starts come from frame-continuous hits.
    - Each window is aligned by blastn (dc-megablast) to every other genome. The best HSP covering query 31–91 at
      ≥ 70% identity is read at the centre.
  - **Predictions:**
    1. Same species (the discovery genome against GCA_944393025.1): the C is kept at ≥ 90% of sites. This is the
       method's positive control; stage 4 found 20/20 at the sequence level.
    2. Between genera within a family: the C is lost (a gap in B) at more than 5× the rate of single-base gaps at
       other codon starts.
    3. Between genera within a family: B gains a base right before A's GGN codons at more than 5× the rate at
       other codon starts, and ≥ 50% of those gained bases are C.
    - If 2 and 3 fail while 1 holds, skip sites are conserved, which points to programmed sites rather than a free
      4-nt Gly codon.
- 16:01 **Turnover test: skip sites are gained and lost between relatives.** Outputs: `scratch/site_turnover.tsv`,
  one row per (window, genome B). There were 4,557 windows (1,519 for each kind), aligned with blastn dc-megablast
  (18 s).
  - **Frame check:** the in-frame "other" windows have 0 stops and 0 CGG at the centre among 1,519.
  - **1. PASS.** Same species (discovery genome vs GCA_944393025.1): the C is kept at 41 of 43 sites (95%). The two
    "losses" are hits at 87% and 95% identity, probably paralogs in the 79%-complete MAG.
  - **2. PASS.** Between genera within a family, B has a gap opposite the C at 12,344 of 22,416 CGG sites (55%).
    Gaps opposite other codon starts: 54 of 22,372 (0.24%), about 230×.
  - **3. FAIL as frozen.**
    - B gains a base right before A's GGN at 74 of 27,260 windows (0.27%), against 26 of 22,372 (0.12%) before other
      codon starts: 2.3×, not 5×.
    - Single-base gains are C 42, A 9, T 6, G 0. The C share of all gain events is 42 of 74 (57%), so that half of
      the prediction passes.
    - **Post hoc, not a test:** the 5× threshold was out of reach by construction, and I only worked this out
      afterwards.
      - A gain can only show up where B has its own skip site. T has 1,519 detected sites against 623,372 in-frame
        GGN codons (not preceded by C) in frame-continuous hits, a ceiling of about 0.24% per GGN.
      - Restricted to single-C gains, it is 42 of 27,260 before GGN against 2 of 22,372 before other codon starts.
  - **Net length is a better reading than the centre column.** B's length minus A's over ±15 nt is −1 at 19,528 of
    22,416 between-genus CGG sites (87%), wherever the aligner put the gap. Most "kept, non-C" calls are misplaced
    gaps (6,811 of 7,200). Net −1 in control windows: 0.56% (other codon starts), 0.58% (GGN).
  - **Share of A's skip sites where B keeps the C** (net 0, C opposite), by relation and window identity:
    - same species: 33/33 at 99–100%, 8/9 at 95–99%;
    - same genus, other species: 36–44%;
    - same family, other genus: 8–11% (at 70–90% identity), 29% at 90–95% (n 126);
    - other family: 2–7%.
    - Single-base indels at control codon starts stay at 0–1.5% across all bins.
  - **Reading:** in TANB77 a C before GGN behaves like a near-neutral change rather than a frameshift. Skip sites
    are shared within a species, mostly not between genera, and new ones arise as C insertions before GGN. Possible
    sources of the turnover include near-neutral drift, or a mildly costly site at mutation–selection balance. This
    data can't separate those. Programmed sites (like prfB) would be conserved, and these are not.
  - **Positive control for conservation:** prfB. 14 windows are centred on the TGA of CTT-TGA +1 breaks on RF2
    hits (`scratch/turnover_prfb_windows.py`, output `scratch/site_turnover_prfb.tsv`). The site is kept (net 0, no
    gap) in 15 of 15 comparisons:
    - species 2/2, genus 2/2, family 10/10, order 1/1.
    - Skip sites between genera within a family are kept in about 9%. The method does see a conserved programmed site.
  - The GGN ceiling count is in `scratch/turnover_ggn_pool.py`. Windows, the blastn database and raw hits are in the
    session scratchpad; they can be rebuilt with `site_turnover.py windows` (seed 1).
- 16:03 **CAG-698 test (Bacilli_A, order ML615J-28), predictions fixed before any CAG-698 genome was scanned.**
  - **Why this family:** it is a second lineage where the CGG reader is lost in only part of the family, so it is a
    within-family natural experiment. Another agent's genus pooling (genetic-code-frontier README) found unread CGG in its
    genus RGIG8164.
  - **Set:** `data/cag698_set.tsv`, all 22 usable CAG-698 species reps (census rules of 15:58). 7 have no CCG or
    TCG tRNA (RGIG8164, CAKOMG01 ×2, RGIG3506 ×2, CAG-698 ×2) and 15 keep one. Genus CAG-698 has both kinds.
  - Scanned on the production path. `scratch/cgg_test.py` gives the per-genome call (cgg ≥ 3 and share ≥ 0.25).
  - **Predictions, if reader loss leads to skipping as in TANB77:**
    1. Pooled over the 7 no-reader genomes, the CGG share of clean +1 breaks is ≥ 0.25. Pooled over the 15
       reader-keepers, it is ≤ 0.05.
    2. Calls fall only in no-reader genomes (one-sided Fisher).
    3. If there are ≥ 10 CGG breaks in all, the homolog residue opposite is Gly at ≥ 50% of them.
  - **The other outcome, reported as such:** a no-reader share near 0, as in Mollicutes (0 of 28), means the reader
    was lost and CGG dropped without skip sites.
- 16:08 **CAG-698 result: P1 FAILS. The reader is lost there without skip sites, as in Mollicutes.**
  - No-reader genomes: 1 of 33 clean +1 breaks at CGG (3%). Reader-keepers: 3 of 106 (2.8%). 0 calls, so P2 and P3
    are not testable (`scratch/cgg_test_cag698r.tsv`).
- 16:08 **In-frame CGG against skipped CGG, by reader status** (`scratch/inframe_cgg.py`: in-frame codons of
  frame-continuous hits, ≥ 20 codons from each end; reader status from codon-hunt's ARAGORN screen).
  - CGG per 1,000 in-frame codons, median per genome:

    | group | no reader | reader |
    |---|---|---|
    | TANB77 | 0.013 (n 45) | 0.384 (n 2) |
    | CAG-698 | 0.035 (n 7) | 1.67 (n 15) |
    | Mollicutes family reps | 0.17 (n 7) | 0.96 (n 3) |
    | low-GC Clostridia (C50) | 0.17 (n 6) | 0.39 (n 44) |

  - Every lineage that lost the reader also dropped in-frame CGG. So avoiding CGG is not what sets TANB77 apart.
  - What sets it apart is where its remaining CGGs sit. Pooled over genomes, skipped CGG (clean +1 breaks) against
    in-frame CGG:
    - TANB77 no reader: 1,500 against 262 (5.7 : 1). This counts 45 genomes with a screen call plus 18 with none in
      the screen, all in no-reader families.
    - TANB77 with a reader: the 2 RGIG9151 genomes, which keep tRNA-CCG, have 0 skips against 139 in-frame.
    - CAG-698: no reader 1 against 271; reader 3 against 12,788.
    - C50: no reader 0 against 1,303; reader 2 against 8,450.
    - Mollicutes no reader: 0 of 28 breaks against 276 in-frame.
  - So only in TANB77 are most CGGs in genes skips rather than codons.
  - Caveat: the two counts come from different hit sets (broken against frame-continuous hits), so the ratio is
    descriptive, not a rate.
- 16:08 **Stage 6 (prior art), from a subagent search**
  (scratchpad `lit/`: KACI table `lit/kaci/kaci_tanb77.tsv`, Codetta rows `lit/codetta_supp/codetta_tanb77.tsv`).
  - **The signal is already in public data.**
    - KACI's Zenodo calls (Melnykov, doi:10.5281/zenodo.19323489, v2 of 25 Jun 2026) give CGG = G in 941 of 2,763
      TANB77 assemblies (CAG-508 926). All of them carry KACI's own "is interesting" flag.
    - The discovery genome GCA_963931445.1 and its same-species MAG GCA_944393025.1 are both called CGG = G
      (checked).
    - Every CGG = G call in the file that is in GTDB R232 is TANB77.
    - Codetta (Shulgina & Eddy supplement) leaves CGG at '?' in 96 of 113 TANB77 assemblies, the highest rate of any
      Clostridia order.
  - **The interpretation is not published.**
    - Read in full: Codetta 2021 and 2023, the Shulgina thesis, the KACI preprint and Fakih NAR 2026. None reports a
      +1 skip at CGG, or CGG = Gly in TANB77.
    - The KACI preprint only notes uncertain CGG calls in Clostridia, with no tRNA-UCG, and leaves them "for a
      separate investigation".
    - The final PNAS version of KACI (doi:10.1073/pnas.2610659123) could not be opened (403/429) and is **unchecked**.
    - No NCBI table (gc.prt v4.6) touches CGG.
  - **Precedents:**
    - lineage-wide first-base skips in organelles only: *Perkinsus* mitochondria AGG/CCC (Masuda 2010); dinoflagellate
      mitochondria (Gornik 2022);
    - Ty1 CUU-AGG-C slip at a hungry Arg codon (Belcourt & Farabaugh 1990);
    - E. coli +1 shifts at CGG and AGG under tRNA shortage (McNulty 2003, Gurvich 2005);
    - FMR1 CGG-repeat +1 shift from the Arg frame to the Gly frame (Wright 2022).
  - **Reading:** a frameshift-blind codon-reassignment tool scores the skipped C plus GG as "CGG aligned to Gly".
    So KACI's G calls are very likely this same signal, read as a reassignment. What is new here:
    - it is a skip and not a reassignment (breaks, read-checked C, ~0 in-frame CGG);
    - the mechanism clue (CTT before CGG);
    - turnover of the sites;
    - reader loss without skipping in Mollicutes, CAG-698 and other orders.
  - Stage 6 stands as **partial**: data-level prior art exists, and the interpretation looks new pending the KACI
    PNAS text.
- 16:11 **The stop-codon screen's Arm A CGG flags are this signal** (read-only, `work/new-stops-prok/results/gtdb/armA_flags.tsv`).
  - 826 of Arm A's 919 CGG flags are TANB77, out of 1,017 TANB77 genomes Arm A scored. No other order has more than
    11.
  - Arm A's statistic is the codon's count within 100 codons of a Prodigal gene end, against the interior.
  - A skip produces exactly this: after the skipped C, frame 0 soon reaches a stop. In GCA_963931445.1, 18 of the
    21 clean CGG skip sites lie within 100 codons of an atlas pyrodigal gene end (median 12 codons). That compares
    with 6 of 11 other clean breaks (median 79).
  - Arm A counts n3 = 102 CGG there against 7 in the interior.
  - That screen's Arm B did not pass CGG in any of the 826 genomes.
  - So the same sites show up as "CGG = Gly" (codon-hunt, KACI), as "CGG uncertain" (Codetta, 85% '?') and as
    "CGG = stop-like" (Arm A).
- 16:11 **No quadruplet tRNA in the discovery genome.**
  - ARAGORN (atlas) finds 36 tRNAs.
  - Arg: ACG, CCT, TCT.
  - Gly: GCC, and TCC ×2 (74 nt, no extra anticodon-loop base visible by length).
  - Nothing reads CGG, and no tRNA-Gly could read CGGN as a quadruplet (that would need an anticodon NCCG).
  - Caveat: ARAGORN may not report tRNAs with 8-nt anticodon loops as such. (Withdrawn 2026-11-12: pyaragorn does report them, see the 2026-11-12 tRNA entry.)
  - With no dedicated tRNA, the hungry-codon skip (P-site slip or out-of-frame A-site read) remains the simpler
    model.
- 16:14 **Skip sites avoid ribosomal proteins (exploratory, not pre-registered).**
  - Set: TANB77 T, 65 scans. Genes are assigned by the Swiss-Prot mnemonic of the hit.
  - Ribosomal proteins (RL*, RS*): 3.0% of aligned residues, but 0 of 1,500 CGG skips. Other clean breaks: 16 of
    1,189 (one-sided Fisher p 2e-6).
  - Other core genes (rpoA/B/C, aminoacyl-tRNA synthetases, EF-G, gyrB, dnaG and similar): 8.6% of aligned residues,
    23 of 1,500 skips, against 35 of 1,189 other breaks. Combined with ribosomal proteins: p 1e-5.
  - So skips are not barred from essential genes: 4 are in rpoB/C and 16 in tRNA synthetases. In the discovery
    genome one sits in its single rpoC hit (Swiss-Prot RPOBC_HELHP, residue 2,097), so translation must pass that
    site.
  - But skips are rarer in highly expressed genes than other breaks are. That fits a mild cost, as the turnover
    test suggested (16:01).
  - Caveat: the denominator is aligned residues, not GGN opportunities.
- 16:17 **Stage 5 PASS: blinded re-check with a different method** (subagent; `recheck_cgg/VERDICT.md`, written
  16:12; scripts in `recheck_cgg/scripts/`). Nine genomes were given as G1–G9 with no roles. The method was miniprot
  `-S` against 8 standard-code Firmicutes isolates, with every CGG and stop residue masked to X. It did not use
  DIAMOND or `sprot_p50`.
  - **Unblinded** (the key is the label-to-accession table in the agent prompt):

    | label | genome | lineage | CGG reader | my pipeline (CGG / clean +1) | blind verdict |
    |---|---|---|---|---|---|
    | G1 | GCA_946496695.1 | TANB77 CAG-508 CALWZV01 | none (agent's tRNAs) | 129 / 170 | +1, Gly (G127 of 134 consensus) |
    | G2 | GCF_012317385.1 | Peptostreptococcaceae *Romboutsia_E* | none in screen | 0 / 11 | Arg, low confidence (26 sites; tRNA set incomplete) |
    | G3 | GCA_900553755.1 | TANB77 CAG-508 UMGS1901 | none | 2 / 10 (not called) | +1 (13 of 14), residue undetermined |
    | G4 | GCA_963931445.1 | TANB77 CAG-508 CALXVA01 (discovery) | none | 21 / 31 | +1, Gly (G25 of 26) |
    | G5 | GCF_005473905.2 | *Ruminiclostridium herbifermentans* | CCG | 0 / 6 | Arg |
    | G6 | GCA_964631135.1 | TANB77 RGIG9151 | CCG | 0 / 4 | Arg (no CGG frameshift at 87 sites) |
    | G7 | GCA_009785865.1 | TANB77 CAG-508 WQVL01 | none | 149 / 270 | +1, Gly (G134 of 148) |
    | G8 | GCF_900187165.1 | *Clostridium_G cochlearium* | TCG | 0 / 6 | Arg, moderate confidence |
    | G9 | GCA_949283575.1 | TANB77 UBA1234 | none | 24 / 32 | +1, Gly (G14 of 16) |

  - **All 9 frame verdicts match my pipeline's expectation.** The 5 TANB77 genomes with no reader shift +1 at
    80–98% of their aligned CGG. The TANB77 genome with a CCG tRNA (G6) and the 3 outgroups do not. The residue is
    Gly in 4 of the 5; G3 has too few sites. No genome shows CGG as a stop, an in-frame reassignment, or a −1.
  - **Controls passed in the blind run.**
    - The prfB +1 at CTT-TGA-C was found in all 8 genomes where RF2 aligns.
    - Planted +1 and −1 edits in G6: 149 of 150 each.
    - Planted "C before a Gly codon": 146 of 150 detected, residue G128.
    - Planted "CGG plus one base at an Arg codon": 99 of 100 detected, residue R84.
    - So the readout separates Gly from Arg at single sites.
  - **Extra support from the blind run.** The base after CGG follows the third-base usage of Gly codons (G1 A/T/C/G
    .58/.30/.05/.07 against .61/.26/.05/.08), which fits "C, then a GGN codon".
  - **Discordance: G3 shows my census under-calls.** My pipeline did not call UMGS1901 (2 CGG of 10 clean +1 breaks;
    the rule needs ≥ 3 and share ≥ 0.25). Miniprot finds a +1 at 13 of its 14 aligned CGG, and 14 of its 18
    frameshift events contain CGG. The clean-break filter loses sites, so the "54 of 64 TANB77 genera call" count is a
    lower bound. The blind run makes no false positive among the 4 non-skip genomes.
  - **Not separable from DNA** (the agent's caveat 1): skipping C and reading GGN with tRNA-Gly, against decoding
    CGGN as a quadruplet. Both give Gly. ARAGORN finds no quadruplet anticodon (NCCG) in any skip genome.

## 3. Log, 12 Nov 06:05: tRNA check in the pig genomes

Run while the PXD077676 spectra download; files in `scratch/trnascan_pig/`.
- **Callers.** tRNAscan-SE 2 (`-B`, bacterial model, secondary structures) and ARAGORN (pyaragorn 0.3.0, code 11,
  linear) on the 24 pig TANB77 representatives (`data/proteomics/mgnify_pig/`).
- **Positive control first.** The real MGYG000302328 tRNA-Gly(TCC) with a G inserted after C36 (an 8-nt loop that
  could read GGGN), in random flanks (`control_8loop.py`). ARAGORN calls it with anticodon `tccg`, length 4,
  (Gly, Arg). tRNAscan-SE draws the 8-nt loop. So both callers can show an 8-nt loop. This withdraws the 16:11
  caveat above.
- **No 8-nt anticodon loop.**
  - tRNAscan-SE: 744 calls; no 8-nt loop in any structure.
  - The 21 tRNA-Thr(CGT) drawn with a 9-nt loop are an A31·C39 drawing artifact. Their loop TTCGTAA is standard.
  - ARAGORN: 695 calls with a 7-nt loop and one 6-nt loop (MGYG000302432, `tc`, Glu/SeC). None has an 8-nt loop.
- **No CCG or TCG anticodon of any type.** tRNA-Arg anticodons: ACG in 24 genomes, CCT in 22, TCT in 21.
- **Intron-containing tRNAs.** tRNAscan-SE calls 15 of 149 to 320 bp, in 10 genomes. Excision by brute force
  (best ARAGORN energy, `intron_trnas.py`) gives Ser(AGA) ×9, Thr(AGT) ×4, His(ATG) ×1 and Trp(CCA) ×1. None can
  read CGG.
- **Short Undet calls.** The 12 short calls (54 to 115 bp) have no anticodon arm in tRNAscan-SE's drawing.
  ARAGORN, given each with 150 nt flanks (`undet_short.py`), calls nothing in 11. The 12th is the Glu/SeC call
  above.
- **82 chicken TANB77 genomes** (`data/proteomics/mgnify/`, ARAGORN only, `aragorn_tally.py`).
  - 2,077 calls: 2,074 with a 7-nt loop, two with 8-nt loops (Cys/Met `ACAT` and Stop/Ile `TTAT`, neither NCCG)
    and one with a 6-nt loop (Leu/Pro `AG`).
  - The only tRNA-Arg(CCG) is in MGYG000309316, already marked as a CGG reader in
    `results/proteomics/genomes.tsv`.
- **Reading.** Apart from the known reader genome, nothing in these genomes decodes CGG in frame. No tRNA-Gly could
  read CGGN as a quadruplet. The hungry-codon skip (no cognate tRNA, so the ribosome slips +1 and reads GGN) remains
  the simplest model. Not tested: modified or misacylated tRNAs, which sequence cannot show.


## 4. Log, 15 Nov 03:27-03:28: the CUU context and protein data

- **G37 tRNA-Gly(UCC) lead: killed** as an explanation of the skip (`results/trna_gly37/README.md`).
  - Not needed: WQVL01 has 149 skips and only A37, by two callers.
  - Not enough: it is common in Actinomycetota that keep a CGG reader.
  - It stands as an observation only: one old paralog in TANB77.
- **CUU before skipped CGG: replicated** by a rule fixed in advance on 94 chicken genomes (`results/skip_context_chicken/`):
  7.0-fold against the GGN comparison.
  - An independent re-derivation showed that comparison is too generous: 3.9 against CGA/CGC/CGT, and 2.8 with related
    sites counted once; in-frame-like 1.1.
  - The excess is in CAG-508 and UBA1234, and absent in CAG-465.
  - The post's exploratory line overstated the size the same way. It is rewritten, with two fact-checks
    (`post_draft/rna/factcheck5`, `factcheck6`), all fixes applied.
- The chicken RNA chain is queued behind a load gate.

### 2026-11-15 03:28 — Protein evidence: no public dataset has the power (checked from the 12 Nov survey's tables)
- **PXD077676, all 64 runs:** the submitters' table has 1,596 strict TANB77-only PSMs. The three runs already searched
  hold 1,037 of them (`scratch/proteome_survey/out/pxd077676.runs_bins.tsv`). Searching the other 61 runs would scale
  the re-check's yardstick, about 1.2 skip peptides under a full skip, to about 1.8 (1.2 × 1,596/1,037). The rule
  needs S ≥ 5, so it is not worth the downloads.
- **Every other public gut metaproteome** in the 12 Nov survey fails the abundance gate
  (`scratch/proteome_survey/REPORT.md`). The laying-hen caecal set PXD077115/PXD077163 gives about 40 TANB77-only
  PSMs in its best 3 runs.
- So protein evidence for the skip needs new data, for example deep metaproteomics of TANB77-rich chicken caecal
  content, or a culture. The post's closing line already says so. Nothing more is planned on proteins.
- Still queued: the chicken RNA replication, behind the load gate (load 56-62 at 03:27).
