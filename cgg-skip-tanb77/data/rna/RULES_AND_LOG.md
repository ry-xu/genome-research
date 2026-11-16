# Is the C of CGG in TANB77 messenger RNA? (RNA test of the skip's mechanism)

Main analysis agent. Rule written 2026-11-12 07:00, before any RNA read was mapped to any TANB77 genome.

## Question
TANB77 skips the C of CGG and reads GGN as Gly. Three mechanisms fit the DNA (post, "Not known"): the ribosome
slips +1 at CGG, a tRNA reads the four bases CGGN, or the C is removed from the messenger RNA before translation
(editing; transcriptional slippage would look the same). The first two act in translation, so the mRNA keeps the
C. Under the third, most transcripts lack it at CGG sites, and keep it at other Arg codons. Reads from public pig
gut metatranscriptomes can tell these apart. They cannot tell the two translational mechanisms apart.

## Sites (fixed now; `src/rna_cgg.py sites`)
- Genomes: the 24 TANB77 species representatives of MGnify pig-gut v1.0 (`data/proteomics/mgnify_pig/`), as in the
  PXD077676 build. Primary genomes (no CCG or TCG tRNA-Arg anticodon; `results/proteomics_pig/genomes.tsv`): 22.
- Each catalogue CDS is walked from its start in the skip reading, with the same walk as `src/proteomics_db.py`: an
  in-frame CGG skips the C, then the walk continues in the +1 frame until a stop.
- **Set A (the test):** the C of every in-frame CGG met by a walk in the 22 primary genomes.
- **Set C (control):** the C of every in-frame CGA, CGC or CGT met by the same walks. These are Arg codons read in
  frame by tRNA-Arg(ACG). Same first two bases, same genes, same reads.
- **Set A0 (reported only):** the in-frame CGGs of the 2 non-primary genomes, which have a CGG-reading tRNA.
  *Correction (2026-11-12 12:24, fact-check 2, `post_draft/rna/factcheck2/REPORT.md`): wrong. None of the 24 has a
  CGG-reading tRNA (`cgg_reader` 0 for all). The 2 are non-primary because their tRNA-Arg calls are incomplete,
  so the absence of a reader cannot be shown: in the screen behind `genomes.tsv`, MGYG000301655 has no tRNA-Arg
  called and MGYG000301930 only CCT. (The 06:05 check with tRNAscan-SE and ARAGORN finds ACG in all 24, and CCG or TCG in none.) The
  same applies to "no CCG or TCG tRNA-Arg anticodon" above: primary also requires the CGN reader ACG to be called.
  Set A0 was reported only, and no result used it.*
- **Set B (reported only):** the Set A sites of MGYG000302328 within 3 nt of a clean +1 break
  (`results/pig_check/`, the family-scan rule).
- Positions are deduplicated per genome (walks that share an end), and sites within 30 nt of a contig end are
  dropped.
- **Site run:** the maximal run of identical bases on the genome's plus strand that contains the C's position. For
  a minus-strand gene the C shows on the plus strand as a G. An aligner can place a 1-nt deletion anywhere in such a
  run.

## Reads (fixed now)
- Mapped with minimap2 -x sr -a --sam-hit-only to the 24 representatives' contigs together (one index). Reads are
  streamed from NCBI SRA with `fastq-dump --stdout --split-spot --skip-technical`, as another agent's `rna_stream.py`. Only
  mapped reads are kept.
- A read counts if it is a primary alignment (not secondary or supplementary) with MAPQ ≥ 20 and identity
  1 − NM / (aligned columns) ≥ 0.97.
- **Observation at a site:** the read must align ≥ 20 reference nt on each side of the site run, with no soft clip
  inside those 20 nt. It is then classed as:
  - **keeps:** no insertion or deletion within the run or 5 nt either side;
  - **lacks:** exactly one 1-nt deletion inside the run, and no other indel within the run or 5 nt either side;
  - **other:** anything else (reported, not used).
  - The base at the C is not quality-filtered, so neither class is filtered more than the other.
- One observation per spot (both mates of a pair can cover a site). If the mates disagree the spot is "other".
- **D** for a site set = lacks / (keeps + lacks), pooled over its sites.

## Decision (fixed now; on Set A against Set C, pooled over all streamed runs)
- **Powered** if Set A has ≥ 200 observations over ≥ 30 sites and Set C has ≥ 1,000.
- **Transcripts keep the C**, if the upper end of D_A's 95% interval is < 0.10. The interval comes from a bootstrap
  over sites (10,000 resamples). This would mean at most a tenth of TANB77 transcripts lack the C at CGG, so
  editing and transcriptional slippage are not how the skip is made; it happens in translation.
- **Transcripts lack the C** if D_A ≥ 0.25 and the lower end of D_A's 95% interval is ≥ max(5 × D_C, 0.05).
  - That is RNA-level evidence, but a genomic deletion of the C in the sampled strains would look the same.
  - Before any claim, the same sites are checked in DNA reads from the same samples. That rule is written before the
    DNA is looked at.
- Anything else: **inconclusive**, with D_A, D_C and the bootstrap interval reported.
- Reported, not part of the decision:
  - D for A0 and B;
  - D per genome;
  - the sites with ≥ 10 observations and their D;
  - the "other" class;
  - D_A divided by the injection test's ratio of call rates (below).

## Known-answer test, before any real read is counted (fixed now; `src/rna_cgg.py inject`)
- Simulated 150-nt reads from both strands of the 22 primary genomes, at 2,000 Set A and 2,000 Set C sites drawn at
  random:
  - the site run starts at a random offset between 25 and 115;
  - half the reads lose one base of the site run;
  - every base is substituted with probability 0.005.
- They go through the production path: FASTQ → minimap2 → `rna_cgg.py count`.
- **Pass** if, among reads that make an observation:
  - ≤ 1% of keep-reads are called "lacks" and ≤ 1% of lack-reads are called "keeps", in each set;
  - the share of lack-reads that make an observation is ≥ 0.8 × the share of keep-reads, in each set.
- If it fails, the filters are fixed and the test rerun before any real read is counted. The change is logged.

## Run selection (fixed now; only abundance is used, never a site observation)
- Candidates: the 1,330 ENA runs with library_source METATRANSCRIPTOMIC and pig gut metagenome or pig host
  (`ena_pig_metatx.tsv`, queried 2026-11-12 06:59).
- **Screen:** for each study, its largest run (by read count). The first 1,000,000 spots are streamed and mapped.
  - Score = reads that pass the read filters on the 22 primary genomes, per million spots. `count` is not run on
    screen reads.
  - Then the 5 top-scoring studies get every run screened the same way, up to 20 runs per study (largest first).
- **Main set:** the 10 top-scoring runs over all screened runs, streamed whole, at most 4 at a time.
  - If Set A has < 200 observations after that, the next 10 runs by score are added before D is computed. Only
    observation totals are looked at for that step (`count --totals-only`).
  - The cap is 30 runs.

## Amendment 1 (07:07, before any real read was streamed or mapped)
- **Short-read runs only.** The 239 MinION runs (PRJNA784703, PRJNA813149) are dropped from the candidates, leaving
  1,091 runs in 42 studies. Reasons: `minimap2 -x sr` is a short-read preset, and nanopore errors are mostly
  homopolymer deletions, which is the "lacks" class itself. Set C would only partly absorb that, since the error rate
  depends on run length and context.
- **Spots, not reads.** "Per million spots" uses the spots streamed: min(1,000,000, the run's ENA read_count, which
  counts spots).

## Log
- 07:07 Sites built (`sites.tsv`): A 2,871, C 32,501, A0 238, B 30 (from 32 clean +1 CGG breaks). Every
  row's codon checked against the genome (0 wrong). Site runs: A 2,549 of length 1, 295 of 2, 27 of 3 or more.
- 07:07 Known-answer test (`data/rna_cgg/inject/`, seed 1, production path): **PASS**.
  - Set A: keep-reads 1,949 of 1,993 observed, 0 called lacks; lack-reads 1,919 of 1,994 observed, 1 called keeps
    (0.05%), 21 other. Call-rate ratio 0.984.
  - Set C: keep-reads 1,941 of 1,995, 0 called lacks; lack-reads 1,917 of 1,993, 1 called keeps (0.05%), 14 other.
    Ratio 0.989.
- 07:20 Screen, step 1: 42 study-largest short-read runs (`screen_runs.tsv`), first 1,000,000 spots each, 0
  failed. Scores are in `screen_scores.tsv`.
  - Top 5 studies: PRJNA615303 (96,496 filtered reads per M spots), PRJNA595469 (30,720), PRJNA1203610 (17,171),
    PRJNA1371775 (12,467), PRJNA1441222 (6,923).
  - 99% of hits in the screen SAMs had MAPQ < 20 (rRNA multimappers). From now on `rna_stream.py` keeps only primary
    alignments with MAPQ ≥ 20, the only ones `rna_cgg.py` reads, so nothing used by the rule changes.
- 07:20 Screen, step 2: all runs of the top 5 studies, up to 20 each, largest first: 60 runs
  (`screen2_runs.tsv`).
- 07:41 Screen step 2 done: 60 runs, 0 failed; 97 runs screened in all (`screen_scores_all.tsv`).
- 07:41 **Main set** (`main_runs.txt`), the 10 top-scoring runs, streamed whole into `data/rna_cgg/main/`
  (4 runs at a time, 2 streams each):
  - PRJNA615303: SRR11426316, SRR11426305, SRR11426303, SRR11426304, SRR11426322, SRR11426327, SRR11426338
    (3.0 to 4.5 M spots each);
  - PRJNA595469: SRR10686660 (72,200 spots);
  - PRJNA1203610: SRR31853479, SRR31853486 (82,303 and 93,646 spots).
  - Seven of the ten are from one study. The rule does not cap runs per study, so this stands; D per run is reported.
- 07:41 Deleted `data/rna_cgg/screen/` (3.6 GB of screen SAMs; scores are in `screen_scores_all.tsv`, the .n/.err files in `screen_meta.tgz`; regenerable by re-streaming).
- 07:42 Independent re-check launched (fresh agent, alignment-free exact probe matching on the same 10 runs, own site build; brief `recheck_brief.md`, output `recheck/`). It is told nothing about any count.

## Amendment 2 (08:07; written after seeing only observation totals, before any keeps/lacks split)
- **What was seen.** `count --totals-only` on the 10 main runs, as the rule allows: Set A 1 observation at 1 site,
  Set C 552 at 16 sites. Then abundance by genome feature, which uses no site and no class:
  - SRR11426316: of its filtered reads, 873,049 fall in LSU and 255,218 in SSU rRNA features.
  - The old main set's runs have 0 to 34 filtered reads per M spots that overlap a CDS and no ncRNA feature.
  - No keeps/lacks split for any set or site has been looked at.
- **Why the score changes.** The screen score counted every filtered read, and in these runs that is almost all rRNA
  (unique to one genome, so MAPQ ≥ 20). It does not measure TANB77 mRNA, which is what makes observations.
- **New screen score:** filtered reads on the 22 primary genomes whose alignment overlaps a catalogue CDS and no
  ncRNA feature, per million spots streamed (`rna_cgg.py score`, column `score_cds`).
- **Selection redone with the new score**, the same steps as before:
  - step 1: the 42 study-largest runs, first 1 M spots;
  - step 2: every run (up to 20, largest first) of the 5 top-scoring studies;
  - main set: the 10 top runs, whole, plus 10 more if Set A has < 200 observations (totals only), cap 30.
- The old main set's runs count only if the new score selects them. The old rule's "add the next 10 by score" step
  is not run.
- **Not powered** stays a possible outcome, and is reported as such.

## Amendment 3 (08:19; after the re-check's site report, before any keeps/lacks split)
- **What was seen.** Nothing new about keeps or lacks. The re-check (`recheck/RECHECK.md` §7) found:
  - a predicted CDS inside a 23S rRNA (MGYG000300925_01647), whose Set C sites collect rRNA reads;
  - from the stopped old streams, by location only, that most Set C observations sat at those sites.
- **Sites in ncRNA are dropped.** A site is dropped if its observation window overlaps any `ncRNA` feature of the
  catalogue GFFs (rRNA, tRNA, tmRNA, RNase P...). The window is the site run plus 20 nt each side. Reads there are not
  mRNA.
  - This drops 5 Set C sites, all in MGYG000300925_01647, and no Set A, A0 or B site.
  - The sets are now A 2,871, C 32,496, A0 238, B 30. `sites_v1.tsv` keeps the old file.
  - A first version dropped every site of a CDS overlapping an ncRNA. That also dropped one gene (rluD,
    MGYG000301514_01167) for an 8-nt overlap with RNase P RNA far from its site. It was replaced before any count.
- **Set C, described correctly.** Set C sites are met by the same walks in the same genomes, not mostly the same
  genes. 29,430 of the 32,496 are met only by walks with no CGG, so D_C measures the genomes' background rate.
  - Column `walk_has_cgg` marks the 3,066 met by a walk with a CGG. It agrees with the re-check's build at every
    site.
  - D_C for those 3,066 is reported, not part of the decision.
- **Added to the reported extras:**
  - per run, the share of observations whose read 1 is sense to the gene (strandedness). A share near 0.5 in a
    stranded library points to DNA reads, which keep the C at both sets.
- The Decision is unchanged.

## Amendment 4 (08:20; after the CDS screen scores of step 1, before any site observation from the new runs)
- **Main-set runs are streamed up to their first 20,000,000 spots** (whole if smaller), not whole.
  - Reason: compute. The studies the CDS score picks have runs of 40 to 290 M spots. Ten whole runs would take hours
    to stream here and more than half a day for the re-check's exact matcher (about 6,000 spots/s per process).
  - The cap was chosen from abundance alone: at the step-1 scores (up to about 2,000 CDS reads per M spots), 20 M
    spots give tens of thousands of TANB77 CDS reads per run.
- Everything else stands: the "10 more runs" step if Set A has < 200 observations (totals only), cap 30 runs, the
  same cap of 20 M spots for added runs.
- `main_runs.txt` gives the spots streamed per run (min(20 M, the run's spots)) in its third field.

## Log (continued)
- 08:20 Sites rebuilt under Amendment 3 (`sites.tsv`): A 2,871, C 32,496, A0 238, B 30. Columns 2–10
  match `sites_v1.tsv` row for row except the 5 dropped Set C rows.
- 08:20 Screen redone, step 1 (`screen_scores_cds.tsv`): 42 study-largest runs, first 1 M spots each, 0
  failed.
  - Top 5 studies by `score_cds`: PRJEB86261 (2,075 CDS reads per M spots), PRJNA1441222 (2,000), PRJEB100724 (343),
    PRJNA1432673 (314), PRJNA1129603 (199).
  - PRJNA615303, which gave seven of the old main set's ten runs, scores 28.
- 08:20 Screen step 2 started: the other runs of those 5 studies, up to 20 per study, largest first.
  That is 67 runs (`screen3_runs.tsv`), streamed into `data/rna_cgg/screen/`.
- 08:42 Screen step 2 done: 67 runs, 0 failed. All 109 runs screened are scored in
  `screen_scores_cds_all.tsv`.
- 08:42 **Main set** (`main_runs.txt`; the old list is `main_runs_v1.txt`): the 10 top runs by
  `score_cds`, first 20 M spots each, into `data/rna_cgg/main2/`, 4 runs at a time, 2 streams each.
  - PRJEB100724: ERR15738679 (4,409 CDS reads per M spots).
  - PRJEB86261: ERR15873938, ERR15873922, ERR15873927, ERR15873933, ERR15873934, ERR15873936, ERR15873944,
    ERR15873926 (2,024 to 3,150).
  - PRJNA1441222: SRR37726778 (2,000).
  - Eight of the ten are from one study. D per run is reported.
  - No old main run is re-selected; the best old one scores 28.
- 08:42 Deleted `data/rna_cgg/main/` (4.5 GB, the old main set's filtered SAMs; not used under Amendment 2)
  and `data/rna_cgg/screen/` (screen SAMs of the 109 runs; scores in `screen_scores_cds*.tsv`, .n/.err files in
  `screen_meta_v2.tgz`). Both can be regenerated by re-streaming.

## DNA check (08:49; fixed before any RNA observation of the main set has been counted)
The rule's Decision asks for this before any "lacks" claim. A C missing from the genome would look the same in RNA.
It is fixed now, while the main-set streams and the re-check are still running and no count has been run on them.
- **Reads.** PRJEB86261 sequenced caecum digesta of 24 pigs (J001–J024), from the same day, as both metagenome
  (DNA) and metatranscriptome (RNA). The pairs come from ENA sample attributes (`host subject id`,
  `host body site` Caecum_S2) and are listed in `prjeb86261_digesta_pairs.txt`.
  - The DNA runs are the 8 that pair by pig with the 8 PRJEB86261 main-set runs (`dna_pairs.txt`).
  - First 20 M spots each, the same path as the RNA: `rna_stream.py` into `data/rna_cgg/dna/`, then
    `rna_cgg.py count` with the same `sites.tsv`.
- **Statistic.** D_A and D_C in DNA, with the same bootstrap over sites. It is powered if DNA Set A has ≥ 200
  observations over ≥ 30 sites.
- **If the RNA decision is "lacks":**
  - **RNA-level** (editing or transcriptional slippage) if the upper end of DNA D_A's 95% interval is < 0.05, and the
    "lacks" criterion also holds on the 8 paired RNA runs alone.
  - **Genomic** if DNA D_A ≥ half of RNA D_A on the paired runs (point estimates). The C is then missing from the
    sampled strains' DNA, and the RNA says nothing about the mechanism.
  - **Mixed or unresolved** otherwise, or if DNA is not powered. No mechanism is claimed.
- **If the RNA decision is "keeps" or "inconclusive":** DNA D_A and D_C are reported, not used. They give the
  pipeline's floor where the true D is 0, unless strains lack the C.
- **Reported in every case:** per-site D in DNA and RNA at sites with ≥ 10 observations in both.

## Log (continued)
- 08:53 Spot counts from SRA (`runinfo`): ENA's `read_count` for these ERR runs is reads, not spots. For
  ERR15873922 it is 36,935,238, exactly 2 × its 18,467,619 spots. The re-check found this.
  - Only ERR15873922 is under 20 M spots, so its line in `main_runs.txt` now reads 18,467,619 (the stream with
    -X 20,000,000 stops at the run's end).
  - The screen scores are not affected: they divide by the spots actually streamed (fastq-dump's "Read N spots").
  - All 17 ERR runs checked (the 9 ERR main runs and the 8 DNA runs) have an ENA count of exactly 2 × their spots,
    so a study's "largest first" order is very likely unchanged.
  - The 8 DNA runs have 43.9–65.7 M spots each.
- 09:04 Main set streamed (`data/rna_cgg/main2/`): 10 runs, 0 failed. Every run streamed exactly the spots in
  `main_runs.txt` (fastq-dump's "Read N spots").
- 09:04 `count --totals-only`: Set A 6,262 observations over 705 sites (73 other); Set C 62,251 over 12,111
  sites (18 other); A0 25 over 12; B 10 over 7.
  - **Powered** (≥ 200 over ≥ 30 sites, and C ≥ 1,000), so the "10 more runs" step does not apply.
- 09:04 The DNA runs (`dna_runs.txt`) are streaming into `data/rna_cgg/dna/`.
- 09:08 **Full count and decision (original method; nothing claimed until the re-check agrees):**
  `count` → `main2/main.obs.tsv`; `decide main2/main.obs.tsv 1`:
  - Set A: D 0.0000 (95% 0.0000–0.0000), lacks 0 of 6,262 observations over 705 sites.
  - Set C: D 0.0000 (95% 0.0000–0.0000), lacks 0 of 62,251 over 12,111 sites.
  - **Decision by the rule: transcripts keep the C.**
- 09:08 Reported extras (original method):
  - A0: 0 of 25 over 12 sites. B: 0 of 10 over 7 sites.
  - D_A divided by the call-rate ratio 0.984 is still 0.
  - With 0 lacks the bootstrap interval is degenerate. One-sided 95% upper bounds:
    - per observation: 0.00048 (6,262 observations);
    - share of Set A sites whose transcripts lack the C in every read: 0.0042 (705 sites, treated as independent).
  - D per run is 0 in all 10 runs (`main2/<run>.totals.txt`, `main2/<run>.obs.tsv`).
    - Set A observations per run: ERR15738679 436, ERR15873922 426, ERR15873926 658, ERR15873927 731,
      ERR15873933 370, ERR15873934 306, ERR15873936 271, ERR15873938 633, ERR15873944 738, SRR37726778 1,693.
    - The per-run counts sum to the pooled totals.
  - D per genome is 0 in all 21 primary genomes with Set A observations. The top four hold:
    - MGYG000302328: 1,424 over 51 sites;
    - MGYG000301353: 1,381 over 93;
    - MGYG000302241: 1,039 over 63;
    - MGYG000302837: 932 over 88.
  - Sites with ≥ 10 observations: A 138 sites, 4,650 observations, 0 lacks; C 1,453 sites, 32,149 observations,
    0 lacks.
  - D_C at walk_has_cgg = 1 sites: 0 of 7,352 over 856 sites. At walk_has_cgg = 0 sites: 0 of 54,899 over 11,255.
  - "other": A 73, of which 67 are at one site.
    - s28605 (MGYG000302328_2:169937, minus strand): all 67 reads come from SRR37726778.
    - Each carries a 1-nt deletion at 169943, the first base of an 8-T run that starts 5 nt past the C. That
      deletion falls inside the ±5 window.
    - So it is a poly-T length variant in that sample's strain, not a change at the C. C "other": 18.
  - Background check, added after the count and not part of the rule:
    - The rate of 1-nt deletions in the flanks (6–20 nt from the run) of the reads the counter observed is
      5.3e-5 per nt at A and 2.1e-5 per nt at C.
    - That rate, times the run lengths observed, predicts 0.39 lacks at A and 1.46 at C from error alone.
    - Observed: 0 and 0 (P(0 | Poisson 1.46) = 0.23).
    - So the counter is not hiding real 1-nt deletions; real ones pass the filters, as the s28605 reads show.
  - Strandedness: from the re-check (my SAM path has no mate numbers).
- 09:14 **The independent re-check agrees** (`recheck/RECHECK.md` §6). It used alignment-free exact
  probes, its own site build, and its own streams of the same 10 runs.
  - **Decision: keeps.** Set A 0 lacks of 5,489 observations over 615 sites; Set C 0 of 51,787 over 10,675.
    D_A = D_C = 0, bootstrap [0, 0]. D = 0 in every run.
  - Per-site join (its file is 1-based, mine 0-based):
    - Set A: 600 sites observed by both. It has 0.98× my depth (median, 130 sites with ≥ 10 of mine).
    - Set C: 10,409 shared sites, 0.96×.
    - Sites seen only by mine: A 105 (423 observations). They come from its probe-uniqueness drop and exact-match
      losses.
    - At s28605 it has 9 keeps against my 10. The 67 poly-T-variant reads cannot match an exact probe.
- 09:14 Strandedness (re-check, read-1 sense share; reported extra):
  - The 8 PRJEB86261 runs sit at 0.43–0.58: unstranded, or DNA. This measure cannot tell the two apart.
  - ERR15738679: 0.23.
  - SRR37726778: 0.0095, i.e. stranded, read 1 antisense (dUTP). This run is clearly RNA.
    - On its own it has 0 lacks: 1,693 Set A observations over 91 sites in my count, 1,581 over 82 in the re-check's.
    - That clears the power thresholds by itself, so the "keeps" result does not rest on runs that might be DNA.
    - Its Set A observations are concentrated (about 19 per site).
  - This subset reading is an extra, not part of the rule.
- 09:15 **Are the unstranded runs RNA? (an extra; criterion written before any ribosomal-protein count)**
  - Among filtered primary reads (the `screen` read filters), CDS reads are 64–84% in the 8 PRJEB86261 runs and
    93% in ERR15738679, against 30% in the stranded SRR37726778.
  - So the CDS-score ranking may have picked rRNA-depleted libraries, or DNA. Strandedness cannot tell the two
    apart in an unstranded library.
  - The test: E = (filtered primary reads in ribosomal-protein CDSs ÷ filtered primary CDS reads) in the RNA run,
    divided by the same share in the paired pig's DNA run (`dna_pairs.txt`).
    - Ribosomal-protein CDSs are those whose GFF product matches `^(30S|50S) ribosomal protein`.
    - The two runs not from PRJEB86261 have no pair; for them the denominator is the share pooled over the 8 DNA runs.
  - Call: E ≥ 2 → RNA-like; E ≤ 1.5 → DNA-like; otherwise unclear.
    - Ribosomal-protein mRNAs are among the most abundant in growing bacteria, while in DNA their share follows gene
      length.
    - As a check of the measure, SRR37726778 (stranded, so RNA) should come out RNA-like, and any DNA run against
      the pooled DNA share should give E near 1.
  - However this comes out, the decision stands. It changes only how the post may describe the 8 runs.

## Follow-up S: stranded libraries only (09:17; written before any of these runs was streamed past the 1 M-spot screen)
- **Why.**
  - The rule's decision (keeps) pooled 10 runs.
  - 8 of them are unstranded, and by ribosomal-protein share 4 of the 5 checked so far look DNA-like (log above).
  - The one stranded run, SRR37726778, carries the RNA reading alone (0 lacks of 1,693 over 91 sites), but it is one
    sample.
  - This follow-up asks the same question on other stranded libraries, new samples only.
- **Stage 1, strandedness probe.**
  - Candidates: every run in `screen_scores_cds_all.tsv` with score_cds ≥ 400 that is neither in PRJEB86261 nor
    in the main set. That is 13 runs: SRR37726779, SRR37726777, SRR37726776, SRR37726773, SRR37726775, SRR29669627,
    SRR29669633, SRR29669628, SRR38252817, ERR15738683, ERR15738673, ERR15738675, ERR15738674. (Count corrected at 09:17 from 12; the list itself is unchanged.)
  - Stream: the first 2 M spots, `fastq-dump --split-spot --readids`, the same minimap2 and filters.
  - Read-1 sense share = among filtered primary CDS reads (the `screen` definition) from read 1, the share on the
    CDS's sense strand.
  - Call: **stranded** if the share is ≤ 0.10 or ≥ 0.90 over ≥ 200 such reads.
  - The ribosomal-protein share is recorded too.
- **Stage 2, count.**
  - Each stranded run: the first 20 M spots (or all, if fewer), same `rna_stream.py` path, same `sites.tsv`, same
    `count`.
  - Decision: the rule's Decision section, applied to the stranded new runs pooled (SRR37726778 excluded), with the
    same thresholds and seed 1.
  - Then the same with SRR37726778 added, reported.
  - D per run and per study are reported.
- If no candidate is stranded, or Stage 2 is not powered, that is reported as is, and the post's RNA sentence rests on
  SRR37726778 alone, with its numbers.

## DNA check: result (09:20)
- All 8 DNA runs streamed exactly 20,000,000 spots (fastq-dump's "Read N spots"), 0 failed parts.
- Counted with the same `sites.tsv` and `count` → `results/rna_cgg/dna/dna.obs.tsv`; per run
  `dna/<run>.obs.tsv`.
- RNA said "keeps", so by the DNA rule this is reported only.
  - **DNA Set A: D 0.0018 (95% 0.0000–0.0058), 11 lacks of 6,199 observations over 948 sites.**
  - **DNA Set C: 1 of 113,117 over 15,925 sites.** A0: 0 of 34; B: 0 of 6.
- **All 11 Set A lacks are at one site, s25050** (MGYG000302139_12:43185, minus strand, CDS MGYG000302139_00958).
  - There: 10 keeps, 11 lacks, from 4 pigs.
    - ERR15873817 (J015): 5 keeps, 1 lack. ERR15873832 (J018): 0, 6. ERR15873836 (J023): 0, 1.
      ERR15873889 (J017): 1, 3.
  - Some strains in these pigs have a genome without this C. The deletion of a skipped C leaves the protein as it
    was.
  - RNA at s25050: 3 keeps, 0 lacks (too few to compare).
  - The one Set C lack is a single read at s32307 (ERR15873827).
- Per-site comparison (the rule's always-reported table) at sites with ≥ 10 observations in both:
  - Set A: 65 sites; DNA 0 lacks of 2,467, RNA (10 runs) 0 of 2,214.
  - Set C: 1,052 sites; DNA 0 of 31,392, RNA 0 of 23,757.
  - On the 8 paired runs only: A 51 sites, DNA 0 of 2,280, RNA 0 of 1,765; C 714 sites, DNA 0 of 26,408,
    RNA 0 of 17,356.
  - The 8 paired RNA runs alone: Set A 0 lacks of 4,133.
- **Ribosomal-protein check (criterion above), all pairs:**
  - Share of CDS reads in ribosomal-protein CDSs: DNA runs 0.0181–0.0222, pooled 0.02035 (21,349 / 1,049,234).
  - E (RNA ÷ paired DNA):
    | pig | RNA / DNA run | E | call |
    |---|---|---|---|
    | J023 | ERR15873938 / ERR15873836 | 1.72 | unclear |
    | J015 | ERR15873922 / ERR15873817 | 1.28 | DNA-like |
    | J022 | ERR15873927 / ERR15873827 | 1.42 | DNA-like |
    | J003 | ERR15873933 / ERR15873839 | 1.23 | DNA-like |
    | J007 | ERR15873934 / ERR15873890 | 1.29 | DNA-like |
    | J018 | ERR15873936 / ERR15873832 | 1.54 | unclear |
    | J017 | ERR15873944 / ERR15873889 | 1.86 | unclear |
    | J001 | ERR15873926 / ERR15873891 | 1.57 | unclear |
  - Against the pooled DNA share: ERR15738679 1.12 (DNA-like); SRR37726778 2.82 (RNA-like, as required of the
    measure).
  - Each DNA run against the pooled share: 0.89–1.09, near 1 as required.
  - **None of the 9 runs other than SRR37726778 is RNA-like by this measure; 5 are DNA-like.** (8 are unstranded;
    ERR15738679 has a read-1 sense share of 0.23.)
    - They are probably mixes: some mRNA, plus genomic DNA that survived the RNA prep.
    - So their 0 lacks says little about mRNA beyond what DNA already says.
  - The RNA reading rests on the stranded run: 0 of 1,693 over 91 sites. Follow-up S tests other stranded runs.
- 09:20 Follow-up S stage 1 tool: `src/rna_strand.py`.
  - Known-answer controls, run alongside the 13 candidates: SRR37726778 must come out stranded with share ≈ 0.01, and
    ERR15873938 must come out unstranded (≈ 0.5). These are the re-check's values from its own counter.
  - If either control fails, the stage-1 calls are not used.
- 09:21 SRR37726778 Set A by genome:
  - MGYG000302328: 1,418 observations over 49 sites;
  - MGYG000302241: 226 over 21;
  - 8 other genomes: 49 over 21.
  - Per site: median 3, max 338.
  - So the stranded evidence is mostly two genomes.
- 09:30 **Follow-up S stage 1 done** (`strand_probe.tsv`; 2 M spots streamed in each of the 15 runs, 0 failed).
  - **Controls pass** (both get the right call). SRR37726778: read-1 sense share 0.030 (1,998 reads), stranded.
    ERR15873938: 0.522 (3,063), unstranded.
    - SRR37726778's 0.030 is higher than the re-check's 0.0095. The read sets differ: the re-check counted reads
      matching site probes over the whole run; this counts every filtered CDS read in the first 2 M spots, which
      includes antisense transcription and reads at overlapping or misannotated genes. Both are far below 0.10.
  - In a dUTP library genomic DNA gives read-1 sense and antisense reads equally, so DNA is at most twice the sense
    share of CDS reads: about 6% in SRR37726778 and at most 9% in any of the stranded runs below.
  - **Stranded (8):** SRR37726779, SRR37726777, SRR37726776, SRR37726773, SRR37726775 (PRJNA1441222, the study of
    SRR37726778; 0.020–0.045) and SRR29669627, SRR29669633, SRR29669628 (PRJNA1129603; 0.017–0.046). All are
    read-1 antisense (dUTP-type), over 416 to 1,770 reads each.
  - **Unstranded (5):** SRR38252817 (0.70), ERR15738683, ERR15738673, ERR15738675, ERR15738674 (0.18–0.23).
    The four ERR runs sit between 0.10 and 0.5. That may be a partly stranded protocol, or a DNA share on top of a
    stranded library; by the rule they are unstranded and are not used.
  - Ribosomal-protein share of CDS reads: 0.023–0.066 in the stranded runs, 0.022–0.047 in the others.
- 09:31 **Stage 2 started:** the 8 stranded runs, first 20 M spots each (all have 33.8–106.0 M spots, SRA
  runinfo), 4 runs at a time, 2 streams each, into `data/rna_cgg/strandS/` (`followS_runs.txt`).
- 09:34 Re-check of stage 2 started: a fresh agent runs the re-check's alignment-free counter (`recheck/`
  code, its own `sites_recheck.tsv`, its own streams) on the same 8 runs, first 20 M spots each, into
  `recheck_S/`. It is told not to read the minimap2 counts or this log past the Follow-up S section.
- 09:54 Stage 2 streamed (`data/rna_cgg/strandS/`): 8 runs, 0 failed, each exactly 20,000,000 spots (fastq-dump's
  "Read N spots", both parts).
- 09:55 **Follow-up S result (original method).** `count` → `followS/S8.obs.tsv`; `decide followS/S8.obs.tsv 1`
  (`followS/S8.decision.txt`):
  - Set A: D 0.0000 (95% 0.0000–0.0000), **lacks 0 of 7,498 observations over 286 sites.**
  - Set C: D 0.0000 (95% 0.0000–0.0000), lacks 0 of 17,358 over 2,414 sites.
  - A0: 0 of 970 over 17. B: 0 of 28 over 14.
  - Powered (≥ 200 over ≥ 30 sites; C ≥ 1,000). **Decision by the rule: transcripts keep the C.**
  - One-sided 95% upper bounds: 0.00040 per observation; 0.0104 for the share of Set A sites whose transcripts lack
    the C in every read (286 sites, treated as independent).
  - **D per run is 0 in all 8** (Set A observations / sites): SRR37726779 2,832 / 120, SRR37726777 1,264 / 95,
    SRR37726776 439 / 53, SRR37726773 865 / 70, SRR37726775 693 / 87, SRR29669627 721 / 58, SRR29669633 329 / 57,
    SRR29669628 355 / 47 (`followS/<run>.obs.tsv`).
  - **Per study:** PRJNA1441222 (5 runs) A 0 of 6,093 over 193 sites, C 0 of 11,636 over 1,802; PRJNA1129603
    (3 runs) A 0 of 1,405 over 117, C 0 of 5,722 over 742.
  - **By genome** (Set A): 19 genomes. MGYG000302241 4,405 over 42 sites; MGYG000302328 1,470 over 60;
    MGYG000301895 1,000 over 42 (all from PRJNA1129603); MGYG000299987 208 over 28; 15 others 415.
    So the stranded evidence now covers 4 genomes with ≥ 200 observations, against 2 before.
  - Per site: median 2, max 1,756 (s27352, MGYG000302241); 76 sites with ≥ 10. Without s27352, 0 of 5,742.
  - Set A "other": 46, of which 45 at s28605, the poly-T length variant already explained above.
- **With SRR37726778 added** (reported): `followS/S9.obs.tsv`: Set A 0 of 9,191 over 304 sites; Set C 0 of 22,389
  over 2,591; decision keeps. Upper bounds 0.00033 per observation, 0.0098 per site.
- The re-check (`recheck_S/`) is still counting; nothing is claimed until it reports.
- 09:55 Samples (ENA filereport): the 9 stranded runs are 9 different BioSamples. PRJNA1441222: pig colonic
  digesta, 6 samples (SAMN56646515–21, SRR37726778 included), NovaSeq 6000. PRJNA1129603: Yorkshire pigs YK1-1, YK1-2,
  YK1-6, collected 2021-10-05, MGISEQ-2000RS. Both list library source METATRANSCRIPTOMIC.
- 10:42 **Re-check of stage 2 (`recheck_S/RECHECK_S.md`; fresh agent, alignment-free exact probes, its own streams of
  the same 8 runs, spots 1–20 M, exact-spot check passed): decision keeps, powered.**
  - Set A: 3 lacks of 6,731 observations over 251 sites, D_A 0.00045 (95% 0–0.0018). Set C: 0 of 14,577 over
    2,084, D_C 0.
  - Per study: PRJNA1441222 A 2 of 5,412 (165 sites); PRJNA1129603 A 1 of 1,319 (106 sites). C 0 in both.
  - With SRR37726778 added: A 3 of 8,312 over 271 sites, D_A 0.00036 (0–0.0014); C 0 of 19,045. Its SRR37726778
    counts match the first re-check exactly.
  - All 8 runs read-1 antisense (A 0–0.022, C 0.013–0.097).
- **The 3 lacks, all at one site** (re-check MGYG000299934:MGYG000299934_7:8068:+ = my s74, in ffh, at the P-loop
  GKTT; no keeps there, one lack each in 3 runs from both studies, read 1 antisense in all three).
  - My pipeline has no MAPQ ≥ 20 alignment within 7,900–8,070 of that contig in any stage-2 run, so s74 has no
    observation at all in `followS/`. Synthetic 151-nt reads at s74 map at MAPQ 60, with the C (151M) and without
    it (75M1D76M), so the 1-nt deletion itself is not what hides them. The real reads must differ from
    MGYG000299934 elsewhere along their length.
  - The 40-nt lack probe occurs exactly in the DNA of 6 of the 278 TANB77 GTDB species representatives on disk
    (`tanb77_noreader_reps.tsv` with a local .fna): GCA_013333035.2, GCA_934167965.1 (reverse strand), GCA_934167095.1,
    GCA_947363765.1, GCA_900772235.1, GCA_948919095.1. The keep probe occurs in none of them.
  - So these 3 reads are most likely transcripts of relatives whose ffh has no C at this point (a site that "comes
    and goes", as in the turnover test), not transcripts that lost it. The re-check's probe-uniqueness filter looked
    only at the 24 genomes, so it could not catch this. Without the site: 0 of 6,728.
  - The re-check's other notes: one site (my s27352) holds about a quarter of Set A; Set C has more read-1-sense
    reads (0.033) than Set A, concentrated at 20 sites, all keeps (antisense transcription at a few loci).
- **Follow-up S conclusion:** both methods give "keeps" on the 8 new stranded runs. With SRR37726778: 9 samples,
  2 studies. The post's RNA line rests on these.

## Fact-check finding and Amendment 5 (10:54; post hoc: written after the fact-checker reported subset counts)
- **Finding (fact-check agent, post_draft/rna):** Set A is every CGG a skip walk meets from the start of a catalogue
  (Prodigal) CDS, so it also holds CGGs that the gene reads **in frame** (the 262-of-1,762 class of the main result,
  and spurious or read-through CDS). Most Set A coverage is there: in S9 the three most-covered sites (1,851, 503 and
  390 read pairs) are CGGs at codons 55, 12 and 28 of "hypothetical protein" CDS of 134, 230 and 204 codons, where
  the skip would end the protein 9, 3 and 1 codons later. So the S8/S9 Set A totals mostly do not bear on skip sites.
  The decision rule was applied as written; the problem is what Set A measures.
- **Pre-registered subset that does bear on it:** Set B (Set A sites of MGYG000302328 within 3 nt of a clean +1
  break, reported only): S9 **0 lacks of 37 over 17 sites**.
- **Amendment 5, post hoc subset (the fact-checker's definition, adopted unchanged):** for each Set A site, from the
  C's codon on,
  - **skip-like:** the original frame reaches a stop within 30 codons after the CGG and the skipped (+1) frame runs
    on for at least 60 codons;
  - **in-frame-like:** the reverse; **unclear:** anything else.
  - Its check on known skips: 27 of the 30 Set B sites come out skip-like, 3 unclear, 0 in-frame-like.
  - The fact-checker's S9 counts (its script, scratchpad skipcheck2.py): skip-like 304 read pairs over 68 sites in 9
    genomes, 0 lacks (minimap2); re-check 278 keeps, 3 lacks (all at the ffh site s74). In-frame-like 6,282 over
    166 sites, unclear 2,602 over 68; 0 lacks in both.
  - To do on resume (paused now): implement the subset in `src/rna_cgg.py` myself, recompute it for both
    pipelines and S8/S9, and give D with the site bootstrap. Also check whether the top in-frame-like CDS sit in
    rRNA or other ncRNA (descriptive).
- **Consequence for the post:** the claim rests on the skip-like subset (about 300 read pairs, clustered: 132 at one
  site) plus Set B, not on 9,191. Upper bound on lacks at skip-like sites is about 1% per read pair (0 of 304,
  ignoring clustering). Headline and conclusion are to be hedged ("most likely during translation").

### Amendment 5 result (12:01; `src/rna_skiplike.py`, outputs in `skiplike/`)
- My own code (written without the fact-checker's script, then set to its conventions, which are the adopted
  definition: "within 30" is <= 30 codons before the stop, and the +1 frame skips again at any later CGG).
  - All Set A sites: skip-like 1,330, in-frame-like 782, unclear 759. Set B: 27 skip-like, 3 unclear (as the
    fact-checker found).
  - My first, stricter reading (< 30; +1 frame read plainly) gives S9 skip-like 280 over 61 sites, 0 lacks
    (`skiplike_strict/`). The conclusion does not depend on the convention.
- **Skip-like subset, decision rule of the Decision section applied to it (`skiplike/decide.txt`):**

  | | original (minimap2) | re-check (exact probes) |
  |---|---|---|
  | S8 (8 new stranded runs) | 0 of 280 over 63 sites | 3 of 259 over 60 sites, D 0.012 (0-0.061) |
  | S9 (+ SRR37726778) | **0 of 305 over 69 sites** | **3 of 281 over 66 sites, D 0.011 (0-0.054)** |
  | decision | keeps | keeps |

  - The S9 count matches the fact-checker's in-frame-like (6,282 / 166) and re-check skip-like (281) counts exactly.
    Its 304 over 68 is 1 short of mine because its script puts a site whose +1 frame reaches the contig end after
    more than 60 codons in a separate class.
  - The 3 re-check lacks are the s74 (ffh) site, explained by relatives' DNA (Log 10:42).
  - S9 skip-like reads by genome: MGYG000302241 156, 302328 79, 299987 40, 301895 17, five others 13. One site
    (s26997, a "hypothetical protein") holds 132; the other 68 sites hold 173 (median 2 per site). Genes include
    ackA, ftsW, dapA, polC, sodA, rny, asnO, mreB, relA.
  - DNA at skip-like sites: read-1 sense 6 of 281 in the re-check, so DNA is at most about 4% there (about 1% over all
    Set A sites, 41 of 8,312).
  - Set B (pre-registered, reported only): 17 observed sites, 37 read pairs, 0 lacks (15 skip-like with 30, 2 unclear
    with 7).
  - In-frame-like: 0 of 6,282 over 166 sites; unclear: 0 of 2,604 over 69 (both pipelines 0 lacks).
- **Reading.** At CGGs where the gene needs the skip to be made in full, every read pair from the stranded libraries
  keeps the C (305; 278 in the re-check besides the 3 explained ones). If 1 transcript in 50 lacked it, about 6 of 305
  would show it (0.98^305 = 0.002). The subset was defined after the count, which the post says.
