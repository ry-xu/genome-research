# Pre-registration: sense-to-stop screen of GTDB R232 (Basil)

Written 2026-11-10 18:02 EST, before any GTDB genome is scored. GTDB profiles saved so far have been used only to
build the Arm A null (the per-GC-bin medians below); no GTDB genome has had its own codons scored. Sections marked
PENDING were filled from control runs only, before any GTDB genome other than the 40 Arm B null genomes goes through
Arm B. Any later change to a frozen threshold is logged at the bottom with its reason and date, and results are
reported under both the frozen and the changed rule.

## Data
- All GTDB R232 species representatives (bacteria and archaea), streamed from
  `gtdb_proteins_nt_reps_r232.tar.gz`; GTDB's own Prodigal nucleotide gene calls (`src/gtdb_scan.py`).
- Calling table per genome: 4 if fewer than 0.2% of complete genes end at TGA, else 11 (GTDB's calls use one of the
  two). The codons tested are the sense codons of that table (61, or 62 with TGA under table 4).
- A table-4 genome in which TGA is less than half of all TGA + TGG codons is a *suspect table-4 call* (GTDB called a
  standard-code genome with table 4; the audit found 4 of the first 75 table-4 genomes like this, GC 0.48-0.58, TGA
  25-34% of Trp codons against a median of 90%). Its flags belong to the chimeric-call artifact class; it still goes
  to Arm B, called with table 11.
- Genomes with fewer than 300 complete genes are not scored by Arm A (they can still enter Arm B through the
  stratified sample).

## Arm A (screen)
Definitions are those in `src/score.py` and `src/codonpos.py` (n3 = count at d <= 100 codons from the gene end; nI =
count at d >= 151, first 20 codons excluded; exp3 = nI * (N3/NI) * rho; excess = n3 - exp3; lr = log2((n3+1) /
(exp3+1)); log10p = Poisson upper tail with expectation floored at 2; depl = interior rate / GC-bin median interior
rate).

**Null.** rho and the interior-rate medians come from every scanned GTDB genome with calling table 11, >= 300
complete genes and >= 1,000 sense codons in each region, per GC bin of 0.025 from 0.15 to 0.80 (bins with < 5
genomes take the nearest filled bin). Real stops that are sense under the calling table (TGA under table 4) are scored with rho = 1. The tested
genomes are part of the null; a median over ~10^5 genomes is not moved by a handful of positives.

**Flag rule (frozen).** A (genome, codon) pair is flagged when
- log10p <= -10 and excess >= 10, and
- T1: lr >= 1.5, or T2: lr >= 1.0 and depl <= 0.2.

**Calibration behind it** (all before any GTDB genome was scored). The rule was chosen with a null from the first
8,000 scanned genomes, which are all archaea (the tar streams archaea first). The numbers below are the re-check
(`scratch/recheck_armA.py`, 2026-11-10 18:25) with a null from the first 34,000 genomes (33,892 in the null, all
archaea plus ~24,000 bacteria) and the 0.15-0.80 GC bins; target recall was the same under both nulls.
- Natural controls (`data/controls.tsv`). Genes called with a forced table so a real stop looks like a sense codon
  to the method:
  - ordinary genomes under table 4: TGA flagged in 29/36; under table 16: TAG flagged in 33/36;
  - *Methanosarcina acetivorans* (TAG mostly a stop) under table 16: TAG flagged; *Mycoplasma capricolum* and
    *M. genitalium* (TGA = Trp, TAG a stop) under table 16: TAG flagged in 2/2.
  - Genomes that really read the codon as sense must not flag it: TAG readers under table 16, 0/10; TGA readers under
    table 4, 0/8; code-34 genomes under table 16, 0/2.
  - Under GTDB's own rule (the condition of the real scan): 1/36 ordinary genomes flagged anything (*Bordetella
    pertussis*, CTA and AGA, an IS481-rich genome), 0/10 TAG readers, 0/2 code-34 genomes. For the 8 TGA readers
    the rule picks table 4, so their GTDB-rule calls are the forced table-4 calls above (0/8).
  - Non-target codons also flagged under a forced table: 1/36 genomes under table 16 (*Thermus*: TGT, CAT, AGT);
    5/36 under table 4 (*E. coli* AGG; *Thermoplasma* TGT; *Aquifex* CGA, CGG; and *S. coelicolor* and *M. luteus*,
    high-GC genomes whose table-4 calls are chimeric, 7 and 12 codons). *M. capricolum* under table 16 also flags TGG, its other Trp
    codon (its TGA Trp codons end the forced calls).
- Injections (`src/inject.py`, `src/calibrate.py`: a fraction s of real gene ends recoded to X, and a fraction q of
  X's other occurrences recoded to synonyms; q = 1 is a clean reassignment): recall of T1 or T2 at q = 1 is 0.46 /
  0.68 / 0.86 / 0.89 / 0.98 / 1.00 at s = 0.005 / 0.01 / 0.02 / 0.03 / 0.05 / 0.1; at q = 0.9, 0.23 / 0.25 / 0.37 at
  s = 0.03 / 0.05 / 0.1; at q <= 0.5, <= 0.11. Injection null panel: 2/78 runs flag anything (both *B. pertussis*,
  4 codons); 39/1,166 injected runs also flag a non-target codon.
- The final null (all GTDB) changes these slightly; the calibration is recomputed with it and reported, and the rule
  is not changed.
- Arm A therefore covers clean reassignments; an ambiguous codon (read as sense part of the time) is left to the
  stratified Arm B sample.

**Known artifact classes**, written down before scoring so they are not discovered post hoc: IS-element-rich
genomes (*B. pertussis*); AT-rich rare codons piling up in run-on tails of genes called past a frameshift or
pseudogene; chimeric calls under the wrong table. A flag is a reason to run Arm B, never a result.

## From Arm A to Arm B (budget)
- If <= 2,500 genomes carry a flag, all go to Arm B.
- Otherwise, in this order until 2,500: (1) codons flagged in >= 2 species of the same GTDB genus, or >= 3 of the
  same family (clade recurrence); (2) the rest by most negative log10p. At most 20 genomes per genus at each step.
- Stratified sample, independent of Arm A: one genome per GTDB family (largest number of complete genes; ties by
  accession), up to 1,500 families chosen at random with seed 20261110; Arm B scored on every sense codon. This is
  the only coverage of ambiguous reassignments.
- Arm C extra entry (`src/rf_motifs.py`): a genome goes to Arm B when, for a release-factor class it carries, no
  copy has the canonical motif at the reference positions (RF1 P.T at E. coli RF1 188-190; RF2 SP[FY] at E. coli
  RF2 205-207; aRF1 NIKS at M. jannaschii 56-59 and Y.C...F at 121-127). Copies aligned over < 60% of the reference,
  or with a gap in the motif, are missing data, not a motif change.

## Pre-registered lead: Coriobacteriales RF2 (added 2026-11-10 18:09 EST, before any codon data of these genomes was looked at)
The motif screen on the first ~18,000 scanned genomes found RF2 with SPT in place of SPF (E. coli RF2
`VRKSPFDSG`, Coriobacteriales `VRISPTDDK`) in every Coriobacteriales genome with an RF2 hit (200 genomes, all five
families scanned so far), and no other RF2 copy in them. In the tripeptide-anticodon model the third residue reads the
third codon base (RF1 T: A or G; RF2 F: A only), so an S-P-T RF2 could read UGG (Trp) as well as UGA. The test:
- every Coriobacteriales GTDB representative goes through Arm B (the Arm C entry) with the frozen Arm B rule; TGG
  passing in any of them is a candidate as usual;
- clade contrast, whatever the per-genome result: TGG ratio (Arm B) and TGG lr (Arm A) in Coriobacteriales vs the
  other Actinomycetota orders (one Arm B genome per family, the stratified sample), one-sided Mann-Whitney at
  p < 0.01, effect reported as the median difference;
  - unit (clarified 19:28, before any codon data of these genomes was looked at): one genome per GTDB family on both
    sides, chosen by the stratified-sample rule (most complete genes, ties by accession). Coriobacteriales: all its
    families. Other side: every other Actinomycetota family (all of them, not only those drawn into the random
    1,500; those not drawn are added to Arm B for this). Families whose chosen genome has < 20 Arm B events are
    left out of the Arm B comparison and counted. Reported, not tested: the same with every Coriobacteriales
    genome, and a split by RF2 motif (SPT vs SP[FY] carriers) across Actinomycetota;
- a negative result is reported as such (the motif change then does not make TGG a stop at a detectable level).

## Arm B (test)
Definitions in `src/term_codon.py`: genes called with the genome's table (Pyrodigal single mode); DIAMOND blastp
vs Swiss-Prot; an event is a gene whose best hit covers >= 60% of the Swiss-Prot protein, ends within 3 residues of
its C-terminus, and whose protein continues >= 5 residues past the alignment end; O = events with X in the window
e*-1 .. e*+3; O_fam = distinct families among them (events sharing a Swiss-Prot mnemonic or with proteins >= 50%
identical over >= 50% of both, joined; Changes 19:59); hits >= 95% identical are dropped first; E from the larger of X's tail and interior
frequency in the same genes; ratio = O / E; log10p Poisson upper tail.

**Pass rule (frozen 2026-11-10 19:20 EST, from control runs only; `scratch/armB_eval.py`):** a (genome, codon) pair
passes when O_fam >= 10, ratio >= 3 and log10p <= -10. The codon need not rank first; its rank is reported.

Calibration (`results/armB2/`, fixed code, Changes 19:59; the first run, `results/armB/`, is superseded):
- Natural positives (ordinary genomes called with table 4 or 16, so TGA or TAG is a "sense" codon that really ends
  genes; *M. acetivorans* under table 16): 69/73 pass for the target codon, and no other codon passes in any of
  them. The four misses: *M. luteus* under both tables (12 events under table 4, O = 6; O = 2 under table 16),
  *C. difficile* under table 4 (O_fam 17, ratio 6.3, log10p -8.35) and *M. acetivorans* TAG (O_fam 9, ratio 9.5,
  log10p -6.16; its TAG is partly read as pyrrolysine).
- Negatives: TAG readers under table 16 and under GTDB's rule (10 + 10), TGA readers under table 4 (8), code-34
  genomes (2 + 2), *B. pertussis* (1), and 40 random high-quality GTDB R232 representatives (`data/armB_null.tsv`,
  seed 20261110): 0/73 pass any codon. The most extreme of the 4,473 null tests is log10p -8.2 (CCC, ratio 2.45).
- The null tail is heavier than Poisson: tests at log10p <= -2 / -3 / -4 / -5 / -6 / -7 number 73 / 25 / 10 / 4 / 2
  / 1, falling ~0.4 log10 units per unit. Extrapolated, a run of ~4,000 genomes x 61 codons gives on the order of
  1-3 false passes; the candidate checks below (clade, raw reads, independent re-check) are what remove them.
- Known artifact class (written 19:26, before any GTDB genome other than the 40 null genomes went through Arm B):
  C-terminal composition. The window covers the homolog's last residues, so amino acids common at C-termini raise
  their codons there; lysine codons are 8 of the 23 null tests at log10p <= -3, with ratios up to 2.98 (AAA)
  (fixed code: 8 of 25, AAA up to 2.98, AAG up to 3.96). A
  composition artifact raises all codons of the amino acid; a stop raises X alone. Hence check 2 below.
- Injections (12 per cell; s = fraction of gene ends recoded to X, q = fraction of X's other occurrences recoded
  to synonyms): s = 0.1: 7, 7 and 11 of 12 at q = 0 / 0.5 / 1; s = 0.03: 0, 2 and 3 of 12; s = 0.01: 0, 1 and 1 of
  12 (the first run, inflated by self-hits: 10, 11, 12 / 2, 2, 5 / 2, 2, 2). The 6 passes at s <= 0.03 come from 5
  of the 12 host genomes (*B. subtilis* 3 of 6). No injected run passes a non-target codon. Arm B therefore finds a
  codon that ends ~10% of genes, best when X is no longer read as sense elsewhere; Arm A finds rarer clean stops
  (s ~ 0.02-0.03), which Arm B mostly cannot.

## From an Arm B pass to a result
A (genome, codon) pair that passes Arm B is a candidate. A candidate becomes a result only if all of:
1. **Clade:** the same codon passes Arm B in at least one other assembly of the same species or genus (any GTDB
   genome, not only representatives). A single-assembly signal is reported as unconfirmed.
2. **Synonym contrast** (added 19:26, see Changes): X's ratio >= 3 x the pooled ratio (sum O / sum E) of X's
   synonymous codons in the same genome. Not applicable to ATG and TGG under table 11 (no synonyms; reported). On
   the controls (fixed code): natural positive passes 69/69 at >= 4.7 x; injected passes 32/32 at >= 3.8 x; the
   null's composition rows ~2 x (AAA 2.98 vs AAG 1.46); null rows with a single-codon excess can exceed 3 x (ATA 16.8 x
   at log10p -5.4, GGA 11.0 x at -6.6), so check 2 removes composition artifacts but not codon-specific noise.
3. **Raw reads:** at >= 10 Arm B event sites, reads (SRA) mapped to the assembly carry X in >= 90% of reads at
   depth >= 5.
4. **Independent re-check:** a fresh subagent, with a different method (tblastn of conserved single-copy proteins
   on the raw contigs, no gene calls), finds the same homolog-end pile-up at X.
5. **Machinery** (reported, not required): release-factor motif, tRNA readers of X (seed-archive tRNA screen),
   codon usage of X across the clade.

## Phages and viruses (written 2026-11-10 19:31 EST, before any phage or virus genome was scored)
A virus is translated by its host's ribosomes and release factors. It can end its genes at a sense codon only if
its host already does (then GTDB finds the host) or it brings its own release factor; some recoded gut phages
encode RF1/RF2 (Borges et al. 2022, Nat Microbiol). The viral search therefore starts from the machinery. Most phage
genomes are below 300 genes, so Arm A does not apply, and Arm B is pooled.

**Data.** INPHARED 14Apr2025 (complete phage genomes from GenBank, with host genus and ICTV taxonomy); MGV (189,680
human-gut virus genomes >= 50% complete, 54,118 vOTUs, its own Prodigal proteins and stop-readthrough calls); NCBI
RefSeq viral release (eukaryotic viruses, including giant viruses that encode translation factors).

**V-C, machinery (all viruses).** Every 280-480 aa protein with GGQ (INPHARED: Pyrodigal meta calls; MGV: its
proteins; RefSeq: its proteins) is searched against `markers.dmnd` exactly as in the GTDB scan (e-value 1e-20,
best hit RF1/RF2/aRF1), and motifs are read by `src/rf_motifs.py`. Every viral release factor is reported with its
motif and host. A viral RF without the canonical motif (Arm C's definition, same missing-data rule) is a lead.

**V-B, pooled Arm B (viruses that carry an RF).** Pool = one genome per vOTU (MGV) or species (INPHARED, RefSeq)
among all viral genomes whose RF is >= 50% identical to the lead RF (DIAMOND), i.e. the clade that shares the
machinery; a giant-virus genome with >= 300 complete genes is also tested on its own with the GTDB rule. Genes are
called on the whole pool at once (`term_codon.py` single mode trained on the concatenated pool, table 11), and
events are kept once per Swiss-Prot family per pool (the first by accession and position), so near-identical genomes
cannot repeat an event; O = O_fam. Every RF-carrying pool is tested on every sense codon, plus the synonym contrast.

**V-B calibration, before any RF pool is scored** (pools of the same size distribution as the RF pools, one genome
per vOTU/species; the pass rule starts as the GTDB rule and is frozen from these alone):
- positives: pools of ordinary phages (INPHARED genera with no readthrough call) called with forced table 4 (TGA)
  and 16 (TAG) must pass the target;
- negatives: MGV pools of genomes that MGV calls as reading TAG (or TGA) through, called with forced table 16
  (or 4), must pass nothing; 30 random pools under table 11 must pass nothing;
- injections into pools of ordinary phages: X recoded at 3% and 10% of gene ends, q = 1 and q = 0.
A V-B pass then goes through the same result checks as a GTDB candidate (clade = a second pool or vOTU passing the
same codon; raw reads where the virus came from sequenced reads; independent re-check; synonym contrast).

**Reported, not tested:** phages of any GTDB candidate's host genus; viruses whose RF was not found by the
e-value 1e-20 search are invisible to V-C (a stated limit).

## Changes after freezing
- 2026-11-10 18:25 EST, after the independent audit of the Arm A code and before any GTDB genome was scored: GC
  bins extended from 0.20-0.70+ to 0.15-0.80 (the old top bin pooled GC 0.675-0.754, where AT-rich codon interior
  rates vary ~9-fold and made T2 easy for AT-rich codons in the highest-GC genomes); `gc_bin` now rejects values that
  are not GC fractions; suspect table-4 calls labelled; calibration text corrected (non-target counts, injection null
  2/78 runs not 4/78, archaeal calibration null). Flag thresholds unchanged.
- 2026-11-10 19:26 EST, before any GTDB genome other than the 40 Arm B null genomes went through Arm B: C-terminal
  composition written down as an Arm B artifact class, and the synonym contrast added as check 2 for a candidate to
  become a result (prompted by Clover's 18:58 note that her nuclear homology-end null was anti-conservative for
  K/N/Q/D/E codons; `scratch/armB_syn.py`). The Arm B pass rule is unchanged.
- 2026-11-10 19:28 EST, before any codon data of Coriobacteriales was looked at: the unit of the Coriobacteriales
  contrast written down (one genome per family on both sides; all Actinomycetota families on the other side;
  >= 20 Arm B events). It was ambiguous as first written.
- 2026-11-10 19:31 EST: phage and virus section written (V-C, V-B, its calibration) before any phage or virus genome
  was scored. V-B thresholds are frozen from its calibration and logged here.
- 2026-11-10 19:59 EST, after the independent audit of the Arm B code (`scratch/audit_armB/`) and before any GTDB genome other than
  the 40 Arm B null genomes went through Arm B: (1) hits >= 95% identical are dropped before the best Swiss-Prot hit
  is chosen (top 10 targets kept, was 5), for every genome, because control genomes that are themselves in Swiss-Prot
  matched their own entries (*B. subtilis*: all events; *M. tuberculosis* 61%), which emptied the background and
  inflated injection recall at s = 0.01-0.03 (2 high-self genomes passed 11/12 there vs 4/60 for the other 10);
  (2) families are now unions of events sharing a mnemonic or with proteins >= 50% identical over >= 50% of both
  (numbered paralogs and IS copies such as INSI1/INSI3 had counted as separate families); FAMDEDUP uses the same
  families; (3) a forced table now always means single mode trained with that table, and inputs < 20 kb are
  refused (meta mode below 100 kb had ignored the requested table and labelled the result 11); (4) `inject.recode`
  no longer recodes start codons or codons overlapping another gene's start or stop (it had changed up to 12% of
  starts for X = TTG at q = 1 and turned a few stops into non-X codons). All Arm B controls are rerun with the fixed
  code (`results/armB2/`); the pass rule is re-checked on them alone and any change is logged here. Arm A's injection
  figures above were made with the old `recode` (affected: X = GTG/TTG at q > 0, about 1 in 30 runs, and the small
  overlap effects); that grid is rerun when cores are free and the figures are replaced, flag thresholds unchanged.
- 2026-11-10 20:03 EST, before any phage or virus genome was scored: V-B pools are built by `src/viral_pools.py`. RF clusters =
  single linkage of viral RFs >= 50% identical over >= 50% of both. Calibration genomes come from MGV, not INPHARED
  (INPHARED has no readthrough calls, and the RF pools are mostly MGV gut phages, so MGV matches them better):
  ordinary = Complete/High-quality MGV genomes with no readthrough call; negatives = MGV TAG (or TGA) readthrough
  genomes under forced 16 (or 4); one genome per vOTU, never an RF-carrying vOTU; each calibration pool copies an RF
  pool's shape (same number of genomes, each drawn at random from the 5 unused genomes nearest in length), seed
  20261110. RF pools are split by MGV readthrough call (one pool per cluster and call; the longest genome per unit),
  and a TAG (TGA) pool is called with table 16 (4), as its calibration negatives, instead of 11, since calling a
  TAG-reading phage with 11 cuts its genes at every TAG. All V-B runs use FAMDEDUP=1 and the fixed `term_codon.py` (Changes, 19:59).
- 2026-11-10 20:48 EST, before any GTDB genome other than the 40 Arm B null genomes went through Arm B: Arm B controls
  rerun with the fixed code (`results/armB2/`, evaluated by `scratch/armB_eval2.py`, synonym contrast by
  `scratch/armB_syn.py`). The frozen pass rule (O_fam >= 10, ratio >= 3, log10p <= -10) is kept unchanged: natural
  positives 69/73, negatives and nulls 0/73, no off-target pass in any positive or injected run; injection recall at
  s <= 0.03 dropped (self-hits had inflated it). Calibration text, the O_fam definition and the check-2 figures above
  are updated to the fixed-code numbers; the old numbers are kept in parentheses where they were quoted.
- 2026-11-10 20:51 EST, before any giant-virus genome was looked at: GVMAGs V2 (Vasquez et al. 2025, bioRxiv 10.1101/2025.09.26.678796;
  18,727 Nucleocytoviricota and Mirusviricota genomes, portal.nersc.gov/cfs/nelli/gvmagsV2/) is added to the V-C
  data, since RefSeq holds few giant viruses and most known virus-encoded eRF1s are in metagenome-assembled ones. The
  same screen is used (Pyrodigal meta calls unless the release has proteins; then its proteins). The paper reports only
  stop-to-sense codes in giant viruses (TAA/TAG to Gln, TAG to Tyr), chosen by GVClass among codes that do not add stops;
  its per-genome code (`ttable`) is reported for every RF-carrying genome and used as that genome's table in V-B. A
  giant virus with a non-canonical RF goes to V-B as written above (pool per RF cluster; on its own with the GTDB rule
  if it has >= 300 complete genes).
- 2026-11-10 20:55 EST, after the giant-virus RF screen (`results/gv_rf_motifs.tsv`: 46 RFs, 5 eRF1-type copies non-canonical, 2 with
  the motif missing) and before any codon data of any giant virus was looked at: (1) for viruses an RF is a lead
  per copy (V-C wording), not per genome, since a virus can also use its host's factor; (2) before any lead genome
  is scored, the per-genome GTDB rule is checked on giant viruses: 30 GVMAGs V2 genomes drawn with seed 20261110 from
  those with >= 400 genes and no RF hit, run under table 11 (null: must pass nothing) and forced 4 and 16
  (positives: TGA / TAG must pass); if any null genome passes, the per-genome rule for giant viruses becomes
  log10p <= (most extreme null log10p - 2), logged here before lead genomes are scored; (3) lead genomes are run
  under GVClass's `ttable` (codemeta -> 11), as already written, and also under 11; (4) the lead eRF1 contigs are
  checked for being viral (the contig's other genes' best hits) before any result is claimed.
- 2026-11-10 21:04 EST, after the giant-virus controls of the 20:55 entry and before any codon data of any lead genome
  was looked at: Arm B against Swiss-Prot has no power on giant viruses (`results/gv/`: 30 null genomes, median 3
  events, at most 12, none passes; forced 4 and forced 16 positives 0/30 each), because few giant-virus proteins have
  a Swiss-Prot homolog covering >= 60% of it. Giant-virus genomes are therefore tested against a reference of
  giant-virus proteins instead (`src/term_codon.py` with `ARMB_DB`; code otherwise unchanged, Swiss-Prot runs
  unaffected): the GVMAGs V2 release proteins of all genomes whose `ttable` is code1, code11 or codemeta (standard
  stops), minus the 30 control genomes, the 5 lead genomes and the 2 genomes whose eRF1 motif is missing
  (`data/giantvirus/gvref.dmnd`). Each reference protein is its own family key (no mnemonic), joined by the same
  >= 50% identity / >= 50% coverage rule; hits >= 95% identical are dropped as before. The same 30 null genomes
  (table 11) and 30 positive genomes (forced 4, forced 16) are rerun against it; the pass rule (O_fam >= 10,
  ratio >= 3, log10p <= -10) is kept if no null passes and is otherwise tightened as in the 20:55 entry (2) — logged
  here before any lead is scored. If fewer than 20/30 positives pass in either forced table, the giant-virus test is
  reported as underpowered and no giant-virus lead is called a candidate from Arm B alone.
- 2026-11-10 21:13 EST, after the V-C screen of all three viral sources (`results/viral_rf_motifs.tsv`: 42 RFs in 8
  clusters, all canonical) and before any phage or virus pool was scored: (1) RF-pool genomes from INPHARED and
  RefSeq have no MGV readthrough call, yet one cluster holds the 14 INPHARED Lak megaphages next to MGV genomes
  called TAG. Their call is made the usual way for phages: Pyrodigal single-mode coding density under tables 11, 4
  and 16; the call is TGA (TAG) if density under 4 (16) is >= 10 percentage points above density under 11 (the
  larger gain if both), else none; the genome then joins that call's pool and table as in the 20:03 entry. Computed
  in `src/viral_pools.py clusters` before pools are built. (2) The six RefSeq Marseilleviridae (cluster with
  canonical eRF1s) are giant viruses: besides their V-B pool, each is tested on its own against the giant-virus
  reference under the giant-virus rule once that is frozen (21:04 entry). (3) INPHARED units stay genus or
  accession (as coded since 20:03); near-identical genomes in one pool (the Lak phages) cannot repeat an event
  because events are kept once per family per pool.
- 2026-11-10 21:31 EST, after the V-B calibration against Swiss-Prot and before any RF pool was scored: pooled Arm B
  against Swiss-Prot has no power on phages (`results/vb/`: 30 null pools mean 2.7 events, at most 18; forced 4 and
  16 positives mean 3.2; nothing reaches O_fam >= 10), for the same reason as giant viruses. V-B therefore uses a
  phage reference instead (`ARMB_DB`, as in the 21:04 entry): MGV proteins (MGV's own Prodigal calls, complete genes
  only, partial=00) of one genome per vOTU (the longest), for every vOTU none of whose genomes MGV calls as reading
  TAG or TGA through, minus every vOTU with a genome in any V-B pool (RF or calibration)
  (`data/phage/phageref.dmnd`, built by `src/phage_ref.py`). The six calibration sets are rerun against it; the
  pass rule (O_fam >= 10, ratio >= 3, log10p <= -10) is kept if no null, negative or off-target codon passes, and is
  otherwise tightened to log10p <= (most extreme of those - 2), frozen and logged here before any RF pool is scored.
  If fewer than 20/30 positive pools pass in either forced table, V-B is reported as underpowered and no RF pool is
  called a candidate from V-B alone. The Marseilleviridae pool (rf8, giant viruses) is run against the giant-virus
  reference instead, under the giant-virus rule.
- 2026-11-10 21:34 EST, giant-virus rule frozen, before any lead genome was scored (`results/gv2/`, giant-virus
  reference): 30 null genomes (table 11) pass nothing (events per genome median 22, 6-48; most extreme test AAC,
  log10p -3.95); forced-4 positives pass TGA 19/30, forced-16 positives pass TAG 21/30; no off-target codon passes
  in any positive. The pass rule (O_fam >= 10, ratio >= 3, log10p <= -10) is kept. By the 21:04 criterion (TGA
  19 < 20) the giant-virus test is reported as underpowered: a lead that fails Arm B is not ruled out (a whole-codon
  stop is missed about 1 time in 3 in a 400-2,400-gene genome, more often in smaller ones), and a lead that passes
  is not called a candidate from Arm B alone.
- 2026-11-10 21:56 EST, V-B rule frozen against the phage reference, before any RF pool was scored (`results/vb2/`,
  `scratch/vb_eval.py`): 30 null pools (table 11) pass nothing (events per pool median 6, at most 148); the 30
  negative-16 and 10 negative-4 pools (forced tables on unrecoded genomes) pass nothing; no off-target codon passes in
  any set. The most extreme non-target test is GAG in positive-16 pool pos_11 (O_fam 74, ratio 2.79, log10p -13.49;
  fails the ratio), so the pass rule (O_fam >= 10, ratio >= 3, log10p <= -10) is kept unchanged. Recall: forced-4
  positives pass TGA 18/30, forced-16 positives pass TAG 11/30, injected sense-to-stop pools 7/30 (s = 0.03: 2/15,
  s = 0.1: 5/15). Both positive sets are below 20/30, so V-B is reported as underpowered: an RF pool that passes is
  not called a candidate from V-B alone, and one that fails is not ruled out.
- 2026-11-10 23:01 EST, after the independent audit of the final-run scripts (`src/armA_run.py`, `src/armB_select.py`,
  `src/rf_motifs.py`; `scratch/audit_final/`) and before either script was run on the full scan. Disclosure: to test
  its findings the auditor ran the scripts on early chunks and so saw partial outcome data (70 flagged genomes among
  the first 4,000 scanned, two of them Methanomassiliicoccales flagged at CGA, and median TGG lr among Actinomycetota
  in a subset of chunks); testing the fixes I saw that 120 of the first 8,000 genomes carry a flag (no codons or
  genomes looked at), and the Arm C entry load on the first 21,793 and 150,841 scanned genomes (counts by motif
  state). None of the rules below was chosen by which codon or genome is flagged. Flag thresholds, the Arm B pass
  rule and the Coriobacteriales test are unchanged. Fixes:
  (1) `armA_run.py` reads GTDB metadata from `data/gtdb_meta/` (copies, checked identical, of the symlinked
  stop-codon-census files) and refuses to score unless `results/gtdb/scan.log` has its `done:` line and the profiles
  hold exactly the 199,923 R232 species representatives (189,801 bacterial + 10,122 archaeal, `gtdb_representative`
  = t). Unscanned representatives go to `armA_missing.txt`; `ALLOW_PARTIAL=1` overrides, and any use of it in the
  final run will be logged here with the count.
  (2) `armA_rows.tsv.gz` keeps the TGG row of every scored genome whatever its lr (the Coriobacteriales contrast
  compares TGG lr between families, and the lr >= 0.5 cut kept only a biased few); values written to 6 significant
  digits (were 4).
  (3) `armB_select.py` writes each genome's Arm B table: 11 for a suspect table-4 call (Data section), else the Arm A
  calling table. Arm B is run with that table, not re-derived.
  (4) Budget step 2 ranks every flagged genome not chosen in step 1, including clade-recurrent genomes left out of
  step 1 by the per-genus cap (it had excluded every recurrent genome, which could leave budget unused while
  flagged genomes waited). The cap of 20 per genus counts within each step, so a genus can reach 40 in all.
  (5) `rf_motifs.py` writes `ref_cover` to 6 decimals (was 2: a copy covering 0.5962 of the reference printed as
  0.60 and counted as present, against the < 60% missing-data rule). The three viral motif tables have no row
  printed 0.60, so V-C is unaffected.
  (6) Arm C cap (new). The Arm C entry was unbudgeted: it entered 1,051 of the first 21,793 scanned genomes and
  2,854 of the first 150,841 (motif table made before fix (5); 4 rows at 0.60), likely ~3,500 in all, more than
  Arm B can run alongside the flag budget. It now keeps one genome per
  (GTDB genus, RF class, set of non-canonical motif residues over that class's non-canonical copies), the one with
  most complete genes, ties by accession: 436 and 763 entries on the same sets. A motif change that turns a codon
  into a stop is expected to be shared within a genus, so one test per motif state per genus keeps the question
  each entry asks. All Coriobacteriales representatives already go through Arm B on their own.
  (7) Literal readings kept and written down: a copy with a gap in either aRF1 motif is missing data even if the
  other motif is non-canonical; an X residue counts as non-canonical; steps 1 and 2 order a genome by its most
  negative log10p over all its flagged codons.
  (8) Known, unchanged: `score.poisson_sf_log10` overflows to +inf when n3 is far below its expectation (ratio series
  with k << mu). It cannot touch a flag (log10p <= -10 needs k > mu); such a value can only appear in a TGG row.
  Tests: `scratch/test_armB_select.py` (extended: table column, step 2, Arm C cap and its key), the auditor's edge
  cases on the fixed code (`scratch/audit_final/f4_fixed.py`; E2 now expects 40 from one genus, by (4)), and
  `armA_run.py` on 4 chunks in a scratch directory (refuses without `done:`; with ALLOW_PARTIAL=1, TGG rows for all
  8,000 genomes, no inf or nan written).
- 2026-11-11 00:29 EST, data handling, no rule changed. (1) Six Coriobacteriales representatives have no sequence at NCBI under
  the GTDB accession (Datasets API "no fna"; the FTP directories hold reports only). Each was run from a stand-in
  accepted only if its total length matched GTDB's genome_size exactly and its scaffold count and N count matched
  the NCBI .1 assembly stats (criterion set before download): the .2 versions of GCA_902364775, GCA_902364885 and
  GCA_902373375 (same length and contig count as .1), and MGnify UHGG v2.0.2 MGYG000000156/207/266 for
  RS_GCF_902363845.1, RS_GCF_902364355.1 and RS_GCF_902364945.1 (UHGG MGYG-HGUT-00156/207/266; length, scaffolds and
  N all equal). Map: `data/corio/substitutes.tsv`; rows `results/corio/armB_subst.tsv`. Arm B results seen after the
  run: none passes. All 1,934 Coriobacteriales representatives have now gone through Arm B. Seven GenBank
  representatives were run from their RefSeq twin (GCF_x for GCA_x), as `src/armB_plan.py` does for the GTDB set.
  (2) Two Coriobacteriales genomes first run with table 4 were rerun with table 11 (`results/corio/armB_rerun.tsv`,
  replaces): GCA_030536315.1 (calling table 11) and GCA_045641605.1 (suspect table-4 call). (3) `src/calibrate.py`
  takes `INJECT_OUT`; the Arm A injection grid with the fixed `inject.recode` goes to `results/inject2/` (same
  panel, specs and seeds), scored by `src/inject_eval.py` against the final null; the old grid stays in
  `results/inject/`.
- 2026-11-11 00:57 EST, Arm A injection calibration recomputed (replaces the Arm A section's injection figures; flag rule
  unchanged). Final null (198,932 GTDB genomes with calling table 11 and >= 300 complete genes), frozen flag rule,
  grid rerun with the fixed `inject.recode` (`results/inject2/`, `src/inject_eval.py`; output
  `results/inject2/eval_*.txt`). Recall at q = 1: 0.65 / 0.73 / 0.88 / 0.87 / 0.98 / 1.00 at s = 0.005 / 0.01 / 0.02
  / 0.03 / 0.05 / 0.1 (n 46 / 59 / 43 / 45 / 44 / 58); at q = 0.9: 0.11 / 0.21 / 0.39 at s = 0.03 / 0.05 / 0.1; at
  q <= 0.5: <= 0.09. Null runs: 2/78 flag anything (both *B. pertussis*, CTA and AGA, as before). Injected runs also
  flagging a non-target codon: 33/1,168. The old grid under the final null gives the old figures to within one run
  (0.46 / 0.68 / 0.86 / 0.89 / 1.00 / 1.00; 40/1,166 off-target), so the change is the recode fix, not the null.
  Why s = 0.005 rose from 0.46 to 0.65 (`scratch/inject_diff.py`, `scratch/inject_tail.py`): (a) the fix draws the
  q coin only for X codons it may recode, so almost every realization differs (1,141 of 1,170 injected runs; at s = 0.005, q = 1,
  12 runs turned to recalled and 3 the other way); (b) X codons overlapping a neighbouring gene's start or stop are
  now left as X. At q = 1 they are the only X left near gene ends besides the rewritten stops, mostly the last sense
  codon of a gene in a tight operon, so they add to the end signal (q = 1 runs: X in the last-10-codon bin median 29
  vs 23 with the old code, although fewer stops were rewritten, median 52 vs 60). A real clean reassignment would
  likely keep such an X (it is a stop one codon early), but recall at the smallest s depends partly on it, so the
  s = 0.005 figure is reported with this caveat.
- 2026-11-11 02:27 EST, independent audit of the evaluation scripts (`src/corio_contrast.py`, `src/armB_gtdb_eval.py`,
  `src/armB_plan.py`; `scratch/audit_eval/`, synthetic fixtures, no GTDB Arm B output read), during the GTDB Arm B
  run and before any of its output was evaluated. On the real label and plan files nothing changes a number: the
  family pick equals armB_select's, all 1,934 Coriobacteriales labels resolve to distinct GTDB accessions, the run
  plan covers the 4,599 selected genomes once, and an independent recomputation of the Arm A contrast matches. Fixed
  (pre-fix copies `scratch/*_pre_audit.py`; fixtures pass after, `scratch/audit_eval/fx_test_post_fix.out`):
  (1) `armB_gtdb_eval.py` reads the Coriobacteriales files first, so a Coriobacteriales genome rerun in the GTDB
  run with another table is scored on its new rows (none today); (2) Arm A flags whose codon has no Arm B row (TGA
  flags in suspect table-4 genomes run with table 11) are counted as untestable, not as unconfirmed; (3) a genome
  with several armc reasons counts once per reason class; (4) the evaluation refuses unless every selected genome
  has all its codon rows (61, or 62 under table 4) and was run with its selected table; (5) `corio_contrast.py`
  reads only `results/gtdb/armB_<i>.tsv` (it also picked up `armB_passes.tsv`) and (6) refuses the Arm B side
  while any family's genome lacks Arm B rows. Noted, unchanged: the pass rule is applied to the rounded output
  (ratio 3 decimals, log10p 2); passes within rounding of a bound are listed by the evaluation. `armB_plan.py` does
  not consult the stand-in map or look for GCF->GCA twins; no selected genome needs either.
- 2026-11-11 06:08 EST, check 3 (raw reads) made operational, before any read was mapped (`src/reads_check.py`). Sites: every Arm B
  event of the candidate genome whose window holds X, one site per event, at X's window offset nearest e*+1 (ties:
  the earlier); its genome position comes from the gene id and is verified against the assembly. Control sites: in
  the same events, the codon of query residue floor(ql/2) (inside the homolog-aligned part). Reads: the candidate's
  WGS runs from `data/reads/plan.tsv` (amplicon runs and runs with no fastq are skipped), streamed from ENA into
  minimap2 (`-x sr` Illumina, `-x map-ont` Nanopore; mates mapped as single reads; `--secondary=no`), never stored;
  only alignments overlapping a site are kept. Counted: primary alignments with MAPQ >= 20 that cover all three
  codon bases with no deletion or insertion inside the codon and base quality >= 20 (Illumina) or >= 10 (Nanopore)
  at all three; the codon is read on the gene's strand. A site supports X when depth >= 5 and >= 90% of counted
  codons are X; the candidate passes check 3 with >= 10 supporting sites (the PREREG rule). Reported with it: sites
  with depth >= 5 that do not support X, and the same numbers at the control sites. Streaming stops early once
  >= 95% of all sites have >= 30 overlapping primary MAPQ >= 20 alignments (fastq order is not position order, so
  an early stop is a random subsample); otherwise runs are read in plan order until the end.
- 2026-11-11 06:12 EST, check 3 amended before any read was mapped: ENA refuses fastq downloads from the lab's IP (HTTP 403), so
  reads are streamed from the NCBI SRA mirror with `fastq-dump --stdout --split-spot --skip-technical` in 4
  parallel spot ranges per run (as in `stop-codon-census/src/read_support.py`); runs and spot counts as in
  `data/reads/plan.tsv`. As in the TAG read support, alignments with NM > 3% of aligned bases are not counted
  (reads of relatives in the same metagenome), in addition to MAPQ >= 20; the stop rule counts the same alignments.
- 2026-11-11 06:12 EST, check 3, Nanopore runs (written before any was mapped): the NM limit there is 10% of aligned bases
  (Nanopore error rates are above 3%); everything else as for Illumina.
- 2026-11-11 06:17 EST, check 4 (independent re-check) made operational, before it was run. A fresh subagent gets a blinded list
  of 65 genomes (the 31 Arm B candidates mixed with 34 controls: non-candidate genomes of the candidates' genera and
  families, Polyangiaceae fetched from NCBI, two *Streptomyces*, *E. coli*, *B. subtilis*), under neutral ids, with
  no word of which are candidates or which codons are suspected. It picks a fixed set of conserved single-copy
  proteins per domain (queries from organisms outside the tested families; listed before any genome is scored),
  runs tblastn against each genome's contigs (no gene calls), and keeps per marker the best hit that covers >= 60%
  of the query and ends within 3 residues of its C-terminus. e* = subject codon aligned (by extrapolation) to the
  query's last residue; window = the subject codons at e*+1 .. e*+3 in the hit's frame. Per genome and sense codon
  X of table 11: O = hits with X in the window; f = X's frequency among the subject codons of the aligned parts of
  the same hits; E = max(n_hits x (1 - (1 - f)^3), 0.5); log10p = Poisson upper tail. A codon passes in a genome when
  O >= 5, O/E >= 3 and log10p <= -4; it reports every pass in every genome, plus the share of canonical stops
  (TAA/TAG/TGA) at e*+1. Check 4 confirms a candidate (genome, X) when X passes in that genome. Unblinding and the
  verdict are mine; the key is `results/check4_key.tsv`.
- 2026-11-11 06:21 EST, check 1 (clade) made operational, before any added genome was run. It is met for a candidate when another
  GTDB assembly of its genus already passes X in GTDB Arm B (Palsa-504, Trebonia, JAFARZ01, JBEKHI01, JADLFU01,
  JANWXX01, Methanoprimaticola, UBA168). For the other candidates, Arm B (frozen rule, table by GTDB's 11-or-4 rule)
  runs on every other GTDB R232 assembly of the candidate's species, plus the genus's not-yet-run GTDB assemblies up
  to 20 added genomes per genus, ranked by CheckM2 completeness - 5 x contamination (ties: accession). Check 1 is met
  if any added genome passes X; the count passing / run is reported. Candidates with no other assembly in their
  genus (DASTPP01, JAXLYQ01) stay unconfirmed by check 1.
- 2026-11-11 06:37 EST, check 3, after the independent audit of src/reads_check.py (scratch/audit_reads/). Two candidates had
  been mapped (DASTPP01 GB_GCA_035762175.1: 16 of 16 X sites supporting; JADLFU01 GB_GCA_963169165.1: 10 of 10),
  and two runs had been started and stopped. minimap2 -x sr had paired the interleaved mates, which lifts the MAPQ of a
  read in a repeat through its mate. Reads are now mapped one by one (--frag=no), which can only lower
  counted depth, and all candidates are mapped again; the paired-mode outputs are kept apart
  (data/reads/paired_v1/, results/reads/paired_v1/) and are not used. Also from the audit:
  - a reader-thread error, or a range for which fastq-dump reports fewer spots written than requested, now aborts
    the run;
  - one candidate codon per genome is asserted;
  - runs are de-duplicated.
  Thresholds unchanged.
- 2026-11-11 07:00 EST, check 3: spot ranges now take the spot count from SRA itself (vdb-dump --info, SEQ), not ENA's read_count.
  For ERR10643344 (JADLFU01 GB_GCA_963169165.1), ENA lists 130,832,460 reads but SRA has 65,416,230 spots. Ranges 0
  and 1 covered all 65,416,230 spots (fastq-dump reported each written), and ranges 2 and 3 were empty. The new
  written-spots check aborted the run, as it should. The candidate is mapped again. Rules and thresholds unchanged.
- 2026-11-11 07:09 EST, check 4 unblinded (`src/check4_unblind.py` -> `results/check4_verdict.txt`). Under the frozen rule, 12 of 31
  candidates are confirmed: 4 Polyangiaceae TTA (JADLFU01 x2, JANWXX01 x2), 7 archaeal CGA, *Achromobacter* GGA.
  0 of the 14 Streptosporangiaceae TTA are confirmed. Before unblinding, the blinded subagent reported that the
  rule's E (X's frequency in aligned coding sequence) is too small for windows after a stop, which are non-coding.
  17 of 34 controls carry some passing codon. CGA passes in three non-archaeal genomes where CGA is not the candidate
  codon, one of them a control. Under its downstream-codon diagnostic (E from in-frame codons e*+4..e*+30, computed
  blind), only the four Polyangiaceae TTA keep O/E >= 3.
  **Added for calling a result, stricter only:** a candidate must also keep O/E_down >= 3 in that diagnostic. This
  can only remove passes. The frozen-rule verdict is reported as is.
- 2026-11-11 07:12 EST, correction to the 07:09 entry: CGA passes in **four** non-archaeal genomes where CGA is not the candidate
  codon, not three: control Planotetraspora (g36), Palsa-504 g04 and g54, and JADLFU01 GB_GCA_963169165.1 (g34).
  Nothing else changes.
- 2026-11-11 07:29 EST, independent audit of the check-4 unblinding (`scratch/audit_check4/REPORT.md`).
  - **Checks out:** `src/check4_unblind.py` and the downstream diagnostic (`check4/src/qc_pass_positions.py`). All
    rows reproduce, and the four Polyangiaceae TTA keep O/E_down >= 3. g61 and g34 sit on the 0.5 E floor, where
    the stricter rule cannot remove a TTA. With the 95% upper bound on the downstream count, the four ratios are
    still 9.3, 8.7, 5.0 and 3.9.
  - **Withdrawn:** my post-unblinding comparison `src/check4_vs_armB.py` (07:1x) was not like-for-like.
    - The observed count included TTA after a canonical stop at e*+1, but the prediction did not.
    - Its denominator took every hit without a stop at e*+1, about 5x more than the hits that are Arm B events.
    - Its "Polyangiaceae 43 vs 43.7" and "Streptosporangiaceae 10 vs 48.4, P = 2e-11" are withdrawn.
  - **Corrected reading.** Check 4 had almost no power for the Streptosporangiaceae TTA.
    - Their expected O is about 0.2-0.6 per genome.
    - Arm B has TTA events at only 5 check-4 marker loci across the 14 genomes.
    - Their 0 of 14 is uninformative, not a refutation. They stay unconfirmed, because the frozen chain needs check 4.
  - **Polyangiaceae.** Arm B has TTA at 27 check-4 marker loci, and check 4 finds TTA at 18 of the same loci. Check 4
    also finds TTA at 23 marker loci that Arm B has no event at.
  - Rules and thresholds unchanged.
- 2026-11-11 07:35 EST, an outside re-check offered by Clover (post hoc; not one of the frozen checks, and it changes no verdict).
  - `src/blind_panel.py` writes 16 Polyangiaceae genomes to `blind_panel/`, shuffled with seed 20261111 and with
    contigs renamed.
  - The panel is the 8-tip TTA clade plus the 8 out-of-clade tips nearest the clade's MRCA by branch length.
  - The key is `results/blind_panel_key.tsv`.
  - Clover chooses and freezes the method. I will report its calls against the key as they come.
- 2026-11-11 08:08 EST, adversarial review of the Polyangiaceae TTA leads (post hoc; fresh subagent; `scratch/redteam/REPORT.md`).
  It changes no verdict and no rule. Not broken: on close-relative orthologs (no Swiss-Prot) TTA is the top codon
  after the sister's last residue in all 8 clade genomes (pooled 252 vs about 9.6 expected) and about 0 in the 8
  sisters. Weak points recorded: interior TTA in conserved genes at 1-1.8x the canonical-stop rate (ambiguous), and
  JAXLYQ01 shows the signal with an SPF RF2. The mechanism (release-factor stop vs unassigned codon) is open.
- 2026-11-11 08:26 EST, machinery on the blinded panel (post hoc, reported only, no verdict changes; `notes/candidates.md`).
  tRNA-Leu(UAA) absent in 8 of 8 clade genomes and present in 8 of 8 sisters. RF2 has more clade-derived fixed
  differences than any of 205 shared single-copy genes (10 in 262 sites; outgroup-polarised), including L201F next
  to the SPF motif in all 8, JAXLYQ01 included. JAXLYQ01's only RF2 copy on raw contigs is SPF.
- 2026-11-11 08:30 EST, Clover's blinded re-check unblinded (post hoc, outside the frozen chain, no verdict changes;
  `results/blind_panel_unblind.tsv`).
  - Clover's 6 UUA calls are all clade genomes. The other 2 clade genomes are Clover's only "strong" and "suggestive"
    non-calls, and all 8 controls show no signal (lowest-8 split by pB, chance 7.8e-5).
  - Correction to the 08:26 entry: "205 shared single-copy genes" should read "205 genes with reciprocal-best-hit
    orthologs in all 16 genomes". They were not filtered to single copy. RF2 ranks 1 by count and by rate.

### 2026-11-11 09:16 EST: fact-check of the write-up (post hoc; no verdict changes)
- An independent subagent fact-checked the outbox draft, the figure and RESULT.md (`scratch/factcheck/REPORT.md`).
  The literature framing was corrected (details in `notes/candidates.md`, last section).
- New post-hoc machinery numbers: RF2 SPY is in 9.7% of GTDB genomes, so SPY alone is not specific. F201 + SPY is in
  exactly the 7 SPY clade genomes (`scratch/rf_panel/spy_prevalence.py`). Reported only; no rule changes.

### 2026-11-11 09:43 EST: check 3 data incident, fixed (no rule changes)
- Lane B's queue repeated GB_GCA_020849655.1, which lane D was still mapping. At 09:39 a second reads_check started
  on it and opened lane D's in-progress SRR11294028 part SAMs for writing, truncating them (NUL bytes in 3 of 4
  parts). The final tally re-reads the SAMs from disk, so that run would have been corrupted.
- Fix: both processes killed by PID. SRR11294028 had no .done, so a restart re-streams it from scratch, and the 16
  finished runs (all .done before 09:37) are reused. `src/reads_check.py` now takes an flock per genome
  (`data/reads/<acc>/.lock`), so a second run on the same genome exits.
- Check 3 for JAXLYQ01 GB_GCA_046351315.1 (lane B, not affected): 11 of 11 testable sites support TTA -> PASS. It
  stays unconfirmed under the frozen chain (no other assembly in its genus; check 4 not met).

- 2026-11-11 10:03 EST, check 3 for JADLFU01 GB_GCA_020849655.1 (no rule changes): 18 of 18 testable sites support TTA -> PASS,
  after the 09:39 re-stream (clean SAMs). It meets every frozen check, the third genome to do so.

## 2026-11-11 10:57 EST: the non-representative GTDB members of the clade (no rule changes)
Frozen Arm B rule and checks 1–3 applied, unchanged, to the two non-representative GTDB R232 genomes of the clade's
genera (`src/clade_members.py`).
- GB_GCA_963169955.1 passes Arm B, check 1, check 2 (77.5×) and check 3 (20 of 20). Check 4 not run.
- GB_GCA_033578855.1 passes Arm B, check 1 and check 2 (17.4×); check 3 running.
- GB_GCA_025057345.1: check 3 PASS (16 of 16).
- CAYNEK01 GB_GCA_964396175.1, frozen Arm B: O_fam 5, so it does not pass.

## 2026-11-11 11:37 EST: timing notes and check 3 for GB_GCA_033578855.1 (no rule changes)
- The 10:57 entry was written after those genomes were scored (10:35–10:39). CAYNEK01's frozen Arm B was first run at
  07:16 (`results/spy_check/`); the 10:35 rerun gave the same values.
- GB_GCA_033578855.1, check 3: 15 of 15 testable sites support TTA -> PASS. Same read runs as GB_GCA_025057345.1.
- Reported for completeness: under the frozen check-4 rule, before the 07:09 addition, archaeal CGA in
  Methanoprimaticola GB_GCA_017540825.1 meets checks 1, 2 and 4, and check 3 (33 of 33 sites, 10:18). Only the 07:09
  O/E_down rule (0.77 < 3) keeps it from being called. The write-ups now say so (fact-check 2).

## 2026-11-11 12:17 EST: catalogue genomes outside GTDB R232 (no rule changes)
Frozen Arm B (src/term_codon.py, table g) and checks 1–2 as before, applied unchanged to the 9 non-redundant candidate
clade members from the marker-AAI search of NCBI and MAG catalogues (`results/catalogues/`, searched 10:45–12:10 by a
subagent; its member rule was fixed from controls before any candidate was scored). Arm B run 12:07–12:16
(`results/clade_members/armB_cat*`, `src/catalogue_members.py` -> `cat_members.tsv`, `cat_candidates.tsv`).
- Pass Arm B and check 2 (TTA): GCA_982529355.1 (1-FT3.2, co-culture, Finland; O_fam 102), bin496 (40), bin242 (29),
  MAG.1207 (28), spire_mag_00879468 (21), bin474 (16).
- Below the bar: GCA_057137765.1 (O_fam 9, TTA rank 1), spire_mag_00671464 (8, rank 1), GCA_035762045.1 (2, rank 4;
  56% complete, 78 events).
- Check 1 as frozen (a GTDB genus mate passes) does not apply cleanly to genomes without a GTDB genus; reported with
  the nearest genus by marker AAI.
- Check 3 plumbing (`src/reads_check.py`, backup `scratch/reads_check.py.bak-1212`): reads catalogue candidates from
  `cat_candidates.tsv` with their own events file and genome path; plan `data/reads/plan_cat.tsv`; Element AVITI
  runs are treated as short reads (minimap2 -x sr, the Illumina thresholds unchanged). Site rules and thresholds
  unchanged.

## 2026-11-11 12:34 EST: check 3 counts read pairs; catalogue metadata (no threshold changes)
- An independent audit (`scratch/audit_reads_cat/REPORT.md`) found that check 3 counted both mates of a pair when they
  overlapped a site (minimap2 `--frag=no` maps mates separately). Depth now counts one alignment per read name per site
  (`src/reads_check.py` tally, from 12:25; jobs started before then used the old count). All finished genomes were
  re-tallied from the kept SAMs (`src/reads_retally.py` -> `results/reads_pairs/`, `compare.tsv`): every check 3
  verdict is unchanged; two genomes lose one supporting site each. Thresholds unchanged.
- 1-FT3.2 (GCA_982529355.1) check 3: 104 of 104 (pairs: 104 of 104).
- The 1-FT3.2 completeness in `results/catalogues/members.tsv` was a marker proxy; replaced by the deposited CheckM2
  values (81.9 / 0.38). No verdict depends on it.
- Archaeal CGA: a close-relative ortholog test (red team's method, codon as parameter) was run after the verdicts,
  and is reported as post hoc; it does not change the 07:09 rule or the CGA verdict.

## 2026-11-11 12:57 EST: blinded CGA panel (archaea), rule fixed before launch
- Why: the CGA red team (`scratch/redteam_cga/REPORT.md`) found that the 07:09 O/E_down ≥ 3 rule fails TAA, a real
  stop, in 4 of 7 CGA-signal genomes when computed on the same close-relative events (downstream CGA is 2–3% of
  codons in non-coding DNA). So check 4 cannot decide CGA either way (the frozen rule is too loose, the added rule
  too strict for a codon common in non-coding DNA). A blinded close-relative re-check by a fresh agent is the test.
- Panel: the 16 genomes of `scratch/cga_ortho/panel.tsv` (9 Methanomethylophilaceae flagged for CGA, 7 sister-genus
  controls), renamed and shuffled; key in `scratch/blind_cga_key/key.tsv` (not shown to the agent). The agent is not
  told the codon, the lineage or which genomes are controls; it scores all 61 sense codons in every genome with its
  own method and freezes its own tiers before scoring.
- My prediction (from `scratch/cga_ortho/`): CGA in the 6 strong signal genomes (GCA_017540825.1, 017519045.1,
  015063135.1, 964445265.1, 016281925.1, 028561875.1); weak or nothing in 006954465.1; nothing in 000350305.1,
  006954405.1 or the 7 controls.
- Pass ("CGA confirmed blind"): the agent calls CGA (any positive tier) in ≥ 5 of the 6 strong signal genomes and in
  none of the 7 controls, and calls no other sense codon in any of the 16 genomes. Also reported: the rank of the 6
  strong genomes among all 16 by the agent's CGA score (top 6 exactly: P = 1/8008).
- If it passes, CGA is reported as a second case confirmed by a blinded re-check designed after check 4's rule
  failed, with the frozen chain's own status stated as is. If it fails, CGA stays a candidate.

## 2026-11-11 13:24 EST: blinded CGA panel unblinded
- The agent's reply was saved (`scratch/blind_cga/REPORT.md`, outputs `scratch/blind_cga/work/`) before I opened
  the key (`scratch/blind_cga_key/key.tsv`). Its rule was frozen at ~13:05 in its FROZEN.md, before scoring.
- Strong CGA calls (rank 1 of 61, P 1e-28 to 1e-72): G03 = GCA_964445265.1, G04 = GCA_016281925.1,
  G06 = GCA_017519045.1, G11 = GCA_015063135.1, G13 = GCA_028561875.1, G15 = GCA_017540825.1: exactly the 6 predicted.
  G10 = GCA_006954465.1 borderline, not called (predicted weak or nothing). G02 = GCA_006954405.1 and
  G09 = GCA_000350305.1 "untestable" (CGA nearly absent in genes). The 7 controls: no CGA. No other codon called at
  any of its call tiers (four "weak, not a call" hits; its ±6 shifted null gives 1 weak).
- By the rule above: **pass** (6 of 6 strong genomes, 0 of 7 controls, no other codon). The 6 are exactly the top 6 of
  16 by its CGA score (P = 1/8008 = 1.2e-4 under random order).
- Its own disclosed changes: a strand bug in its partial-gene flag, fixed mid-run (calls unchanged); a TAG-reading
  rule that flagged nothing, plus a table-15 sensitivity run (calls unchanged); its TAA positive control failed
  (power for a stop used in <10% of genes is limited). It asks for an independent audit of its scoring code.
