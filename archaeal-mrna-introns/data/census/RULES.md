# Frozen rules: archaeal protein-gene intron census (Mintaka)

Written 19 Nov 2026 before any census cut or control score was read. Changes after this point go in a dated
"Amendments" section at the bottom with the reason, and the original text stays.

## Data
- Genomes: every GTDB (R232) archaeal genome in the codon-atlas cache (3,388), run 1; plus the GTDB Asgard species
  representatives not yet cached, fetched into the atlas cache, run 2. `data/run1.tsv`, `data/run2.tsv`.
- Reader: reading-ledger v0.7 exactly (`work/reading-ledger/src` checked byte-equal to `snapshots/v0.7` before the
  run), `--batch 40 --threads 2`, reference = `data/anchors_v0` as the ledger uses it.
- BHB scorer, removal finder and helpers: the checker's `bhb.py`, `splice_dp.py`, `common.py` from
  `agents/vega/kills/asgard_cuts/src/`, copied unchanged into `src/checker/` (md5 in `notes/md5.txt`).

## R1. Cut (the unit)
A ledger join with 15 <= d <= 300, -5 <= dcol <= 10, min(bitsA, bitsB) >= 40 and (shift != 0 or gap_stops not
empty) -- inserts.py's rule with D_MIN 15 (to admit the 19-nt M. sedula cbf5 intron) and D_MAX 300. One row per
(acc, contig, strand, a_end, b_start). Excluded, and counted: holds another locus (inserts.py rule), N in the stretch,
jtrunc, removal closer than 100 nt to a contig end, family in the fused-reference list
{PGKT, RPOBC, ISPDF, HDRA, SECDF, MOACB, URE23, GSIA} (reading-ledger RESULT section 13).

## R2. Removal
`splice_dp.analyse(s, a_end, b_start, prot, a_se, b_ss)` with defaults; candidates = strict mode, L >= 9 (as the
checker's 02_spliced.py). Best = highest score; near set = candidates within 5 score units of best. A cut is
placed if a best removal exists and its score is above both the best frameshift model (L in -2..2, not 0) and the
best read-through (permissive, L = 0). Unplaced cuts are counted, not scored.

## R3. BHB label
`bhb_score(s, x, y)` on the strand sequence (x = first removed nt, y = first nt after). A cut is BHB if any removal in
its near set is strict or relaxed (`is_strict or is_relaxed`). Strict-only is reported alongside.

## R4. Microhomology k
k at the best removal = how far (x, y) slides left plus right with the spliced product unchanged
(s[x-1-t] == s[y-1-t] leftward, s[x+t] == s[y+t] rightward). Headline counts use k <= 3; k >= 4 reported apart.

## R5. Decoys (same shape)
For a placed cut with best removal length L, the whole near set is moved as a block by
delta in {+-(L+30), +-(L+90), +-(L+150), +-(L+210)} nt along the same strand sequence (8 decoys), each kept only if
all its windows lie >= 100 nt from contig ends and contain no N. Decoy rate of the cut = share of its kept decoys that
are BHB by R3. Expected BHB count of a set of cuts = sum of their decoy rates.

## R6. Excess
Per genome and per GTDB class: obs = BHB cuts (placed, kept, k <= 3); exp = sum of decoy rates of the same cuts;
also the number of distinct host families among BHB cuts (copies of one gene are one event). Class test: one-sided
Poisson P(X >= obs | exp), Bonferroni over classes with >= 5 placed cuts. Class confidence interval by bootstrap over
genomes (2,000 draws, seed 1).

## R7. Recurrent sites
Site = (family, reference column of the last upstream residue +-3, phase). A site is recurrent if BHB cuts at it come
from >= 3 genomes in >= 2 classes. Reported with decoy-rate and whether the stretches are alike (12-mer sharing).

## Controls and bars (fixed now)
- C1 scorer: tRNAscan-SE verified introns (atlas `trnascan` product) in the census genomes, near set = the
  slide-equivalent positions (R4). BHB share >= 0.60; same decoys (R5 offsets with L = intron length) <= 0.08.
  If C1 fails the census BHB labels are not read.
- C2 blind cbf5: in the RNA survey's cbf5 control genomes present in run 1 (S. solfataricus P2 GCF_000007005.1,
  A. pernix GCF_000011125.1, S. acidocaldarius GCF_000012285.1, M. sedula), a placed BHB cut in family TRUB in
  >= 3 of the 4 species. Fewer = the census misses known introns and its zeros mean nothing.
- C3 L02: a cut in RBLL of GCA_025839675.1 whose removal lies within 10 nt of 4487605, length 31 +-3. Reported
  whether placed and BHB; a miss is reported as a recall failure, not hidden.
- C4 recall: taken from the ledger benchmark (13/76 archaeal 25-40 nt plants read as cuts, v0.7); counts are
  reported as observed, with that recall stated, and not scaled up.

## Census verdict (DNA only)
- Asgard excess: pooled Asgard obs >= 20, obs/exp >= 3, and >= 2 Asgard classes each passing the R6 class test.
- Otherwise: "no Asgard excess above decoys at this recall" and the per-class table.
- Either way the census is a list of DNA candidates; "spliced" needs RNA (shortlist #3 RNA step, after Vega's L02
  verdict, Vega's junction rule, templates not fragments; bar >= 5 sites in >= 2 lineages removed in RNA, 0 joined
  in DNA). B35 and B-36 RNA are Vega's and are not touched.

## First kill test (K1)
Census BHB cuts outside TRUB in the RNA survey's genomes are looked up in the survey's kept read tables (no new
mapping). A cut the survey's RNA reads straight through, with no removal, is killed; the survey's whole-genome
verdict (no new protein-gene intron in 81 archaea) is the prior these cuts must face.

## Amendments
- 2026-11-19 15:28, before K1 was run. K1 numbers: a census cut, converted to forward-strand [start, end), matches a
  survey removals row on the same contig if the starts differ by <= k + 3 and the lengths by <= 3. Matched row with
  templates >= 2: "removed in RNA" (a lead for the RNA step). No row and >= 10 survey through-reads (cov.npz, ok runs,
  the smaller of the two boundaries): "read through", killed. Otherwise "not covered". Applies to every placed cut,
  BHB or not, so the BHB and non-BHB removed rates can be compared.
- 2026-11-19 15:29, after C3 and K1 on batches 0-1 and before any per-class or Asgard statistic was computed.
  S. acidocaldarius cbf5 (GCF_000012285.1) is a join with dnt 1, dcol -7 (d 22): its pieces overlap into the intron,
  so R1's dcol >= -5 drops it. Added a separately reported "wide" tier: same rule with -10 <= dcol <= -6.
  Headline counts, C2 and the census verdict stay on the R1 tier. The wide tier is reported beside it, with its
  own decoys. C2's S. solfataricus P2 (GCF_000007005.1) is not in the atlas cache; the cached GTDB representative
  GCF_900079115.1 stands in for it.
- 2026-11-19 15:34, before run2 (which holds these genomes) was started and before any of their reads were fetched.
  RNA step (Vega's L02 verdict at 15:24: L02 stands, step open). Units: P. syntrophicum MK-D1 (GCF_008000775.2,
  survey's mapped DRR199588, read through K1 as is) and two Hodarchaeales cultures (Nature 2026, doi
  10.1038/s41586-026-11170-7): HC1 Ca. Margulisarchaeum peptidophila GCA_055385335.1 (AP029002.1) with RNA
  DRR620298, and SC1 Ca. Flexarchaeum multiprotrusionis GCA_055394625.1 (AP029247.1) with RNA DRR620299. HC1/SC1
  are not GTDB representatives: they run in ledger run2 (data/run2_all.tsv) but stay out of the census verdict
  pool and every census statistic; they are reported only in the RNA step. B35/L. ossiferum RNA is Vega's and is
  excluded. Baldrarchaeota PRJNA779200 is METAGENOMIC/PCR, 31-86 Mb runs: not RNA-seq of a genome, excluded.
  Mapping: the survey pipeline as gated (gapcount.py v3 md5 9afc46903c44e535be1ab76d7cc5ecbb and annot.py md5
  9b7fa0ba149962d86111bada8f4bc5c8, run unchanged from their folder; survey_map.sh md5 7497fecd9ba1fdad3b2d1f52a95be45b
  copied to src/rna_map.sh with only the ref and out roots changed), 8,000,000 spots per run (the survey default),
  2 threads. Templates, not fragments: the survey's template column. A run counts only if .done exists, .n equals
  2 x spots read, mm2.log ends with "Real time", and median identity >= 0.995 (annot.py's own gate).
  Words: K1 (amendment 1) on the placed census cuts of these genomes, both tiers, reported apart.
  Lineage = genus (MK-D1, HC1, SC1 are three genera in two classes). Site bar: >= 5 distinct census cut sites
  "removed in RNA" in >= 2 genera. DNA half: no DNA reads are public for HC1/SC1 (SRA and ENA searched by biosample
  and bioproject at 15:33, 0 runs), so for them "0 joined in DNA" cannot be tested; a site counted there is
  "removed in RNA, DNA unchecked (closed assembly carries the stretch)" and the bar can be "met in RNA" at most,
  never a full pass. MK-D1: PRJDB9032 holds only the RNA run (ENA, 15:34); same DNA status.
  Run check (positive control on these runs): tRNA-class gaps with >= 2 templates per genome reported; if none,
  the run's "not covered" words are marked "no control" (libraries are often size-selected against tRNA).
- 2026-11-19 15:46, after an independent code audit (scratchpad audit/, findings in notes/AUDIT_1.md). Written knowing the
  partial batch 0-2 census (my trial summary at ~15:30, Asgard far above the bar there; the auditor also computed a
  per-phylum table it did not report). Every change below either makes code match the frozen wording or makes a
  pass harder; none can help the Asgard verdict.
  (a) R6 class test: decoy rates are estimates (<= 8 decoys per cut) and a class whose decoys are all negative gave
      exp 0 and p 0 automatically. The class test becomes one-sided Fisher exact: BHB share among the class's placed
      kept k<=3 cuts against the BHB share among those cuts' decoys. Bonferroni as before. exp used in obs/exp
      (class and pooled verdict) = n_cuts x (decoy BHB + 0.5) / (decoys + 1). Cuts with 0 kept decoys are left out
      of obs and exp and counted.
  (b) Bootstrap: classes in sorted order with one generator per class (seed 1 + class index), resampling genomes;
      the CI is for obs/exp with the same pseudo-count exp.
  (c) R7: "reference column" is implemented as stated. Each family's reference = the anchor carrying the most BHB
      cuts. A cut's anchor column maps to the reference by the best ungapped match (BLOSUM62) of its anchor's
      31-aa window centred on the column, with a score of at least 30; cuts whose window does not map are reported
      unmapped. Sites are clustered on reference column +-3 and phase. Every site is kept, without deduplication by
      counts. The 12-mer "alike" share is reported.
  (d) C1: R wording is "tRNAscan-SE verified introns". tRNAscan-SE calls noncanonical (NCI) introns with a BHB
      model, so scoring them is circular. C1 uses tRNAs whose introns are all canonical (note "CI" only). Shares for
      every subset are reported.
  (e) K1 as worded: a matched row with < 2 templates is "not covered", not "read through". Output gains tier,
      strand and a flag for L < 15 (the survey keeps only gaps >= 15 nt, so such cuts cannot read "removed").
  (f) HC1/SC1 (AP029002.1, AP029247.1) are excluded explicitly in summary.py. The 20-nt BHB flank padding (edge
      check 120 nt rather than 100) applies to cuts and decoys alike and stays. Recorded as a departure.
- 2026-11-19 16:00, found while checking which accession MK-D1 runs under; written knowing the partial census
  (25 of 85 batches). The atlas gtdb table lists some assemblies under both their GenBank (GCA) and RefSeq (GCF)
  accessions, so 16 assembly numbers appear twice in run1+run2. In 15 the version matches: the twins have equal
  locus counts and identical cuts in every pair already read (contig names aside). The 16th is T. litoralis
  GCF_000246985.2 / GCA_000246985.3, two versions of one assembly. R1's unit is a genome, so counting one twice is a
  defect. Rule: where an assembly number occurs more than once, the GCF accession is kept and the GCA is dropped from
  every count (R6, R7, C1, C2, K1). The list is data/twins_excluded.txt (16). Only one twin is Asgard: MK-D1
  GCA_008000775.2 is dropped and GCF_008000775.2 (run1, rep, also the survey's accession) is kept. Dropping it
  removes a second copy of an Asgard genome from the Asgard counts, so the change cannot help the verdict.
- 2026-11-19 16:02, reporting only, not part of the verdict; frozen before the full census exists. MAG quality
  could make cuts (chimeric or misassembled loci). summary.py will also report, for the main tier, Asgard pooled
  obs / exp (pseudo-count exp) and the same for all non-Asgard genomes in three strata: (i) CheckM2 completeness
  >= 90 and contamination <= 5; (ii) the rest; (iii) genomes with <= 10 contigs. Reading frozen now: if Asgard
  obs/exp in stratum (i) is under 3, or under a third of stratum (ii)'s, the write-up says the excess leans on
  lower-quality MAGs. Stratum (iii) is shown as it falls.
- 2026-11-19 16:42, wording only (Audit 2, notes/AUDIT_2.md); written after partial results were seen, so each reading is
  the one the code already used, and in both cases it is the reading with fewer passes. (a) The Bonferroni family is
  the classes with >= 5 kept, placed, k<=3 cuts with decoys, i.e. the cuts the class test is run on (partial
  cuts1: m 38; counting any k would give m 40 and 9 Asgard classes passing, not 8). (b) "Sorted order" for the
  bootstrap seeds is by (phylum, class). (c) Amendment 6's "stratum (ii)'s" is Asgard's stratum (ii) obs/exp;
  summary.py now prints that reading as a line. (d) C2 is counted on the main tier (amendment 2) with excl empty;
  summary.py prints the count and the >= 3 of 4 bar, and the wide tier beside it. (e) Guards added: each genome
  once, no empty class. None of this changes a verdict on the partial data.
