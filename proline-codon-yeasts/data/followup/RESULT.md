# RESULT: proline-codon yeasts (Saiph, 20 Nov 2026)
Rules: RULES.md (frozen 10:30, md5 08f7b4ef; A1 10:45 md5 004b72da; A2 md5 34269c4c). Pre-analysis choice: Eremothecium
sinecaudum left out of the Ashbya control sisters (its CUN meaning is not established).

## Arm T: which tRNA reads CUU/CUA as Pro
Controls PASS. PC-T1 Ashbya tRNA-Ala(UAG): D(UAG vs Ala) = 0.094, below all 20 Ashbya ortholog calibration values
(0.122-0.419, vs Ala); G:U3:70 unique vs the 4 sisters' A:U. PC-T2 S. cerevisiae tRNA-Thr(UAG): nearest class His 0.725
(fragile: Trp if gap columns are excluded; audit).

Wickerhamiella (W. siamensis BK074538, W. infanticola BK074532, W. verensis BK074542): one UAG tRNA per genome, 73 nt.
- D(UAG vs Pro) = 0.262 / 0.262 / 0.214 (calibration D_X min 0.136-0.157, median ~0.39): the CUN reader is the ortholog
  of the sisters' (BK074528/37/40) UAG tRNA, not a duplicated tRNA-Pro. Same locus (Thr-Ser-[UAG]-Arg-Trp) in targets
  and sisters (audit). Pro is the least similar class for 2 targets, 18th of 21 for the third.
- It gained C1:G72 (acceptor 5'CTAATAG/3'CTATTAG vs sisters' ATAATAG/CTATTAT; length 73 vs 78). C1:G72 is in tRNA-Pro of
  48 of 49 Dipodascales genomes (exception: W. siamensis's own tRNA-Pro, A1:U72) and in 0 of 1,083 tRNAs of other
  classes; the only other UAG tRNA with it is W. versatilis BK074543, which hardly uses CUN (all_aligned CTN 8; 0 CUN at strict Pro columns, 1 CTA at a strict Leu column;
  Clover's pro_leu_wss.tsv / cun_CTA.tsv). On the midpoint-rooted wt.treefile, W. versatilis falls outside the 15-genome
  Wickerhamiella clade that holds the targets and their sisters, so its C1:G72 is a separate gain: C1:G72 in a UAG tRNA
  arose twice, once with CUN = Pro and once in a genome with too few CUN codons to say what it reads.
  No decoy class gained its specific bases on this branch; but 10 of 21 classes changed acceptor bases on the branch.
- Two of the three sequences are identical and the third differs at 2 nt: one event in one clade, not three.
- N73 = A everywhere (no change at the discriminator).

S. pinicola BK074262: no Pro signal. D = 0.105 (0.123 without one bad sister call; 2nd-lowest of 23 calibration
values); G1:U72 and U73 (sisters A73), found in no sister; its one Pro-specific gain (p6 A) is matched by decoy classes.
So the two origins did not converge on the same tRNA change.

T0b: no CUN-reading tRNA missed by ARAGORN in the 4 targets (independent rerun agrees); a probable 2nd Pro(UGG) locus
in S. pinicola (273-345) reads CCN.
Note: the 1:72 tally uses ARAGORN's acceptor strands (tRNA-Pro: 53 C:G, 1 A:T over 54 loci in 49 genomes). An RNAfold
parse of the raw ARAGORN sequences disagrees (C:G 25, U:A 23) because 23 of the 54 tRNA-Pro calls (74-75 nt) start 1-2 nt
before the acceptor strand ARAGORN reports; folded with them, the extra 5' T pairs with A73. Audit a6 (T-loop-anchored) agrees with the ARAGORN tally.

## Arm R: RNA editing
R0: no RNA or proteome exists for any target species/strain (SRA, ENA incl. strain searches, PRIDE, ProteomeXchange).
STAT not scanned (~8.4 GB > 5 GB rule); Logan Search has no API. R1 therefore not run.

PC-R (A1): DRR917134 was organelle-poor (0.07% organelle reads) -> ERR15942850 (Col-0 background, crwn3-1), 4 M spots:
29,613 coding positions at depth >= 20 (3,375 U at codon pos 2); 452 coding RDDs, all C>U; 389 at codon pos 1-2. PASS.

R2 (species check >= 1% of reads on the own mitogenome at >= 95% id; A2 coverage floor):
| run | species / strain | library | mito share | coding pos tested | U>C | verdict |
|---|---|---|---|---|---|---|
| SRR5487602 | W. domercqiae Y-6692 | RT-PCR, PE | 0.0001 | 0 | 0 | not testable |
| SRR30194991 | W. domercqiae PYCC 3067 | polyA SE | 0.0041 (200k precheck) | - | - | not testable |
| SRR5942407 | W. versatilis Y-6652 | polyA PE | 0.0004 | 1,652 | 1 (cp3) | not testable |
| SRR5487603 | S. americana YB-2067 | RT-PCR, PE | 0.0001 | 0 | 0 | not testable |
| SRR30194995 | W. versatilis PYCC 7102 | polyA SE | 0.0381 | 4,456 (606 U at cp2) | 5 (cp1 1, cp2 1, cp3 3) | testable (A2) |
SRR30194995 detail: 42 coding RDDs, every one at alt fraction >= 0.8 (fixed differences, not partial editing); 30 of 42 at
codon pos 3; types U>A 10, G>A 8, A>U 7, U>C 5, C>U 5, A>G 3, G>U 2, A>C 1, U>G 1. DNA filter = WGS of a different strain
(Y-6652, SRR5938403), so these are most likely PYCC 7102 vs Y-6652 strain differences; the one cp2 U>C (cob 12142, AUG->ACG)
is T in 927 of 928 Y-6652 WGS reads. Under the frozen rule: 5 U>C < 10.
WGS baselines: SRR5938400 0.019, SRR5938403 0.134, SRR16989166 0.025 (all pass).

## Verdict
TBD after SRR30194995, Sabik.
