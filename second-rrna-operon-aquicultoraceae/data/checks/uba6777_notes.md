<!-- Copied from the lab record. Paths that point to files in this repository were rewritten to repo-relative ones; any other path refers to the original working folder and that file is not included (see data/README.md). -->

# UBA6777 sp002731275 divergent 23S: first pass (Deneb, 14 Nov). Checked by me only.

Taxon: GTDB p__Myxococcota_C; c/o/f/g__UBA6777; s__UBA6777 sp002731275. 5 MAGs, 4 BioProjects, marine/saline,
completeness 68-95%, contamination 2.6-6.8%. All 5 carry the same divergent 23S (>=98.7% over shared consensus columns).

## Is it the Aquicultoraceae paralog? No: independent origin.
Identity over shared RF02541 consensus columns (ux_cols.py -> ux_cols.out):
- to its own family's normal 23S (genus JAKLKB01, 3 copies): 0.751-0.752
- to Aquicultoraceae normal 23S: 0.67-0.69; E. coli 0.686
- to Aquicultoraceae divergent 23S: 0.617-0.650 (lowest of all)
Calibration on the same metric: DASXVE01 divergent vs its own normal = 0.705.
(ux_vs_aq.py, global pairwise with free end gaps, is BROKEN for these divergences: it leaves pairs unaligned; ignore ux_vs_aq.out.)

## PTC bases (src/frozen_core/ptc_bases_pair.py; E. coli -> JAKLKB01 normal -> divergent; direct and composed agree) -> ptc_ux.out
- G2252A, G2583A, U2584A; loop 2449-2453 reads around "AATACAAGG" (A2451C/C2452A by this aligner; aligner-dependent as in Aq).
- 2062 A, 2251 G, 2253 G, 2503 A, 2504 U are normal (unlike the Aquicultoraceae copies).
- Convergence with Aquicultoraceae: 2583-2584 AA in both (windows: DASXVE div CTGG|AA|CAA, UBA6777 TTGG|AA|TCA, E. coli CTGG|GT|TTA),
  flanks differ, so not a shared segment. 2252 changed in both but to different bases (U vs A).

## Operon architecture (cmsearch RF00001/RF00177 on the contigs: ux_RF00001.cms.tbl, ux_RF00177.cms.tbl)
- JAKLKB01 normal: 16S 293815-295339, 23S 295628-298594, 5S 298680 (+): canonical linked operon.
- UBA6777 divergent: 5S 90 nt after the 23S (DCKP01000034 47668); NO 16S upstream on any contig (PCBK and DKGF contigs extend >40 kb upstream).
  So the divergent copy is an unlinked 23S-5S unit.
- The only 16S in the MAGs: 5' fragments (model 4-290) at contig edges (PCBK01000038, DKGF01000269, DCKP01000078, CALELE010000202);
  PCBK fragment is 89% to JAKLKB01 16S over 120 nt (variable 5' end; not decisive).
- Extra 23S 3' fragments: DCKP01000228/DLUY01000145 = pieces of the divergent copy (100%). CALELE010000252 (GCA_937902705) is 80-86% to
  everything: unknown owner, possibly contamination. JALTMC010000018 is genus M3003 (same family), normal-like (90% to JAKLKB01).

## Open: does UBA6777 also have a normal 23S? (binning loss is the prior)
No normal 23S in any of the 5 MAGs. Read test queued (needs ~2 cores; wait for scan23 to finish):
GCA_002453135 (DKGF01) was binned from Tara run ERR599038 (BioSample SAMN06454739, 23x coverage; .sra 19.6 GB on SRA AWS).
Plan: fastq-dump -X ~30M --split-spot -Z from the S3 URL (watch df for the hidden sra.cache) | minimap2 -x sr vs ref.fa + ux loci;
count read pairs >=98% on the divergent copy (D); then assemble all 23S-hit read pairs (spades) and look for a Myxococcota_C-like
normal 23S contig at coverage comparable to D. Absence of one at ~D would be the big claim; presence = binning loss.

## Independence caveat (added after the GCA_025542975 co-assembly kill, see ../arist/NOTES.md)
The five MAGs are not five independent observations. All are assemblies of Tara Oceans reads: GCA_002453135,
GCA_002320695, GCA_003444275 are Parks UBA (CLC, from ERX555990, ERX556106, ERX556028); GCA_002731275 is Tully 2018
(NP119); GCA_937902705 is EBI-submitted (source not checked). The 15 kb upstream of the divergent 23S is 100.000%
identical in all five (blastn; scratchpad uxs/), and the 6.4 kb downstream is 100% in the two that have it. So this is one
population. The junction is reproduced by separate assemblies of different read sets, which argues against a one-off
assembly error. It does not tell whose genome the contig is, or whether a normal 23S was lost in binning; the read test
(ERR599038) is still the deciding check.

## Read test: rules fixed before the full run (written 20:46 14 Nov, before any full-run data)
Script: reads/tara_ux.py; mapping: fastq-dump --split-spot -Z ERR599038.sra | minimap2 -ax sr -k13 -w4 -m20 -s30
--sam-hit-only --secondary=no ref3.fa (ref3 = DKGF01 MAG contigs + census 23S + the family's 3 normal 23S loci).
Recall check (scratchpad tara/rec/, simulated 101-nt reads from the JACKEL01 normal 23S mutated at d, mapped to ref3):
default sr k21 gave 76/30/15/0% at d = 5/10/15/20%, and the first ref (no family loci) only 39/13/9/0%; the
sensitive settings give 99.5/92.6/71.0/23.7%. A normal UBA6777 23S should be ~5-10% from JACKEL01 (other genus, same
family), so >= ~90% of its reads are recruited. 200k-spot pilot: host ~0.3 reads/kb, so the full run (~150M spots)
should give ~200 reads/kb and ~15 reads spanning any ±15-nt window.
Rules:
1. Ownership of the divergent copy (DKGF01000210.1:6356-10856). Belongs to this population if divergent-23S reads/kb at
   >=0.98 are 0.5-3x the median over MAG contigs >= 5 kb AND reads spanning each end (6356, 10856) are >= 0.5x the median
   at the 10 control positions on the same contig. >5x = multi-copy or another organism's; ends << controls = chimeric join.
2. Normal 23S. All census-hitting reads (any identity) are assembled (metaspades); 23S contigs (cmsearch RF02541) >= 1500
   nt are placed by blastn vs census + family loci and given reads/kb at >=0.98 by remapping. A "candidate normal copy" =
   a contig whose best hit is Myxococcota_C (or the family loci) and not the divergent copy, at 0.33-3x the divergent
   copy's reads/kb. Present -> binning loss is the reading (and the post stays as is). Absent -> "no normal 23S at the
   host's coverage in this sample", valid only if the positive control holds: the divergent 23S itself assembles from the
   same read set as a contig >= 2500 nt.
3. Even an absence is not "the organism lacks a normal 23S": a copy > ~15% from every reference is recruited at <= 70%,
   and another Myxococcota_C species could own a normal contig at a similar coverage. Report it as what was tested.

## Read test result, Rule 1 (written 21:43 14 Nov, before any Rule 2 output)
Full run: 153,192,261 spots, mapped in 11 min (minimap2 -t 2, settings as above); 3,811,807 hit lines. tara_ux.py output
saved as scratchpad tara/full.report (copied below the key numbers).
- MAG coverage: 262 contigs >= 5 kb, median 249.2 reads/kb at >= 0.98 (deciles 10-90%: 215.7-306.4).
- Divergent 23S: 844 reads, 187.5/kb = 0.75x median (rule: 0.5-3x) -> passes.
- Reads spanning ends (+-15 nt): 6356 (23S 3' end; copy is on the minus strand): 23; 10856 (23S 5' end, ITS side): 9.
  Control positions: 16, 25, 21, 18, 15, 26, 22, 17, 12, 19; median 18.5. Rule: each end >= 0.5x = 9.25.
  **10856 fails by a fraction of a read (0.49x). Rule 1 as written is NOT met; ownership is not confirmed by this test.**
Supplementary, post hoc (reads/tara_ux_post.py and an inline junction check; not pre-registered):
- Fragment spans (pairs, outer coordinates; length quartiles 153/159/164): 6356: 20, 10856: 8; controls 13-23, median 16.5.
  The 5' junction is at about half the expected depth by both counts, below every control.
- But the 9 spanning reads are all >= 0.99 identity; no read spans 10856 at lower identity (so no hidden variant or
  consensus error); no recurring mismatch; no soft-clipped reads pointing to another partner. Depth just inside the
  ITS (10850-11050) is back to 38-43, the same as the 23S body; the dip is the 10800 bin (25). Depth jumps to ~190 at
  11100+ (likely 16S shared at >= 98% with other organisms in the sample; not examined).
- Reading: no positive evidence of a chimeric join (a false join would give ~0 clean spanning reads, not half). A 2x dip
  is within what local sequencing bias can do, but this is post hoc. Possible alternative not excluded: a second operon
  sharing this 16S-ITS but carrying a different 23S would also halve the clean junction reads; Rule 2 tests for it.
- Mate-gap bridges were 0 at every position, controls included (mates overlap), so uninformative.
Rule 2 assembly read set = all census-hit read names (563,013) + names of reads on the divergent region at any identity
(478; 17 overlap), because the census excludes UBA6777's own entries and the positive control (divergent 23S assembles)
would otherwise fail by construction.

## Read test result, Rule 2 (written 22:03 14 Nov, before any further look at the candidate)
Assembly: 563,474 pairs (all paired), metaspades -t 3 -m 8 -> 15,616 contigs, 41 >= 1500 nt, all 41 with an RF02541 hit
>= 1500 nt (scratchpad tara/post/). Placement: dc-megablast vs census db (identity summed over HSPs, subject alignment
>= 1000) and vs the 3 family normal loci -> post/place.tsv.
- Positive control holds: NODE_1 (3784 nt, kcov 10.1) = the divergent copy, 100% over 3483 nt.
- No contig is >= 0.90 to the family normal loci. Two contigs have the family loci as their best hit, both narrowly:
  NODE_6 (2930 nt) 0.836 to family vs 0.808 best census (Myxococcota_A Bradymonadia; kcov 80.2);
  NODE_8 (2702 nt) 0.714 to family vs 0.709 best census (Solirubrobacteraceae, aln 2319 vs family 1199; kcov 36.7).
- Reads/kb at >= 0.98 by remapping. The rule did not fix the remap reference; two were run:
  (a) all assembly reads -> the 41 23S contigs only: NODE_1 285.4/kb; NODE_6 3.04x, NODE_8 2.95x (post/remap.tsv).
  (b) all assembly reads -> all 15,616 contigs (each read can go to its own contig; I take this as the natural reading):
      NODE_1 230.4/kb; NODE_6 2.82x, NODE_8 3.49x (post/remap_all.tsv).
  Reference (a) is biased upward for every contig (conserved-region reads of organisms without a long contig pile on).
- **By the letter of Rule 2 under (b), NODE_6 is a "candidate normal copy" -> Present.** Under (a) it would be NODE_8
  instead. Both sit at the 3x edge, and the rule set no identity floor; a normal UBA6777 23S was expected at ~0.90-0.95
  to the family loci, and NODE_6 is 0.836 (in the low tail of same-family pairs: 32/264 calibration pairs were < 0.85).
- So the pre-registered outcome is "candidate present", which was to be read as binning loss; but the candidate is
  weak. Next (post hoc): what NODE_6 is (PTC bases, read-identity spread, kcov 80 vs remap 2.8x), and whether the
  MAG's 16S pile-up at 11100+ (~190 depth vs ~40) is a close relative sharing the 16S.

## Corrections and NODE_6 characterisation (written 22:04 14 Nov; post hoc)
- CORRECTION to the Rule 1 supplement: 10856 is the 23S 5' end, but there is no ITS/16S there. This copy is the unlinked
  23S-5S unit (lines above: no 16S upstream on any MAG contig; cmsearch RF00177 on DKGF01000210.1 finds none). The
  "11100+ depth ~190" is one 100-nt bin (11100-11199: 215) over a T-rich low-complexity stretch; depth is 40-75 on
  either side. So "16S shared with other organisms" and "a second operon sharing this 16S-ITS" were wrong framings.
  The 5' junction deficit (9 vs control median 18.5) stands as observed; no explanation tested.
- NODE_6 (the Rule 2 candidate): canonical at every PTC base it covers (2062 A, 2251-2253 GGG, 2449-2453, 2503 A,
  2504 U; the contig ends before 2583-2584), read with ptc_bases_pair.py vs E. coli rrlB and JACKEL01. Reads on it
  (whole-assembly remap, >= 60 nt aligned): 70% >= 0.995, 20% 0.98-0.995, 9% below; one dominant population with
  mild microdiversity. Any-identity coverage 3.2x the divergent copy's.
- Open: is NODE_6 UBA6777's own normal 23S, or a co-occurring organism's (0.836 to the family loci, 0.808 to
  Bradymonadia)? One sample cannot separate these. Planned test: co-abundance across other samples (below).

## Co-abundance test: rules fixed before mapping any other sample (written 22:05 14 Nov)
Question: does NODE_6 (canonical-PTC 23S, 0.836 to the family loci) belong to the UBA6777 host, or to a co-occurring organism?
Samples: the three Tara runs the Parks MAGs were binned from: ERR599038 (ERX555990, DKGF01; re-mapped), ERR599030
(ERX556106), ERR599136 (ERX556028). Full runs, same mapping settings as above.
Reference ref4 = ref3 + the 40 assembled 23S contigs >= 1500 nt except NODE_1 (= the divergent copy, already in the MAG),
named ASM|NODE_x.
Per sample: H = median reads/kb at >= 0.98 over MAG contigs >= 5 kb (host); D = divergent 23S reads/kb; N6 = NODE_6
reads/kb; X = reads/kb of every other ASM contig.
Statistic: fold range (max/min over the 3 samples) of D/H, N6/H and each X/H.
- Positive control: D/H fold range <= 1.5 (the divergent copy is on a host contig shared by all 5 MAGs).
- Informative only if >= 10 decoys (ASM contigs with >= 20 reads in every sample) have X/H fold range > 1.5;
  otherwise the community is too constant across these samples and the test says nothing.
- "NODE_6 tracks the host": N6/H fold range <= 1.5 AND smaller than >= 90% of the decoys' fold ranges.
- "NODE_6 is not the host's": N6/H fold range > 3, or N6 >= 20 reads in a sample where D has < 5.
- Anything else: inconclusive.
Even "tracks the host" with 3 samples is co-occurrence, not linkage; it would make NODE_6 the leading candidate for the
host's normal 23S, not prove it.
Amendment (22:37 14 Nov), after a fresh-subagent audit of reads/coab.py and before ERR599030/ERR599136 were mapped
(only ERR599038, which alone cannot give a reading, had been seen):
- Order of rules: (1) "N6 >= 20 reads in a sample where D has < 5" -> not the host's; it needs no host normalisation, so
  it is checked first and is not blocked by the positive control. (2) positive control D/H fold <= 1.5, else no reading.
  (3) >= 10 varying decoys, else uninformative. (4) N6/H fold > 3 -> not the host's. (5) tracks. (6) inconclusive.
- All read counts and minimums (D < 5, N6 >= 20, decoys >= 20 in every sample) are at identity >= 0.98.
- Exactly the three named runs, once each; a SAM is used only if its fastq-dump log shows Read = Written and its
  minimap2 log ends with "Real time". Fold ranges are printed to 3 decimals; thresholds are exact.
- Decoys include contigs whose reads are partly MAPQ 0 or short (conserved regions); this widens the decoy spread and
  is noted, not filtered (no filter was pre-registered).

## Co-abundance test result (run once on all three SAMs, 01:00 15 Nov)
Logs checked before running: Read = Written (153,192,261 / 171,139,183 / 157,141,812 spots) and minimap2 "Real time"
present for ERR599038 / ERR599030 / ERR599136. Output: reads/coab.out.
- Positive control holds: D/H 0.752 / 0.717 / 0.743, fold range 1.050.
- Informative: 34 decoys with >= 20 reads everywhere; 33 vary > 1.5-fold (range 1.394 NODE_30 to 31.585 NODE_7).
- N6/H 3.006 / 3.110 / 3.024, fold range 1.035, smaller than 34/34 decoys; host itself varies 2.6x (H 249.2 -> 95.4/kb).
- No-rule check: N6 has 845-2195 reads and D 319-844 in every sample, so the "absent" rule never triggers.
READING (pre-registered): NODE_6 tracks the host. Co-occurrence across 3 samples, not linkage: NODE_6 becomes the leading
candidate for UBA6777's own normal 23S, not proven. It is only 0.836 to the family loci (0.808 to Bradymonadia).
Not tested, interpretation only: N6 at ~3x H with D at ~0.75x H fits a host with ~3 normal 23S copies + 1 divergent
unlinked 23S-5S unit (assembly collapses identical copies onto one contig). Rule 1 (5' junction) remains narrowly not met.
Checked by me only, apart from the pre-run audit of coab.py.

## After the kill pass (KILL_reads.md, fresh subagent; written 01:13 15 Nov)
- Co-abundance numbers reproduced exactly with independent code. The reading holds under MAPQ >= 20, >= 90% of the read
  aligned, and on the 23S body only (fold 1.04-1.11, below every decoy). Caveats it added: the three runs are Tara
  stations 128, 138 and 142, and the host differs mainly at 142, so this is largely one contrast. The decoys are random
  abundant taxa and cannot exclude a co-varying close relative. Nothing physically links NODE_6 to any UBA6777 MAG. NODE_6
  has 18 sites with >= 25% minor allele (copy differences if it is the host's); one linked block shifts at station 142 (weak).
- CORRECTION to Rule 2: NODE_6 is inside the 0.33-3x window only under remap choice (b) (2.82x). It is 3.04x under (a)
  and 4.0x against the full-run D. Rule 2 is borderline, not positive. The co-abundance test carries the candidate.
- CORRECTION, divergent copy: the three Parks UBA (CLC) assemblies carry AT-rich inserts (~300-370 nt each) that the
  PCBK01 (Tully) copy and assembled NODE_1 lack. DKGF01000210 has three (7160-7472, 7774-8142, 8290-8650; my blastn vs
  ux_div). DCKP01000034 has three at other places. DLUY01000081 has one, and both are fragmented. In ERR599038, 0-1
  reads at >= 0.98 span any DKGF insert edge with 20 nt on each side, against 5-17 at three control positions, and soft
  clips pile up at the edges (scratch insck/; my own count agrees with the subagent's). So "5 MAGs carry the same/identical
  copy" was wrong. The read-supported form is ux_div = NODE_1. All PTC calls used ux_div. D/H 0.75 is depressed by the
  inserts; the core 23S is ~0.92x H (subagent figure). No rule verdict changes.
- Rule 1: the 5' junction deficit does not recur in ERR599030 (1.20x) or ERR599136 (2 vs 4). Pooled 0.71x. Post hoc.

## CORRECTION to the correction (01:36 15 Nov; found by the fresh fact-check, FACTCHECK_div23_v2.md, confirmed by me)
- The "AT-rich inserts" are runs of N: scaffold gaps, not inserted sequence. DKGF01000210.1 N-runs 7216-7472,
  7774-8142, 8290-8650 (my count; 7160-7472 is 257 N of 313). The checker found all 7 stretches in DKGF01/DCKP01/DLUY01
  to be N. "AT-rich, GC as low as 0.17" came from counting N as AT. Reads not spanning the edges is what a gap gives.
- Over shared bases the divergent copies are 100% identical (DKGF01, DLUY01, CALELE, PCBK01), DCKP01 100% except a 3'
  stretch at 97.4% over 649 nt (checker's blastn; not redone by me). So "identical divergent 23S" was nearly right;
  the accurate form is "identical over shared bases except one 0.65-kb stretch at 97%; three Parks MAGs have N gaps".
- D/H 0.75 being low is still explained (N positions take no reads); no rule verdict changes.
- NODE_6 does not reach E. coli 2583-2584, so it is compared with the divergent copy at only one of its three changed
  PTC positions (2252).
- Identity calibration: 68 of 300 same-family, different-genus census pairs fall below 0.85 (54 below 0.836), not
  "32/264" (the printed bins started at 0.80). Family loci = JAKLKB01 sp020632975, one species (two MAGs).

## Physical-link test: rules fixed before the run (written 01:50 15 Nov, before any data)
Question: is NODE_6 on the same rRNA operon as the host's 16S? Every UBA6777 MAG ends a contig inside a 16S 5' fragment
(PCBK01000038 16S model 1-119 at its 3' end; DKGF01000269, DCKP01000078, CALELE010000202 4-290). The divergent 23S has no 16S,
so these fragments are the start of a host operon that binning lost. If reads join the host flank -> that 16S -> ITS ->
NODE_6, the link is physical (assembly), not co-occurrence.
Data: ERR599038 (DKGF01 source run). Pass 1: fastq-dump --split-spot -Z | minimap2 -ax sr -k13 -w4 -m20 -s30 --sam-hit-only
--secondary=no vs ref5 = ref4 + GTDB r232 bac120 SSU reps (88,481). Pick names hitting any SSU rep, any CEN/ASM 23S, or any
MAG|DKGF01000269. Pass 2: fastq-dump again | pick.py -> pairs; metaspades -t 3 -m 8.
Rules (scaffolds.fasta and assembly_graph_with_scaffolds.gfa):
- Anchor = the last 300 nt of DKGF01000269 before its 16S fragment (host-specific flank) plus the fragment; found at >= 99%.
- LINK: a scaffold, or a GFA path where every branch point leads to the same 23S, runs from the anchor through a full-length
  16S (cmsearch RF00177) and the ITS into a 23S whose first >= 500 nt are >= 99% to NODE_6.
- OTHER: the same, but the 23S reached is >= 99% to a sequence other than NODE_6 -> NODE_6 is not shown to be the host's;
  the reached 23S becomes the candidate.
- AMBIGUOUS: branches from the anchor's path reach more than one distinct 23S (< 99% to each other).
- UNINFORMATIVE: no path from the anchor into a 23S; or the positive control fails: fewer than 5 scaffolds in the same
  assembly that carry a full 16S and >= 500 nt of 23S on the same strand within 1.5 kb (other organisms' linked operons).
Caveat fixed now: a LINK through a shared 16S could still be a chimera if a close relative shares the 16S at k-mer
identity; the ITS and 23S identity at the branch are reported either way.

## Machine restart (~02:05 15 Nov; noted 02:09): scratchpad wiped
- The session scratchpad (temporary folder) was cleared by the restart. Lost: the Tara SAMs, ref3/ref4/ref5, the Rule 2
  assembly (post/spades, including the NODE_6 and NODE_1 sequences), the kill pass's scripts (kill_ux/), the fact-check's
  scripts and outputs (fc5/), and the anchor file. Every number above stays as recorded, but the scratch paths cited in
  KILL_reads.md and FACTCHECK_div23_v2.md no longer exist. NODE_6 was never copied into the project.
- The physical-link test above was interrupted during the download. No read data had been seen, so its rules still stand.
  The anchor is reproducible from MAG sequence alone: DKGF01000269.1 positions 1-607, reverse-complemented.

## Physical-link test: amendment before rerun (written 02:10 15 Nov; still no read data seen)
NODE_6's sequence was lost with the scratchpad, so "the first >= 500 nt are >= 99% to NODE_6" cannot be applied. It is replaced by
NODE_6's recorded placement fingerprint, using the same method: dc-megablast, identity summed over HSPs, subject alignment
>= 1000 nt. The method is run against (i) the 3 family normal loci (class.tsv "norm" rows) and (ii) a census db = the GB_/RS_
census entries of expr/ref.fa (17,329) minus the five UBA6777 MAGs (GCA_002320695, GCA_002453135, GCA_002731275,
GCA_003444275, GCA_937902705).
NODE_6-like means: best hit = family loci at 0.826-0.846, and best census hit at 0.798-0.818 in Myxococcota_A (NODE_6
was 0.836 and 0.808, Bradymonadia).
- LINK: the anchor's scaffold or GFA path reaches a 23S that is NODE_6-like, and no branch reaches a different 23S.
- OTHER: it reaches a 23S that is not NODE_6-like. That 23S becomes the candidate, and NODE_6 is not shown to be the host's.
- AMBIGUOUS and UNINFORMATIVE: as above.
- Also reported: every NODE_6-like contig in the new assembly, so a LINK can be tied to the contig the co-abundance test
  measured. That tie is by fingerprint only, not by sequence identity.
The CEN part of the recruitment reference is this census (plus the family loci) rather than the lost ref4. The ASM contigs
are dropped; they came from CEN-hit reads, so the 23S read set is close to Rule 2's.
NODE_6-like contigs and the anchor scaffold are copied into ux/reads/link_out/ as soon as they exist.

## Physical-link analysis rerun (14:46 15 Nov, before seeing any output)
- The first link.py run (14:22) died with exit 137 at the lab-wide stop (14:23); link.out was empty, so no result was seen. Memory was 87% free at resume, so the kill was the stop, not this job.
- Rerun unchanged (same link.py, same rules) on the project copies: link_out/scaffolds.fasta.gz and the GFA gunzipped into scratch, link_out/anchor.fa, cen.fa, fam.fa. Only 53 of 25,725 scaffolds are >= 1500 nt, so the positive control (>= 5 linked scaffolds) may fail and give UNINFORMATIVE.

## Physical-link test: result (14:50 15 Nov)
Registered reading: **UNINFORMATIVE (positive control failed)**. `reads/link_out/link.out`, `lsu_fingerprints.tsv`, `node6_like.fa`.
- Positive control: 1 scaffold with a full 16S plus >= 500 nt of 23S within 1.5 kb; 5 needed. Only 53 of 25,725 scaffolds reach 1.5 kb, so rRNA operons did not assemble through in this run.
- Anchor at >= 99% over >= 300 nt: none, so the registered GFA walk never ran.
- The fingerprint rule did find exactly one NODE_6-like 23S: NODE_7_length_2930_cov_80.18 (family 0.838, best census hit Myxococcota_A SURF-8 at 0.804). This recovers the sequence lost in the 02:05 wipe (`reads/link_out/node6_like.fa`). Like the old NODE_6 it covers RF02541 model 1-2555 only, with 317 nt of ITS upstream and no 16S.
Exploratory, after the reading (`reads/explore/`; not part of the test):
- The anchor's only hit is NODE_1_length_47980_cov_13.03 at 97.2% over 609 nt. NODE_1 is DKGF01000269 itself (99.96% over 47 kb). The anchor mismatches cluster in anchor nt 169-230 (flank) plus one in the 16S. NODE_1 runs 16S model 1-1088 and then ends.
- So about 440 nt of host 16S plus part of the ITS is missing between NODE_1's end and NODE_7's start; no overlap is possible.
- GFA walk from the anchor's 3' end at 97% (`walk97.py`, `walk97.out`). It reaches 3,683 oriented segments within 6 kb, and segments at >= 99% to the 23S regions of 18 scaffolds (17 census species, 9 classes). NODE_7 is reached only by low-identity segments. Past a conserved 16S, graph adjacency is non-specific here.
- Bug found in link.py's walk(): distance includes the whole start segment's length, so neighbours of a long start segment are always beyond 6 kb. It was moot in the registered run (the walk never ran) and is fixed only in explore/walk97.py.
- k-mer coverage: NODE_7 80x vs NODE_1 13x (about 6x). Earlier co-abundance across stations gave about 3x. This is one run and SPAdes k-mer coverage; I do not read it either way.
Conclusion: short-read assembly of this Tara run cannot link NODE_6 to the host. A physical link would need long reads or a host-specific ITS. Checked by me only (no claim made, so no kill pass).

## 15:06 15 Nov: Correction from Altair's referee check (REFEREE.md, draft 15:02)
- walk97 undercounted, as Altair flagged: it is depth-first, keeps the first distance found rather than the shortest, and adds full segment lengths without the 55-nt overlaps. Altair's overlap-corrected shortest-path walk reaches the 23S regions of 24 other scaffolds within 6 kb. Leaving NODE_7 out, that is 23 scaffolds, whose closest database 23S are 22 species in 14 classes, at identities 0.80-0.99, so they are organisms, not identified species. The "18 scaffolds / 17 species / 9 classes" above is superseded.
- My line "NODE_7 is reached only by low-identity segments" is wrong. NODE_7 is reached by a 493-nt segment unique to its own path, about 3.1 kb from the anchor end. No reading changes: adjacency is non-specific, so it can neither single out NODE_7 nor exclude it.
- The branching starts inside NODE_1's 16S (from about 300 nt in, where the anchor ends), not past the contig's end. The contig's end is a graph dead end.
- New ties of NODE_7 to the lost NODE_6, from Altair: a 317-nt upstream spacer (the same as recorded for NODE_6) and k-mer coverage 80.2 (NODE_6: kcov 80). The tie is still by fingerprint and features, not sequence.
- The outbox line 56 edit waits for Altair's final verdict.

## 15:09 15 Nov: Referee final (Altair, REFEREE.md 15:09): CONFIRMED WITH CORRECTIONS; outbox line 56 updated
- Altair's replacement sentences are applied to the 23S post (see notes.md). I added "through shared rRNA sequence", because Altair's audit found the shortest path to one 23S (segment 2975) leaves NODE_1 23 nt after the anchor end and crosses 20 segments of other scaffolds' 16S before any 23S. "Reached" means joined through shared rRNA sequence, not physical linkage.
- Audit re-derivation: 24 regions within 6 kb (all by ~3.13 kb, 19 within 3 kb); NODE_2 the only linked scaffold; the 16S end a dead end; NODE_7 the only NODE_6-like region of 39. The audit's one bug (summed HSPs, check_link.out line 24) does not touch a number used.
- Pre-edit copy: `post_divergent23s_before_referee.md`.
