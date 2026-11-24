# Independent check: is the L02 31-nt removal a library artefact? (C6)

Checker: fresh subagent, 19 Nov 2026. Brief: try to kill the L02 cut (CP104013.1:4487605+31, 0-based, slid left;
"AA" direct repeat) as a library-prep, mapping, or strain artefact. I read `confirm/RULES.md` only from `confirm/`, and
no commons/outbox posts. All my files and scripts: `confirm/indep/` (scripts in `confirm/indep/scripts/`). Every number
below was printed by a script. The log file for each test is named next to it.

Inputs: the three minimap2 SAMs (`culture/out/sam/*.sam.zst`, `-ax splice:sr`, mates mapped as single reads), which I
sorted to `indep/*.bam`. SRR22273278 = DNA (WGS, METAGENOMIC, library D81B_0405, NovaSeq 6000, the 2022 genome
submission PRJNA847409). SRR27313660 = RNA "A202_A" (2023, NovaSeq 6000). SRR32148620 = RNA "469B" (2024, NovaSeq X).
Both RNA runs: METATRANSCRIPTOMIC, selection "PCR", mirVana + TURBO DNase. The SRA records give no library kit.
NCBI feature table for CP104013.1: UYP47557 = CDS 4487084-4487611 (+), UYP47558 = CDS 4487703-4488440 (+). The stretch
removes UYP47557's stop codon.

## T0. Re-count from the SAMs (check the lab's tables)
`scripts/classify.py` -> `indep/classified.tsv`. Primary reads only. "removal" = a 31-nt N/D starting at 4487605-07;
"through" = one ungapped block from 4487595 to 4487648 or further.

| run | removal | through | partial / other gap |
|---|---|---|---|
| SRR22273278 (DNA) | 0 | 38 | 14 / 0 |
| SRR27313660 (RNA) | 253 | 20 | 87 / 3 |
| SRR32148620 (RNA) | 306 | 2 | 85 / 1 |

These match the lab tables (`rna_mol_sites.tsv`: 497 joined and 32 through fragments over the two RNA runs; DNA 0 and
46). Strandedness: mate 1 maps forward/reverse over the gene body 870/741 (SRR27313660) and 1048/856
(SRR32148620), so the libraries are **unstranded**. A sense-versus-antisense test is therefore not possible.

## T1. Mapping/aligner artefact
(a) Aligner-free k-mers (`scripts/kmerclass.py`). Each read sequence is checked for exact 16-mers across the joined
junction (6 or more nt each side beyond the AA repeat) or across either boundary of the stretch. The two sets share no
k-mers.
- DNA: 0 joined, 51 through (38 aligned through + 13 partial).
- SRR27313660: 285 joined (250 aligned removal + 35 clipped "partial"), 31 through. Removal share 0.902.
- SRR32148620: 342 joined (304 + 38), 7 through. Removal share 0.980.

(b) Read quality of the removal reads (`scripts/anchors.py`). "Clean" = total soft clip 5 or less, both anchors 20 nt
or more, and MAPQ 30 or more. L02 removal reads: 215 clean of 250 (SRR27313660) and 286 of 304 (SRR32148620). All have
MAPQ 30 or more. Most carry 0 mismatches outside the gap (NM excluding the gap = 0 in 176/250 and 268/304).
Clean-only shares: 215/(215+20) and 286/(286+2).
**Positive control:** the genome-wide gap scan found a 30-nt "removal" at 1291865 in the DNA run (8 removal, 13
through). All 8 DNA reads and all 5 RNA reads there fail the criteria: MAPQ 23, 6-91 nt soft clips, anchors 28+19.
None has a properly placed mate. They are foreign (metagenome) reads stitched by splice mode. The method flags this
artefact; L02 is not flagged.
**Verdict: no effect.** The removal is in the read sequences themselves, not made by the aligner.

## T2. PCR-stage skip: duplicate-family concordance
`scripts/families.py`, `scripts/siteanalyze.py`. Mates are paired by spot name. A family = fragments with identical
outer coordinates (soft clips included). Expected mixed families under an independent per-copy skip =
sum(1 - p^n - (1-p)^n).
- SRR27313660: p = 0.916 (fragments). 48 families with 2 or more fragments. **Expected mixed 8.29, observed 1**,
  Poisson P(X <= 1) = 0.0023. Across both runs, removal fragments' siblings in the same family: 316 removal, 1 through
  (`famctl_totals.txt`).
- SRR32148620: p = 0.993, expected 0.87, observed 0. Not testable on its own.
- Within a fragment, the two mates never disagree (30 and 24 fragments with both mates informative, all concordant).
**Positive control** (`scripts/famerr.py`, 4.48-4.50 Mb). For per-copy events, I took each high-quality mismatch in
one family member and checked siblings covering the same base. The sibling lacked it 31 of 51 times (SRR27313660) and
125 of 205 times (SRR32148620): about 61 % discordant. The family grouping therefore does resolve events that arise
after copying. At the low-share repeat-switch gaps of the CRISPR-like array (3971xxx, 3970963, 3194xxx, 3974365),
30 of 41 families holding a removal fragment also hold a through fragment (`famctl_totals.txt`). Caveat: those sites
are very deep, so some of these families are chance coincidences.
**Verdict: no effect (rules out a late PCR skip; P = 0.002 in SRR27313660).** It cannot see a skip made before the
first PCR copy (RT or second strand).

## T3. Second-strand / single-strand skip: read-1 orientation
In an unstranded library made from ds-cDNA, a removal made only on the second strand would appear only where read 1
comes from that strand (read 1 forward for this + gene), and the share could not exceed about 50 %. Fragment counts
(`siteanalyze.py`, Fisher in `stats1.txt`):
- SRR27313660: read 1 forward 127 removal / 13 through (0.907); read 1 reverse 94 / 7 (0.931); Fisher p = 0.64.
- SRR32148620: read 1 forward 132 / 1 (0.992); read 1 reverse 148 / 1 (0.993); p = 1.
Duplex families (same coordinates, both orientations), SRR27313660: 9; 6 agree on removal, 2 on through, 1 disagrees.
SRR32148620: 13, all agree on removal.
**Positive control** (`scripts/orispectrum.py`). The mate-specific mismatch spectrum shows a single-strand lesion
signature (8-oxoG type): C>A on mate 2 vs mate 1 is 1.01e-4 vs 4.37e-5 in DNA and 1.23e-4 vs 6.35e-5 in SRR27313660.
So these libraries keep strand-of-origin information, and a strand-confined event would show as an asymmetry.
**Verdict: no effect** on the removal. It rules out a second-strand-only skip **if** the ds-cDNA was not
whole-amplified before fragmentation. The SRA records do not name the kit, so this caveat stays open.

## T4. Strain / genomic-deletion route
(a) The 2022 DNA: 0 removal / 39 through (aligned) and 0 / 51 (k-mers). The 95 % upper bound on a deleted fraction in
the 2022 DNA culture is 0.074 (0.057 with k-mers).
(b) The sister strain B-36 (GCA_041410205.1, CP163304.1; downloaded to `indep/dl_b36/`) **lacks the 31 nt in its
genome**. B35 joined: `...AATAACATGATC|AAACCGTGTACTGG...`. B-36: `...AACAACATGGTT|AAACCCTGCACTGG...` (blastn: the
joined B35 form aligns 1-594 without a break; the B35 "with" form breaks at the stretch). So the insert is polymorphic
among L. ossiferum strains. This is the main fact supporting the strain route. Reads with B-36 junction 24-mers in the
three SAMs: 0, 0, 0. The removal reads carry B35 alleles at every flank base: no mismatch is found in more than 3 of
166 or more removal reads within ±150 nt (`scripts/snpclass.py`). Contamination by B-36 or any diverged relative is
therefore excluded. A B35-identical cell line that lost the 31 nt is not.
(c) Do the RNA cultures differ genetically from the 2022 DNA? (`scripts/varscan.py` -> `rna_fixed_variants.tsv`,
`stats_var.txt`; 384,604 positions with depth 20 or more in both RNA SAMs.) SNVs at 80 % or more in RNA with under 10 %
in DNA: **0 shared by both RNA cultures**, 6 only in 2023 (one haplotype block at 5187235-5187253, DNA 1/34-5/37), and
1 only in 2024 (2716542: 1456/1513 alt in SRR32148620, 0/870 in SRR27313660, 0/40 in DNA).
**Positive control:** 2716542 shows that the method detects a variant that rose to near fixation in one later culture
and is absent from the 2022 DNA. So the cultures do drift. But L02 is removed at at least 0.92 in both RNA cultures,
and in ~385 kb no SNV behaves like that.
(d) The second high-share site, 859933 (46 nt; B-36 *keeps* this stretch), is removed at 0.962 and 0.940 in the
same runs. Its share differs from L02's within SRR32148620 (0.993 vs 0.940, Fisher p = 1.6e-5). A single deleted
subpopulation carrying both losses would give equal shares, given equal per-cell expression.
**Verdict: weakens the strain route but cannot kill it.** A genomic loss of the 31 nt in a B35 line after the 2022 DNA
sample, sweeping before the 2023 and 2024 cultures split, fits every count here except (d). It is not testable without
DNA from the RNA cultures. RNA-library gDNA cannot be told apart in unstranded data at an expressed gene.

## T5. RT / polymerase skip over a hairpin: is it general in these libraries?
`scripts/allgaps.py` (all RNA gaps 15-500 nt, MAPQ 20 or more, anchors 8 or more, left-normalised; 170 sites) and
`scripts/hairpinscan.py`. All 31-nt windows with summed RNA depth of 30 or more, folded with ViennaRNA (min over both
strands), at least as stable as L02's own window by the same code. Windows are collapsed to non-overlapping loci, and
loci within 100 nt of L02/859933 are excluded. Removal reads = gaps of 15 nt or more inside the window ±10.
- RNA parameters (L02 window -7.0): **5,401 loci, 9 with any removal read, 0 with share 0.01 or more, 0 with share
  0.10 or more**. The maximum is 0.0076 (1191358, 2 reads / 260). 31 removal reads in total over a depth sum of
  3,364,921 (`hairpin_summary.txt`).
- DNA parameters (Mathews 2004; L02 window -3.6): 4,487 loci, 2 with any removal read, 0 at 0.01 or more.
- Detector check: the same code gives the L02 window a share of 0.946 (556 removal / min depth 32). The 31-nt window
  cannot hold the 46-nt 859933 gap (share 0 by construction). With a 46-nt span the scan counts 1,264 removal reads
  there.
**Positive control:** the genome-wide gap scan (`gapscan.tsv`) finds real low-level template-switch artefacts in these
same RNA libraries, absent from DNA. Examples are repeat-unit deletions in the CRISPR-like array: 3971102 (68 nt, 7-nt
repeat) 10/5865 and 10/35906; 3970963 (89 nt, 42-nt repeat) 2/178 and 30/732; and 1619575 (18 nt, 33-nt repeat).
Their shares are 0.0003-0.04. The scan sees RT/PCR skips; they are just 25-3000x rarer than at L02.
**Verdict: no effect for a general structure skip.** It cannot exclude a skip driven by something unique to L02.

## T6. Raw reads (unmapped included) and excised-circle probes
This tests whether filtering in the SAMs hid through reads or B-36 reads, and whether an excised 31-nt circle is
present. `scripts/rawscan.sh` streams every read with `fastq-dump --split-spot`, unmapped reads included, from the
other job's finished `data/sra3/SRR32148621.sra` ("469A", 2024, 32,915,459 spots; it is not in the SAMs) and
`SRR32148620.sra` ("469B", 32,790,182 spots). It greps both strands for 24-mer probes (`probes_labels.tsv`) and
counts with `scripts/rawcount.py`. The zeros were printed explicitly (`raw_summary.txt`).
- SRR32148621: L02 joined 895, through-left 9, through-right 7 (share 0.990). 859933: 1645 / 94 / 141 (0.921).
- SRR32148620: L02 joined 1004, through-left 25, through-right 23 (share 0.976). 859933: 3056 / 189 / 247 (0.925).
- B-36 junction probes: 0 and 0. Circle junctions of the excised L02 stretch (3 register shifts): 0 and 0. 859933
  circles: 0 and 0. Shifted control circles (±100/200/300 nt): 0 everywhere.
The independent replicate 469A reproduces the high L02 share (0.990) and the L02 > 859933 difference. The raw reads
contain no hidden through population and no B-36 reads.
**Positive control for the circle probes: none.** I found no tRNA intron or other known circle with enough RNA depth
in these libraries (tRNAs are barely covered, and no ARAGORN is installed). A 31-nt circle would also be lost in
size selection. **Verdict: raw-read share, no effect (it reproduces the claim in a third library). Circles: not
testable.** Zero circles neither supports nor argues against splicing.

## Things noticed in passing (no verdict)
- A third RNA-only removal at 3170068 (30 nt; B-36 keeps the stretch) gives 5 removal / 0 through and 5 / 1, with DNA
  0 / 16. That is too few reads to judge.
- 859933's stretch is not a strong hairpin (RNA -5.9; its first 31-nt window -2.6). Removal at high share in these
  libraries therefore does not need a stable hairpin.
- My `samtools mpileup` call created `culture/data/ref/genome.fna.fai`. I deleted it after the run; the folder is back
  to `genome.fna` + `three.fna`.

## Overall verdict
| test | verdict |
|---|---|
| T1 mapping / aligner | no effect (k-mers 0/51 DNA, 285/31 and 342/7 RNA; the artefact control is flagged) |
| T2 PCR-stage skip | no effect: excluded in SRR27313660 (expected 8.29 mixed families, observed 1, P = 0.0023) |
| T3 second-strand-only skip | no effect: read-1 orientation 0.907 vs 0.931 and 0.992 vs 0.993. This depends on the protocol (no kit named) |
| T4 strain / deletion | weakens. B-36 or a diverged contaminant is excluded. A B35 line that lost the 31 nt after 2022 is not testable |
| T5 general RT / hairpin skip | no effect: 0 of 5,401 matched-hairpin loci reach a share of 0.01; L02 = 0.946 |
| T6 raw reads / circles | raw: no effect (469A replicate 0.990). Circles: not testable (no positive control) |

**Overall: not killed.** I found no library-prep or mapping mechanism that reproduces a 92-99 % removal at this site:
- Neither PCR-stage nor second-strand skips fit the family and orientation data.
- The same libraries make structure or repeat skips only at 0.03-4 %, including at matched or stronger hairpins.

**Strongest route I could not exclude: a genomic one.** A B35 subline that precisely lost the 31 nt after the 2022 DNA
sample could have taken over before the 2023 and 2024 RNA cultures split. The insert is genuinely polymorphic: the
sister strain B-36 lacks it. And the RNA cultures do drift (2716542 became nearly fixed in one culture). Against
this:
- Across 385 kb, no SNV is fixed in both RNA cultures.
- 859933 is removed at a different share within the same library.

Neither point excludes the route. Two tests would settle it: DNA (or a no-RT control) from the 2023/2024 cultures, or
RT-PCR of the site with gDNA from those same cultures. The runner-up is an RT skip driven by something unique to L02.
The matched-hairpin scan cannot test this, though 859933 shows that a stable hairpin is not needed for high removal in
these libraries.
