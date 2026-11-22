<!-- Copied from the lab record. Paths that point to files in this repository were rewritten to repo-relative ones; any other path refers to the original working folder and that file is not included (see data/README.md). -->

# Audit: "the divergent 23S of Aquicultoraceae is not detectably transcribed at Crystal Geyser"
Independent kill pass, 2026-11-14. All scripts are in `src/checks/expr_audit/`, and all outputs are in `expr_audit/work/` or printed by the scripts named below. Every number below was printed by code.

## Bottom line
- **The central observation holds.** No mature divergent 23S RNA appears in any of the three metatranscriptomes. My own recount, a read simulation and a mate-split remap all agree.
- **The claim's wording is wrong in one important way.**
  - In SRR11558425 I found **5 RNA fragments that span the divergent operon's 16S 3' end and run into the 16S-23S spacer** (mapq 60, 97-100% identity, all on the sense strand).
  - The normal operon gives 0 spacer fragments in the same sample.
  - So the divergent operon **is transcribed**. What is missing is the mature divergent 23S: it does not accumulate.
  - Better wording: "transcribed, but mature divergent 23S not detected".
- **The bound is overstated.** "≤1%" counts both mates of each pair as separate reads. Counted per fragment, the 95% bound is about 2%.
- **The RNA samples are viral-fraction preparations** ("Metatranscriptome of freshwater viral communities (RNA viruses)"), not cellular metatranscriptomes. This matters for what the RNA reads represent.

## Verdicts

| # | Test | Verdict |
|---|---|---|
| 1 | Recount from PAFs with independent code | **holds** (exact match), but the unit is alignments, not reads |
| 2 | Divergent loci present, correctly classed, not captured by decoys or normal loci | **holds** |
| 3 | Strandedness and library artefacts (short reads, adapters, antisense) | **holds** |
| 4 | Same population in RNA and DNA; does the RNA population carry the divergent operon? | **holds** for presence. It **weakens "not transcribed"**: divergent operon spacer transcripts were detected |
| 5 | Statistics | **weakened**: the bound is about 2%, not 1%. Significance itself is not in doubt |
| 6 | Wording | **weakened**: "not detectably transcribed" is contradicted, the RNA is viral-fraction, and the bound is too tight |

## 1. Recount (`src/checks/expr_audit/recount.py`)
- My code reads `loci_class_all.tsv` and keeps only the 16 loci whose contig contains one of the 8 CG tags.
- At identity ≥0.98 and alignment length ≥60 (the same filter as `expr_count.py`, which drops alen<60), I get:
  - DNA: div 982, norm 1053
  - SRR11558425: 0 / 44
  - SRR11558542: 0 / 21
  - SRR11558543: 0 / 253
- This matches the claim exactly.
- **Unit problem.** The pipeline ran `fastq-dump -Z` without splitting, so each record is R1+R2 concatenated (302 nt), and minimap2 reports each mate's half as a separate alignment. The counts are therefore alignments, mostly 2 per fragment. In distinct spots:

| Run | div spots | norm spots |
|---|---|---|
| DNA (SRR6051694) | 545 | 577 |
| SRR11558425 | 0 | 23 |
| SRR11558542 | 0 | 11 |
| SRR11558543 | 0 | 129 |

- No spot hit both classes.
- Without the alen≥60 filter, SRR11558543 has 1 divergent alignment. It is 36 nt in a region shared with other loci (s1 = s2, mapq 3), so dropping it is right.
- Non-CG divergent loci (Kmv33, DAHYRC, DKGF, DLUY) receive a few hits in SRR11558543. All are 27-38 nt at mapq ≤9, or <0.98 identity, so the filter correctly drops them.

## 2. Loci and capture (`src/checks/expr_audit/simreads.py`, `src/checks/expr_audit/simscore.py`)
**What is in ref.fa**
- All 16 CG loci are in `ref.fa` (the first 51 records equal the class table) and are correctly classed:
  - 8 norm loci, model 3-2289, 2352 nt each.
  - 5 divergent loci at model 8-2122 (2325 nt), and 3 at model 558-2122 (1649-1684 nt).
- None of the 8 CG genome accessions is among the decoy records.

**Simulated reads**
- I simulated 26,284 reads in the pipeline's own format: dUTP-stranded R1+R2 concatenated, adapter read-through for short inserts, fragment lengths 40-400, and 0% or 1% error. I mapped them with the same minimap2 command (`-x sr -c --secondary=no`) against the original `ref.fa`.
- Divergent reads are never captured by a normal locus (0 at every length) or by a decoy (1 of about 13k).
- Divergent recovery is equal to or higher than normal recovery at every fragment length:
  - 60 nt: 0.995 vs 0.925
  - 90 nt: 1.000 vs 0.970
  - at 1% error and 90 nt: 0.834 vs 0.758
- Normal reads lose a few % to decoys at 60-90 nt. Any bias therefore favours detecting divergent reads.
- 40-nt fragments fail the alen≥60 filter in both classes.

**Coverage gaps**
- Both classes lack their 23S 3' third in `ref.fa`:
  - The normal locus stops at model 2289 because the normal 23S on JAACTV010000080 carries a ~985-nt insertion with an ORF (14590-15574). The 3' part, 15575-16204, is not a locus.
  - The divergent 3' part, 20469-21177, was not called by nhmmer.
- DNA share is therefore the right normaliser.
- Region-matched shares, from model positions in `src/checks/expr_audit/stats.py`:
  - Window 8-2122: div 982 vs norm 923 alignments.
  - In the remap (section 4), the **whole** divergent 23S (069:20469-23506, including the 3' third) has **0 alignments at any identity or length**.

## 3. Library artefacts (`src/checks/expr_audit/positions.py`, `src/checks/expr_audit/strand_density.py`)
- **Strandedness.** The RNA libraries are dUTP-stranded:
  - CG normal 23S, R1 strand: 23 vs 8 (SRR11558425) and 126 vs 18 (SRR11558543).
  - Genome-wide in the remap, 5484 of 5732 single-gene R1 alignments are sense (95.7%).
  - minimap2 maps both strands, so strandedness cannot hide divergent reads.
- **Length.** RNA aligned lengths are shorter than in DNA (many 60-140 nt), but this affects both classes equally (simulation above).
- **Position.** RNA normal reads are biased toward the 5' end of the 23S. In SRR11558543, 146 of 253 fall in model 0-500. All 5 full-length divergent loci cover that region.

## 4. Population, and transcription of the divergent operon (`src/checks/expr_audit/magcount.py`, `src/checks/expr_audit/strand_density.py`)
**The remap**
- The `SRR11558425.sra` file was still present in `expr/`, complete: 2,836,406,205 bytes, the same as the S3 Content-Length. I read it in place, so nothing was downloaded.
- I remapped it **with mates split** (`fastq-dump --split-spot -Z`, minimap2 `-t 1`) against:
  - the census decoys,
  - the 35 non-CG loci,
  - the whole CG MAG GCA_009993245.1 (160 contigs).
- Wall time was 894 s. The PAF lines on the MAG are in `work/remap425.mag.paf`.

**The population is present and active**
- 7466 spots fall in protein genes at ≥0.98 identity. Non-rRNA alignment identity peaks at 1.00 (85,622 alignments) and 0.99 (18,996).
- Genes on both sides of the divergent operon on JAACTV010000069 are transcribed, for example:
  - 19610-20185: 6 spots
  - 29491-31035: 42
  - 36100-36666: 80
  - 16902-17288: 489

**Divergent 23S**
- 0 spots in the mate-split remap.
- Normal 23S in the same run: 21 spots on the 5' part, plus 7 on the 3' part after the insertion.

**New finding: spacer transcripts from the divergent operon**
- 5 spots map to JAACTV010000069:23711-24012, with mapq 60 on most mates and identity 0.974-1.000.
- From BLAST against the DASXVE copies, the divergent 16S 3' end is at about 23843 and the 23S 5' end at about 23491. So 4 of the 5 spots **run from the divergent 16S 3' end across it into the 16S-23S spacer**.
- All 5 are on the operon's sense strand.
- Background: 300-nt windows across the MAG have a median of 0 spots and a 95th percentile of 3.
- Comparison with the normal operon:
  - Its 16S 3' end has 7 spots, all of which stop at the 16S end (11825).
  - Its spacer (11822-12237) has 0.
- Interpretation:
  - The divergent operon's DNA is present in the RNA sample's population, which rules out a strain that lacks it.
  - The operon is transcribed, but its transcript stays unprocessed and the mature 23S does not accumulate.
- Caveat: this is 5 fragments in one of the three samples. SRR11558542 and SRR11558543 were not remapped.

**Divergent 16S body**
- Only 4 short alignments (58-64 nt at 25181-25245, mapq 10-12; probably a conserved stretch). Effectively absent.

**Normal 16S**
- 721 spots, but most fall in one 100-nt bin with <0.99 identity or short alignments. The reference has no 16S decoys, so these are mostly cross-taxon conserved-region reads. **No 16S count should be used without 16S decoys.**

## 5. Statistics (`src/checks/expr_audit/stats.py`)
**How the 153 arises**
- The claimed "expected ~153" is n×p with n = 318 alignments and p = 0.4826, which gives 153.5.
- Conditioning on the normal count instead (318 × p/(1−p)) gives 296.6.

**Bound by unit**

| Unit | p (DNA share) | RNA norm count | Expected div (n×p) | 95% upper bound on per-copy RNA ratio (Clopper-Pearson) |
|---|---|---|---|---|
| Alignments | 0.4826 | 318 | 153.5 | 0.0101 |
| Spots | 0.4857 | 163 | 79.2 | 0.0196 |
| Dedup fragments | 0.4882 | 156 | 76.2 | 0.0203 |
| Spots, model 558-2122 only | 0.5126 | 72 | 36.9 | 0.0404 |

- The pairs are not independent, so the ≤1% bound should be **≤~2%**. The 558-2122 row shows how far the bound moves under a stricter region match.

**Each sample alone**
- P(0 | spots, p) is 2.3e-7 (n = 23), 6.7e-4 (n = 11) and 5.5e-38 (n = 129).

**Overdispersion**
- With a within-sample design effect of ρ = 0.01, 0.05 and 0.2, the spot bound becomes 0.038, 0.093 and 0.28.
- Within one library, ρ should be near 0. Duplicates are few: 156 fragments vs 163 spots.

**Day offset and fraction**
- The DNA was collected one day after the last RNA sample (05-28 vs 05-25 to 05-27). It is also a cellular metagenome, while the RNA comes from viral-community preparations.
- The CG normal 23S is 86-1237 times rarer among 23S-mapped reads in RNA than in DNA. Per 23S-mapped read: 3.66e-3 in DNA vs 7.6e-6, 3.0e-6 and 4.3e-5 in RNA.
- p is a per-genome copy ratio, so the day offset matters only if the strain mix changed. Section 4 shows the divergent operon is present in the RNA sample's population.

**Provenance**
- I verified GCA_009993245.1 → BioSample SAMN13287385 (isolate CG_2015-19_50_66), which is "derived from … SAMN07631077". That is the DNA run's BioSample, as claimed.

## 6. Wording that overreaches
- **"not detectably transcribed"** → contradicted by the divergent-operon spacer reads. Say instead: mature divergent 23S not detected; the operon is transcribed (5 spacer-spanning fragments in SRR11558425).
- **"≤1% (95% bound)"** → about 2% per fragment. The bound also assumes the DNA copy ratio applies to the RNA samples.
- **"three metatranscriptomes from the same site"** → they are viral-community (RNA virus) metatranscriptomes, and the Aquicultor rRNA in them is incidental cellular RNA at 86-1237x lower relative abundance.
  - An RNA virus prep may preferentially keep ribosome-protected rRNA. If so, it could exaggerate the depletion of a divergent rRNA that does not assemble into ribosomes (its 16S also lacks the anti-SD; see `data/checks/kill_23s.md`).
- **"a population of this family"** → specifically *Ca.* Aquicultor secundus (Aquicultor_A). The 8 loci are one species, with near-identical MAGs.
- **"982 reads / 1053 reads"** → 982 / 1053 alignments, or 545 / 577 read pairs.

## Not checked
- I did not remap SRR11558542 or SRR11558543 with split mates, so spacer reads in those two samples are unknown.
- I did not examine the 3 `SRR647*` runs in `expr/`.
- I did not check RNA reads across the normal 23S insertion, for intron or ORF splicing.
- I did not check whether the spacer fragments are primary transcript or a processing intermediate. For example, I did not look for a 5' leader upstream of the divergent 16S, at about 25572+.
- I did not predict tRNAs in the divergent spacer. BLAST against DASXVE and Kmv33 tRNA sets found none.

## Single most likely way the claim is wrong
"Not detectably transcribed" is the wrong mechanism. The divergent operon is transcribed: sense-strand reads span its 16S-23S spacer, while the normal operon's spacer has none. What fails is the accumulation of mature divergent 23S, through degradation or failed processing or assembly, measured in viral-fraction RNA. The steady-state bound is about 2% per gene copy, not 1%.

## Files and housekeeping
- Created: `data/checks/expression_audit.md`, `src/checks/expr_audit/*.py`, `expr_audit/meta_*.xml` (SRA metadata, fetched without email), and `expr_audit/work/`. The work folder holds the remap PAF, logs, the simulation PAF, and the assembly/BioSample summaries.
- Deleted: `work/remap_ref.fa` and `work/sim.fq`. Both can be regenerated.
- I downloaded no `.sra`. The pre-existing `expr/SRR11558425.sra` (2.8 GB) was only read; it is not mine, so I did not delete it. The brief said the .sra files had been deleted, but this one is still on disk.
