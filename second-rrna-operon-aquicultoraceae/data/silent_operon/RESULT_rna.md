<!-- Copied from the lab record (work/silent-operon/, 22 Nov 2026). Paths such as out/, src/, recheck/ or frozen-core/ refer to the original working folder; the small tables from out/ are in this folder (see ../README.md), the rest is not included. -->
# C_rna: is "silent" true beyond three RNA samples? (22 Nov 2026)

**Finding.** In two new Crystal Geyser RNA samples, *Ca.* Aquicultor secundus has 0 divergent 23S read pairs against 117
normal ones. One of the two is an undepleted library that is 92% rRNA. The pooled divergent 23S is at most 3.4% of the
normal per gene copy (95% bound). In the 2015 series the "viral communities" label looks like a GOLD title: the paper that
used these runs describes cell filters. No RNA was found for any other carrier species.

## Per run (first N spots; ≥98% identity; both copies plus decoys in one reference; rules in RULES.md)
| run | sample | N spots | rRNA share | 23S div : norm | 16S body div : norm | 16S spacer div : norm | verdict |
|---|---|---|---|---|---|---|---|
| SRR11558439 | CG05, 25 May 2015 (new in 2015 series) | 3 M | 28% (depleted) | **0 : 74** | 1 : 71 | 0 : 2 | SILENT (r23 ≤ 0.054) |
| SRR6241672 | CG13 metaT, library "RT-PCR" | 1 M | **92% (not depleted)** | **0 : 43** | 0 : 42 | 0 : 0 | SILENT (r23 ≤ 0.095; NM=0 tier 0 : 39, ≤ 0.105, borderline) |
| SRR6241673 | CG13 metaT, library "other" | 1 M | 26% (depleted) | 0 : 2 | 1 : 4 | 0 : 0 | uninformative |
| SRR6261448 | CG22 metaT | 1 M | 94% | 0 : 2 | 0 : 0 | 0 : 0 | uninformative (species rare) |
| SRR6261398 | CG12 metaT | 1 M | 98% | 0 : 0 | 0 : 0 | – | species absent |
| SRR7589086 | CG Site 5 dsRNA, Feb 2015 | 1 M | 0.08% | 0 : 0 | 0 : 0 | – | no rRNA (dsRNA prep) |
| SRR33092280/81/82 | AMD sediment metaT, Guangdong/Jiangxi/Anhui 2023 | 1 M each | 15/17/34% | CADDYS01: 0 : 0 | 0 : 0 | – | CADDYS01 absent (≤1 read at 94% per run) |

- Pooled new informative runs: 23S 0 : 117. At the DNA ratio (545 : 577) about 57 should be divergent. P(0) = 1.6e-34;
  r23 ≤ 0.034. With the three earlier runs, every informative Crystal Geyser RNA library gives 0 divergent 23S
  (0 : 280 in total).
- Fresh-subagent recount of SRR6241672 with its own code and its own competitive reference: 23S 0 : 44, 16S 0 : 41,
  rRNA share 92.1%. Simulated divergent reads were all recovered. It also mapped to the target copies alone, with no
  competitors: still 0 divergent pairs at ≥98%, and the best read reached 97.4% (`recheck/`).
- In the undepleted CG13 library neither 16S has spacer reads (0 : 0). In the 3M CG05 run the divergent 16S has 0 spacer
  pairs, against 4 in SRR11558543 earlier. So this data says nothing new on transcription of the divergent operon. The
  1 divergent 16S body pair in CG05 is within the precursor expectation.

## Controls
- **C1 depletion.** CG13 was sequenced twice from one sample. The "RT-PCR" library was not depleted (92% rRNA) and the
  "other" library was depleted (26%). Silence holds in the undepleted library, so rRNA depletion cannot explain the zero.
- **C2 probes** (`src/probe_windows.py`, `probe/`). The fraction of each gene in 50-nt windows that are ≥90% identical to
  E. coli / B. subtilis rRNA (the templates of probe kits):
  - A. secundus: normal 16S 21% / 30% vs divergent 0% / 0%; normal 23S 5% / 15% vs divergent 1% / 3%.
  - CADDYS01: the same pattern (normal 16S 22% / 33%, divergent 0% / 1%).

  Probe depletion would therefore remove the normal copy more than the divergent one, and favour seeing the divergent one.
- **C3 same species.** The normal 23S and 16S of A. secundus are present at 43–74 pairs in each informative run.
- **C4 recovery.** The new counter reproduces T5 on SRR11558543 (16S 5 : 141, spacer 4 : 13, 23S 0 : 125). Simulated
  pairs: ≥99% go to their own class for the A. secundus 23S, the CADDYS01 23S and 16S, and T5's 16S sims, with 0 cross.

## What was searched
- SRA "Crystal Geyser" (729 runs, `search/cg_all.csv`) plus keyword variants. That is all CG RNA: 8 JGI metaT + 1 dsRNA.
- The source BioProjects of every carrier genus (PRJNA666025, PRJNA317932, PRJNA620721, PRJNA983765, PRJNA636628,
  PRJEB36447). Only PRJNA666025 has RNA.
- Branchwater (k=21 containment over ~1.16 M SRA metagenomes) for 8 representative carrier genomes (`bw/`). It finds
  DNA only (WGS/WGA/OTHER); it does not index METATRANSCRIPTOMIC runs. Logan Search has no usable API.

## Caveats
- Heads of runs, not whole runs (the 5 GB cap; 1.8 GB used). CG13's sampling date is missing in SRA. The CG05 filter
  pore size is not stated in the paper.
- CG13 rests on 43 normal pairs (the bound is 0.095; the NM=0 tier is borderline at 0.105). Deeper CG13 would firm it up.
- One species only. A second genus could not be tested: no RNA from any carrier except A. secundus contains its reads.

**Mechanism.** The zero is not a viral-fraction or depletion artefact. In total, undepleted RNA from whole cells the
mature divergent 23S is still at most a few percent of the normal. Either it is not made into stable RNA under any
sampled condition, or it is degraded fast.
