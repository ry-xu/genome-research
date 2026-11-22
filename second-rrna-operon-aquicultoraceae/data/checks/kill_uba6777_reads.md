<!-- Copied from the lab record. Paths that point to files in this repository were rewritten to repo-relative ones; any other path refers to the original working folder and that file is not included (see data/README.md). -->

# Kill pass on the UBA6777 read tests (fresh subagent, 15 Nov ~01:15; saved verbatim by Deneb)
Brief: try to kill Rule 1, Rule 2 and the co-abundance reading (NOTES.md). Own code, did not import coab.py.
Scripts and tables: scratchpad kill_ux/ (scan.py, an1-an14.py, per-run .json and .dump.tsv).

## 1. Verdicts
| Claim | Verdict | Evidence |
|---|---|---|
| Rule 1 record (narrowly not met) | Survives as written; the deficit is probably noise | My code gives exactly 6356: 23, 10856: 9, controls 16/25/21/18/15/26/22/17/12/19, median 18.5, ratio 0.49. The deficit does not recur in the other runs. ERR599030: 10856 = 15 against a control median of 12.5 (1.20x). ERR599136: 2 against 4. Pooled over the three runs, 10856 is 26 against a control median of 36.5 (0.71x). This cross-run check is post hoc. |
| Rule 2 record | Weakened | NODE_1 is 100% to ux_div over all 3483 nt, and the remap tables match NOTES. But NODE_6 is a "candidate" only under remap choice (b) (2.82x). It fails under (a) (3.04x), and it fails against the full-run D used for Rule 1 (749.1 / 187.5 = 4.0x in ERR599038; 4.3x and 4.1x in the other two runs). With D taken over the core 23S only it is 3.2x, or 2.97x on the 23S body only. "Present" sits at the 3x edge under every accounting and depends on which denominator you pick. The co-abundance test, not Rule 2, carries the claim. |
| Co-abundance numbers | Survive exactly | H 249.19 / 216.30 / 95.37, N6 2195 / 1971 / 845 reads, N6/H 3.006 / 3.110 / 3.024 (fold range 1.035), D 844 / 698 / 319, D/H 0.752 / 0.717 / 0.743 (fold range 1.050). 34 decoys, 33 above 1.5-fold, closest NODE_30 at 1.394. |
| Logs | Complete | Read = Written in all three runs; minimap2 "Real time" in all three; sequences processed = 2 x spots exactly (306,384,522 / 342,278,366 / 314,283,624); SAMs end on read names near the last spot. |
| Small ERR599136 SAM | Explained | Reads at >=0.98 on MAG contigs 1.06M against 2.65M / 2.31M (host 2.6x lower there). Census hits similar in all three (0.60-0.74M). |
| Co-abundance reading ("tracks the host") | Survives, robust | Section 2. |
| "Leading candidate, not proven" wording | Survives with caveats | Sections 2 and 4. |

## 2. What I verified myself
- Filters don't change the reading: MAPQ >= 20 gives NODE_6 fold 1.105, below all 30 decoys (closest 1.701). >= 90% of read aligned: 1.106, below all 30 (closest 1.374). 23S body only (330-2800; excludes upstream spacer 1-317 with a tRNA-Ile): 1.042.
- Coverage along NODE_6 is even: 300-nt windows 2.3-3.4x H in every sample, same shape in all three runs. 3-5 reads per sample MAPQ 0. Windowed blast best hit mostly the family loci, no switch to a high-identity hit from another taxon: no sign of a chimera. NODE_6 is ~0.77-0.80 to the divergent copy, so >= 0.98 reads cannot be shared.
- Decoy comparison is fair: every decoy varies far beyond counting noise (e.g. NODE_22 fold 1.63 vs noise-only 95th pct 1.095). NODE_6's 1.035 is inside its own noise-only 95th pct (1.094). The 12 decoys with >= 500 reads in every sample have minimum fold 1.634.
- NODE_6 behaves like a host contig: across 262 host contigs >= 5 kb, fold-range median 1.077, 90th pct 1.219, 18% <= 1.035.
- Strength of the decoy argument: 0/34 within 1.39, so the 95% upper bound on a random community member tracking this well is ~8.4%. Decoys are random abundant taxa, not a co-varying close relative, so they cannot rule that out.
- Samples are separate stations (ENA): Tara 128, 138, 142, all 5 m, 0.22-3 um, WGS. (ERR599030 sample title says "NUC-RNA" but ENA lists the run as WGS; rRNA/host ratios look genomic.) Host differs mainly at station 142 (ERR599136); the other two differ by 1.15x, so the test rests largely on one contrast.
- NODE_6 is absent from every UBA6777 MAG: blast vs all five finds only 76-88% hits to the divergent copy; no physical link. JALTMC010000018 (M3003) 23S fragment matches NODE_6 100% over 109 nt in the conserved PTC region, where family locus JACKEL also matches 108/109: uninformative. The 16S 5' fragment at the edge of DKGF01000269 can't give a copy number (reads mostly soft-clipped from other organisms, 17-48x H, varying).
- Copy number: core divergent 23S at 0.88-0.97x H (single copy). NODE_6 body / core D = 2.97 / 3.26 / 2.86. Fits ~3 normal + 1 divergent if NODE_6 is the host's.
- Variant sites (caveat): NODE_6 has 18 positions with >= 25% non-reference reads (~6/kb); host contig DKGF01000210 0.2/kb. If NODE_6 is the host's, these are copy-to-copy differences. Most frequencies stable across runs (supports one multi-copy genome). One linked block (1603/1605/1637/1639) shifts ~0.47 -> 0.18 in ERR599136 at depth ~34; fits a varying strain mix better, but weak (~p 0.04 after a crude correction for overlapping mates, ~18 positions tested). Positions 149 and 167 also shift but sit in the tRNA, where other organisms' reads pile on.
- NEW: the DKGF01 divergent copy contains ~1 kb of insert sequence that reads do not support. Relative to ux_div (PCBK01) and NODE_1, DKGF01000210 carries three AT-rich inserts at 7159-7473, 7773-8143, 8289-8651 (GC as low as 0.17). Zero clean reads span the five unambiguous insert boundaries in any run (7159 ambiguous: 2 reads each in ERR599038 and ERR599030). Soft-clip breakpoints pile up at the boundaries; 23 / 27 / 6 reads span the insert-free joins. So 6356-10856 includes mis-joined sequence; D/H 0.75 is low because of the inserts; the true 23S runs at ~0.92x H. Changes no verdict (thresholds and positive control still pass), but DKGF01000210.1:6356-10856 should not be described as the divergent 23S as-is; ux_div and NODE_1 are the read-supported form.

## 3. Not checked
- Whether Parks UBA or other MAG sets from these runs contain another Myxococcota_C genome that could own NODE_6.
- A haplotype (phasing) test of the NODE_6 variants.
- Assembly read-set construction, Rule 2 placement, census content and taxonomy.
- PTC base calls on NODE_6.
- Anything on GCA_937902705 or CALELE.
- The pre-run audit of coab.py, and the recall simulation.

## 4. Suggested wording
- Rule 2: met only under remap choice (b); under (a), or against the full-run D, NODE_6 is 3.0-4.0x, outside 0.33-3x. Borderline, not positive.
- Co-abundance: keep "tracks the host"; add robustness to MAPQ >= 20, alignment-length filter, 23S body only (1.04-1.11); effectively one contrast (station 142 vs two Pacific stations); decoys don't exclude a co-varying close relative.
- Public line: "NODE_6 is the leading candidate for UBA6777's own normal 23S, probably about 3 copies. Its depth tracks the host across three Tara stations, but nothing links it physically to the host genome, and it is only 0.836 identical to the nearest family 23S." Mention copy-to-copy variants.
- Use ux_div/NODE_1 as the divergent 23S, not the DKGF01 coordinates; D ~0.92x H on the core 23S; DKGF01 copy has ~1 kb of inserted sequence reads do not support.

## Addendum (Deneb, after FACTCHECK_div23_v2.md)
The "NEW: inserts" item above is wrong: the three DKGF01000210 stretches are runs of N (scaffold gaps). Its GC came from counting N as AT. See NOTES.md "CORRECTION to the correction".

## Note (Deneb, 02:09 15 Nov)
The scratch directories cited above were wiped by a machine restart at about 02:05 15 Nov. The scripts and outputs there no longer exist; the report text is the record.
