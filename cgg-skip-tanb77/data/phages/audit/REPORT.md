# kill1: independent audit of the TANB77 phage CGG-skip claim

Auditor: fresh subagent, 2026-11-15. My code is in this folder. I read the phage agent's src/ and results/ only after tasks 1-2 were done.
Every number below is copied from a script log in `out/` or `t3/`.

## Verdict
**The claim survives every attack I tried.** The decisive test (task 3) is strongly positive: T1 skip sites are real +1 continuations.
- When the upstream ORF is fused to its +1 frame (read from the CGG's 2nd base), a single uninterrupted ORF from an R phage aligns across the junction in 31 of 40 T1 sites.
- The controls do much worse: T1 non-CGG gene ends span in 3 of 40, and R CGG skip sites in 8 of 40.
- The aligned homolog residue at the first +1 amino acid (Gly from GGN) is Gly in 27 of 38 cases. That fits skipping the C and reading GGN as Gly, not reading CGGN as Arg.
- The signature is not a pyrodigal quirk, not an end-of-contig effect, not prophage contamination, and not GGG-homopolymer assembly slippage. It also tracks the host *within* viral families at matched GC.

**Strongest remaining doubts:**
1. **T0 is not a valid "non-skip" control.** Its calls rest on 2/3 (CAG-465) and 2/5 (CAG-492) CGG shifts. Those counts fit the skipper rate (0.61), and against a 2-5% null they lean toward CGG skipping. So T0 phages scoring 0.55 neither confirms nor refutes genus-level host specificity. What the evidence supports is *TANB77 (CAG-508/CAG-465 etc.) hosts vs other Clostridia orders*, not skipper genus vs non-skipper genus.
2. **42 T1 phages at GC >= 32% show no signal** (28/1595 = 0.018, the same as R). Either CRISPR host calls are wrong for roughly 10% of T1, or the trait is limited to low-GC phages. The reported pooled MH (10.8, over bins 24-52) is diluted by these phages. Over bins 24-30 I get about 30.
3. **The low-GC R reference is small and one viral family dominates the target.**
   - R in bins 24-30 is only 1,307 phages, and the MH leans on bins 26/28, where R skip/O is 0.041 and 0.008.
   - vFAM-00025 alone holds 844 TANB77-host phages (1297/2002 skip/O at GC<30). Outside vFAM-00025 the within-family contrast still holds (344/646 vs 11/355).

## Task 1: rebuild groups (`t1_groups.py`, log `out/t1_groups.log`)
Rule used: CheckV quality Complete/High/Medium, genetic_code 11, CRISPR host (r207) in Clostridia, and host genus in TANB77 with call=1 in cgg_test_T.tsv (54 skipper genera, none mixed).
- T1: mine 436, run_ids 436, overlap 436 (exact).
- T0: 137 vs 137, overlap 137 (exact).
- Tf: 664 vs 664, overlap 664 (exact).
- R: my eligible pool (>= 5 kb) is 41,757. run_ids R has 6,307; 6,240 of these are in my pool. The other 67 are under 5 kb (never gene-called).
  - Bins 24 and 26 use the whole pool (51 and 255 of my 51/255; plus one phage under 5 kb).
  - Bins >= 28 are capped at 500.
- GC-bin disagreements: 0 of 7,613 overlapping ids.

## Task 2: skip/O and MH (`call_genes.py`, `skip_sites.py`, `stats.py`; logs `out/stats_{meta,single,orf}_none.log`)
pyrodigal meta mode found 14,929 CDS (>= 300 nt) whose best model used translation table 4. I excluded them (table 11 only), which probably explains why my counts are about 6% below the phage agent's.

**GC bins 24-30, meta mode:**

| Set | skip/O | Rate |
|---|---|---|
| T1 | 490/834 | 0.588 |
| T0 | 164/300 | 0.547 |
| Tf | 1025/1627 | 0.630 |
| R | 115/5835 | 0.020 |

T1 per bin: 24: 4/4; 26: 206/335; 28: 252/437; 30: 28/58. R per bin: 5/70, 12/294, 9/1115, 89/4356.

**MH RR_skip over bins 24-30** (phage bootstrap within set and bin, 1000 reps):
- T1 vs R: 29.8 (17.3-67.0)
- T0 vs R: 21.8 (11.0-81.7)
- Tf vs R: 34.1 (20.0-76.0)

Over bins 24-28 only: T1 30.6, T0 21.9.

**The phage agent's 10.8** (`results/main/score_T1_R.txt`) uses all T1 bins 24-52. That pulls in the signal-free high-GC T1 phages, so it understates the effect rather than inflating it. Its T0 vs R (21.66, bins 24-30) matches mine.

**Other gene callers:**

| Caller | T1 | R | MH T1 vs R | MH T0 vs R |
|---|---|---|---|---|
| pyrodigal single (phages >= 20 kb; 1121 skipped) | 479/793 = 0.604 | 112/5630 = 0.020 | 29.1 (16.6-65.8) | 20.1 |
| Plain longest-ORF | 534/1170 = 0.456 | 159/8127 = 0.020 | 23.8 (15.6-39.6) | 20.3 |

**Other codons** (T1 vs R, bins 24-30, meta):

| Codon | T1 | R |
|---|---|---|
| AGG | 0.035 | 0.018 |
| CGC | 0.029 | 0.006 |
| CGA | 0.020 | 0.007 |
| GGG | 0.017 | 0.010 |

These reproduce the "0.7-4.6%" statement.

**CGG usage:** unstratified pooled CGG per 1000 codons is T1 0.245 vs R 0.293. The phage agent's GC-stratified O/E is 1.066 (CI 0.97-1.20). Either way CGG is not depleted. Only about 2 CGG per phage.

## Task 3: junction-spanning homology (`t3_build.py`, `t3_eval.py`, `t3_junction_res.py`, `run_t3.sh`; logs `t3/eval_R.log`, `t3/eval_all.log`, `t3/junction_res_all.log`)
Setup:
- Seed 11. Queries: 40 T1 CGG skip sites (pool 518); 40 T1 genes ending at a non-CGG codon with the +1 frame open >= 50 codons and no CGG skip site in the gene (pool 3292 genes, one random site each); 40 R CGG skip sites from GC 24-30 (pool 115).
- Query = upstream codons 0..k-1, then the +1 frame from the codon's 2nd base up to its stop.
- DIAMOND 2.2.8 blastp: --sensitive, -p 1, --max-hsps 0, -k 200, e 1e-5. Peak RSS 0.94 / 1.15 GB.
- "Spanning" = one HSP with >= 20 aa on both sides of the junction.

Against **R-phage single ORFs** (388,977 proteins; own phage excluded):

| Class | Spanning HSP | Of those with a left-side hit | Gap run <= 2 within ±10 aa of junction |
|---|---|---|---|
| T1skip | 31/40 | 31/37 | 26 |
| T1ctrl | 3/40 | 3/24 | 1 |
| Rskip | 8/40 | 8/39 | 8 |

Against **all phages' single ORFs** (486,331 proteins; own phage excluded):

| Class | Spanning HSP | Of those with a left-side hit | Clean |
|---|---|---|---|
| T1skip | 38/40 | 38/40 | 38 |
| T1ctrl | 6/40 | 6/33 | 4 |
| Rskip | 8/40 | 8/39 | 8 |

- In spanning T1 HSPs, identity left and right of the junction is similar (per-query list in `eval_R.log`), so the +1 part is genuinely homologous, not just tolerated.
- Homolog residue aligned to the first +1 aa: **Gly in 27 of 38** T1 sites (Rskip: G in 3 of 8).

## Task 4: artefact checks (`stats.py <filter>`, `shared_clusters.py`, `shared_lowgc.py`, `context.py`; logs `out/stats_meta_*.log`, `out/shared_*.log`, `out/context_meta.log`)
All on GC bins 24-30, meta mode.

| Filter | T1 skip/O | R skip/O | MH T1 vs R (95% CI) | Notes |
|---|---|---|---|---|
| (a) one phage per UHGV vGENUS (NA = singleton) | 353/617 = 0.572 | 87/4423 = 0.020 | 30.1 (17.4-64.9) | T0 MH 23.8 |
| (a') one per vFAM (only 239 phages) | 17/31 | 5/243 | 78 (18.3-inf) | |
| (a'') one per host genus (crude) | 4/13 | 3/86 | - | Too few to test |
| (b) drop checkv_trimmed=Yes or host_markers > 0 (4,024 phages kept) | 241/401 = 0.601 | 19/810 = 0.023 | 27.0 (15.6-65.0) | |
| (c) drop CGG sites < 200 nt from a contig end | 490/829 | 115/5828 | 29.8 (18.0-63.7) | |

Extra checks:
- **Same viral family, GC < 30%, different hosts:** TANB77-host phages 344/646 vs R 11/355 (excluding vFAM-00025). The signal follows the host, not the phage lineage.
- **Not homopolymer slippage:** the base after CGG at T1 skip sites is A 308, T 130, C 27, G 25 (n=490). T1 non-skip CGGs share the same A bias (216/344).

## Task 5: T0 control validity (`t5_t0_power.py`, log `out/t5.log`)
Skippers (55 genomes) pool to 1491/2441 = 0.611 of p1 at CGG. Their minimum p1 is 8 and the median 29.

| T0 genus | p1 | cgg | Call | P(cgg <= obs \| 0.611) | P(cgg >= obs \| 2% null) | P(cgg >= obs \| 5% null) |
|---|---|---|---|---|---|---|
| CAG-465 (GCA_934120565.1) | 3 | 2 | 0 | 0.77 | 0.0012 | 0.0073 |
| CAG-492 (GCA_035388585.1) | 5 | 2 | 0 | 0.30 | 0.0038 | 0.023 |

- The share has no power to say "does not skip". If anything, it leans toward CGG skipping.
- The only contrast is a lower p1 per locus (0.0057 and 0.0095 vs skipper minimum 0.0132). That also lowers non-CGG shifts, so it may reflect assembly quality, not the decoding rule.
- Each call rests on one genome.
- The T0 phage result (0.547, MH 21.8) therefore cannot serve as a negative control.
