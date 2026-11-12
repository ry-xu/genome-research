# Red-team review (fresh subagent, finished 2026-11-11 ~08:00 EST)

Saved by Basil from the subagent's final message (verbatim below). Scripts and outputs are in this directory.

I couldn't break the claim. In all 8 clade genomes, genes end at TTA exactly where their close sister orthologs end at a canonical stop. In the 8 sister genomes this almost never happens, and none of the five attacks explains the signal away. What genome data cannot show is the mechanism: a release-factor stop and an unassigned codon that ribosomes stall on leave the same genomic pattern.

All scripts and outputs are in `scratch/redteam/`; paths below are relative to it. The panel is the 8 clade genomes plus 8 Polyangiaceae sister genomes (four DASWDH01, CAYQYG01, CAITLX01, CAYJPW01 and *Sorangium*; `panel.tsv`). Genes are my own Pyrodigal calls, which read TTA as Leu (`call_genes.py`). Orthologs are reciprocal best hits from a DIAMOND all-vs-all of the 16 genomes (`allvall.tsv`), with no Swiss-Prot involved.

## Attack 1: close-relative ortholog test — no problem found (strong support)
**Script:** `ortho_end.py`, giving `ortho/end_codon_summary.tsv` and `ortho/runon_events.tsv`.

**Method.** I took ortholog pairs where both genes are complete. In each pair one gene (the "ender") ends at a canonical stop, and its alignment reaches its C-terminus and covers at least 60% of it. I located the position matching the ender's last residue in the other gene. That gene "runs on" if it continues 5 or more residues past that point. Each gene is used once per reference group. The expected count for TTA comes from the codon mix 4–30 codons further into the run-on tails.

**Clade genomes against sister enders.** "TTA" is the number of run-on genes whose next codon is TTA. TTA ranks 1st of 61 sense codons by observed/expected in every genome.

| clade genome | TTA / run-on genes | expected | log10 p |
|---|---|---|---|
| C_020849655 (JADLFU01) | 21/260 | 0.75 | −22.7 |
| C_025057345 (JANWXX01) | 36/245 | 1.27 | −38.4 |
| C_025453595 (JANWXX01) | 63/328 | 1.39 | −78.8 |
| C_035762175 (DASTPP01) | 24/524 | 2.10 | −16.9 |
| C_046351315 (JAXLYQ01) | 32/208 | 0.96 | −36.4 |
| C_963169165 (JADLFU01) | 24/269 | 1.06 | −23.6 |
| C_964396175 (CAYNEK01) | 6/196 | 0.29 | −4.85 |
| C_964605685 (JANWXX01) | 46/309 | 1.75 | −47.4 |

- **Pooled:** 252 observed against about 9.6 expected.
- **Position:** the TTA sits at exactly the position after the ender's last residue. For example, in C_025453595 the offsets are 63 at that position, 9 one codon earlier, and 12 each one and two codons later.
- **Clade members as enders:** the signal is stronger. Counts are 49/212, 67/235, 113/364, 26/558, 36/239, 56/209, 7/244 and 111/305, with log10 p as low as −173.8.
- **Share of conserved ends:** about 4–8% of the clade's conserved gene ends are TTA (for example, C_964605685 has 765 ends at the same place as the sister and 46 at TTA).
- **Negative control (sister genomes as the runner):** TTA in that position is 0 in every case except S_035318505, which has 2 against clade enders (expected 0.40, log10 p −1.04) and 1 against sister enders.
  - TTA ranks 40–58 among codons in the other sisters.
  - S_964421485 is a high-GC sister that is itself TTA-depleted (interior TTA/CTA 0.022). It still has 0 TTA ends in 336–408 run-on genes, so rarity of TTA alone does not produce the signal.
- **Caveat:** per-genome tests share ancestry, so they are not independent across genomes. Each genome is still significant on its own.

**Is the sequence after the TTA translated?** Script `runon_coding.py`, giving `runon_all.tsv` and `runon_close.tsv`. I compared close clade-ortholog pairs and measured R3, the third-codon-position substitution rate divided by the mean of positions 1 and 2. Coding sequence has R3 well above 1; non-coding sits near 1. Intervals are 1,000 bootstraps over pairs.

| region | all clade pairs: n, R3 (95% CI) | close pairs (75–86% amino-acid identity): n, R3 (95% CI) |
|---|---|---|
| after TTA ends | 54, 0.79 (0.57–1.12) | 27, 0.87 (0.50–1.55) |
| after canonical stops | 8,278, 0.92 | 1,844, 0.81 |
| non-TTA C-terminal extensions | 933, 1.41 (1.32–1.49) | 298, 1.95 |
| before TTA ends | 60, 4.05 | 27, 6.24 |

The sequence after a TTA end looks like sequence after a stop, not like coding sequence. Power is limited because n is small.

## Attack 2: stop-codon context — no problem found
**Script:** `context.py`, giving `context.out`. Compared 252 TTA ends with 6,057 canonical stops at conserved ends and 1.3 million interior Leu codons (CTG, CTC, TTG).

- **Against interior Leu, TTA ends differ strongly:**
  - Second base after the codon: G = 85.8, p = 1.8e−18.
  - Third base after: G = 61.6, p = 2.6e−13.
  - Base before the codon: G = 67.9, p = 1.2e−14.
- **Against canonical stops, differences are small:**
  - First base after: p = 0.14.
  - Second base after: p = 0.026.
  - Third base after: p = 0.0057.
  - Base before: p = 0.0046.
  - These residuals are mild T excesses: 0.127 vs 0.073 at the base before, and 0.187 vs 0.115 at the third base after.
- **Frameshift check.** If TTA ends were really TAA or TAG shifted by one base, the next base would be A or G nearly every time. It is A or G 51% of the time, against 55% after canonical stops, so there is no sign of that.

## Attack 3: contamination or chimerism — no problem for the end signal; one contaminant contig found
**Contig GC** (`contig_gc.py`, giving `contig_gc.tsv`):
- Contigs carrying TTA events have median GC within 0.006 of other event contigs in every genome.
- Mann-Whitney p ranges from 0.022 to 0.81.
- 14 of 131 TTA-event contigs fall outside their genome's central 90% GC range, against about 10% expected.
- They are long contigs (median 10–64 kb).

**Neighbouring genes** (`neighbour_rbh.py`, giving `neighbour_rbh.out`). The share of other genes with an ortholog in another clade genome is:
- 0.730 on TTA-end contigs;
- 0.676 on canonical-end contigs;
- 0.698 genome-wide.

Contigs carrying TTA ends are, if anything, more typical of the clade than average.

**Contaminant found.** Contig CAYNEK010000112.1 in GCA_964396175 has GC 0.366 in a 0.713-GC genome. It is 13.9 kb with 14 genes, each carrying 0–13 TTAs, and has no clade orthologs (one weak hit at 31.8% identity). Tool: `contig_inspect.py`.
- It contributes none of the end events.
- It does inflate that genome's interior TTA count (10 TTA vs 398 CTA). That works against the claim.

## Attack 4: independence — no problem found
- The 137 Arm B TTA events (`tta_events.tsv`) fall in about 100 protein families. The most repeated families appear 4 times (RL19 and DHE3) or 3 times (RS3, MDH, FDX2 and FAC15).
- Only DHE3 repeats within a single genome (C_964605685).
- No contig carries more than 3 events.
- In the ortholog test each gene is counted once per reference group.

## Attack 5: other checks
**5a. Inspecting the interior TTA sites** — `interior_sites.py`, giving `interior_sites.out`. 41 clade TTAs fall inside regions aligned to a sister.
- **What the sister has at that position:** Leu 17, Phe 6, Ala 4, Ile 3, Arg 2, gap 3, and one each of Tyr, Ser, Gln, Trp, Asp and Val.
- **No interior TTA is shared by any clade ortholog.** All are single-genome cases.
- **Five sites are really shared clade ends.** The other clade orthologs stop right there; they only look interior because the sister protein is longer. These support the claim:
  - C_020849655 at residue 438 of 489;
  - C_025057345 at 316 of 375;
  - C_025453595 at 392 of 411;
  - C_963169165 at 396 of 473;
  - C_964605685 at 581 of 597.
- **Strongest genuine interior cases** (residue / protein length):
  - C_963169165, gene CAUJFL010000260.1_15621_19397_-1, residue 600 of 1,258.
    - Phe in the sister and in all 6 clade orthologs.
    - Local identity 0.80–0.95.
    - The contig looks normal: GC 0.727, and its neighbouring genes are 64–89% identical to C_020849655.
  - C_963169165, residue 120 of 403: Phe in the sister and all 7 clade orthologs.
  - C_046351315, residue 666 of 936: Phe or Met in the orthologs.
  - C_964396175, residue 288 of 858: Phe in the sister.
- The lab's read check (`results/reads/`) does not cover these sites. I did not download reads.

**5b. Are interior TTAs more common than premature canonical stops?** Scripts `stop_background.py`, `stop_background2.py`, `opportunity.py` and `event_gc.py`; outputs `stop_background*.out` and `stop_background_tests.out`.
- **Method.** I counted stop-like codons in conserved genes where the protein's homology continues in the same frame after the codon. The threshold was the 99th percentile of shuffled scores (55).
  - Premature canonical stops were searched from both the 5′ and 3′ sides.
  - As a sensitivity check, the same test recovers 78–90% of CTA sites in clade genes.
- **Mutational opportunity.** I estimated how often a single mutation creates TTA rather than a canonical stop, given each genome's codon use. The ratio is 0.13–0.31 and moves by about 20% across mutation models.
- **Result:** 14 interior TTAs with homologous continuation, against 48 such canonical stops (28 without TGA).
  - If TTA behaved like a canonical stop, about 10.4 would be expected (7.2 using the TGA-free count). That gives ratios of 1.35 (p = 0.17) and 1.95 (p = 0.015).
  - 4 of the 14 TTA cases sit on off-GC contigs, including the contaminant, against 1 of 51 canonical-stop cases. Without those 4, it is 10 observed against 10.2 expected (p = 0.56) or 7.2 (p = 0.19).
- **JAXLYQ01 is an outlier.** It has 26 premature canonical stops, 15 of them TGA aligned to assorted residues (not Trp), and the MAG is only 59% complete (`partner_residue_at_x.py`). Without it:
  - before the GC filter: 11 observed against 4.45 expected (2.5×, p = 0.006);
  - after the GC filter: 2.5× drops to 1.6–1.8× (p = 0.10–0.16).
- **Verdict: weakens it slightly, but ambiguous.** Interior TTAs in conserved genes run at about 1–1.8× the rate expected if TTA were a canonical stop. That is not significant once contaminated contigs are removed. There is room for low-level read-through, but the data don't need it.

**5c. RF2 motif.** The lab's `blind_panel_key.tsv` labels JAXLYQ01 (GCA_046351315) as having the normal SPF motif in release factor 2, not SPY. Yet it shows the TTA-end signal as strongly as the SPY members (32/208, log10 p −36.4). The lab's check 4 called this genome "not confirmed" (4 observed, 0.5 expected).
- If the claim ties the mechanism to an SPY release factor 2, this genome is either a counterexample or is missing the SPY gene; at 59% completeness that is quite possible.
- This affects only the mechanism, not whether TTA ends genes.

## Overall verdict
**Not broken.** The ends at TTA are exactly where sister orthologs stop, specific to TTA among all 61 sense codons, and absent in sisters. The sequence after them looks non-coding, and their context looks like stop context, not Leu context. They are not explained by contamination and are spread over about 100 families. The only weakening points are the single-genome interior TTAs (5a, 5b).

**Strongest remaining alternative:** TTA is an unassigned codon rather than a release-factor stop. The Leu tRNA that reads TTA would be lost, so ribosomes stall at TTA and the stalled complexes are rescued (for example by tmRNA) or released inefficiently. That would leave the same genome-level signs: genes ending at TTA, TTA avoided inside genes, non-coding tails, and some tolerated interior TTAs. Separating the two needs protein-level evidence: C-terminal peptides ending at the residue before the TTA, against tmRNA tags or read-through peptides. Short of that, check whether that tRNA is absent in every clade genome.

Two cheaper open checks:
- Read pileups at the interior sites, especially C_963169165 residue 600 (the lab already mapped ERR10643344 but kept no BAM files).
- Whether the lab's interior-depletion numbers for GCA_964396175 include contig CAYNEK010000112.1.