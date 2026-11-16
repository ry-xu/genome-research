# Is the codon before a skipped CGG still CUU-rich in new genomes? (replication of an exploratory count)

Main analysis agent. Rule written 2026-11-15 02:55, before any upstream codon of a chicken site was looked at.

## Why
The post's exploratory line "the codon before the skipped C is often CUU (153 vs 26 expected)" comes from
`scratch/skip_context.py`, run on alignment-shown skips in the 64 frozen test genomes (project README, 11 Nov). CUU
is the P-site codon of the *prfB* +1 shift, where peptidyl-tRNA-Leu re-pairs one base on, so a replicated CUU excess
points to P-site slippage, one of the two candidate mechanisms. That count had a caveat: the clean-break filter drops
shifts near homopolymer runs, which could shape the depletions. This replication uses other genomes and a site set
with no such filter.

## Sites (fixed now)
- Chicken site build (`results/rna_chicken/sites.tsv`, `site_class.tsv`, rule of 12 Nov): Set A sites (in-frame CGG
  on skip walks of CDSs in primary genomes), classed **skip-like** (old frame stops ≤ 30 codons after the CGG and the
  +1 frame runs ≥ 60) or **in-frame-like** (the reverse), by sequence only.
- Genomes in `scratch/cgg_test_T.tsv` (the genomes of the original count) are excluded. How many are excluded is
  reported.
- Upstream codon: the in-frame codon just 5′ of the CGG, on the CDS strand.

## Expected (fixed now)
For each genome: the frequency of each codon immediately preceding an in-frame GGN codon in that genome's CDSs
(pyrodigal calls in `data/rna_chicken/genomes/<acc>.gff`, read in their own frame, internal stops not counted), times
the number of sites of the class in that genome; summed over genomes. That is the same background as the original.

## Predictions (fixed now)
- **P1:** CTT observed/expected ≥ 2 at skip-like sites.
- **P2:** CTT observed/expected at skip-like sites is above the same ratio at in-frame-like sites.
- **P3:** of the five codons enriched in the original (CTT, TTG, TCA, CCA, GAC), at least 4 have observed/expected > 1
  at skip-like sites.
- **Replicated** if P1, P2 and P3 all hold. **Not replicated** if P1 fails.
- Reported: every codon's ratio in both classes; the four depleted in the original (TTT, AAA, GAA, TAT); a bootstrap
  over genomes (1,000 resamples, seed 1) for CTT's ratio at skip-like sites.

## Log
- 02:55 Coordinate check before the run (3,000 Set A sites, codon reproduced only): `cpos` is the 0-based C; plus codon s[p:p+3], minus codon rc(s[p-2:p+1]). The script's first guess (1-based) was wrong and was fixed before any upstream codon was counted.
- 02:56 **Result: replicated** (`src/skip_context_chicken.py`; `result.txt`, `codon_ratios.tsv`).
  - 100 Set A genomes, 6 excluded as part of the original count's set, 94 kept: 6,521 skip-like sites (93 genomes)
    and 3,923 in-frame-like sites (94 genomes). All 10,444 site codons were reproduced.
  - **CTT before the CGG:** 695 observed at skip-like sites against 99.2 expected, ratio **7.00** (bootstrap over
    genomes, 95% 5.63-8.83). In-frame-like sites: 98 against 59.6 (1.65). P1, P2 and P3 all hold.
  - The original's other enriched codons, skip-like (in-frame-like): TTG 2.70 (3.52), TCA 2.11 (0.75), CCA 2.20
    (0.41), GAC 2.43 (1.88).
  - Its depleted codons: TTT 0.55, AAA 0.47, GAA 0.34, TAT 0.27 at skip-like sites. This site set has no clean-break
    filter, so the depletion is not that filter's doing.
  - Also high at skip-like sites (not predicted): TGC 7.47 (91/12.2), TCG 5.17, CCT 4.48, CTC 3.77. Lowest: GAT 0.17,
    TAT 0.27, GAA 0.34, AGA 0.35.
  - Caveat: these classes are by sequence (the downstream frames), not alignment-shown skips. In one genome 27 of 30
    alignment-shown skips class as skip-like (rule of 12 Nov).
  - Next: a fresh agent re-derives the CTT numbers with its own code before any text uses them.
- 03:13 **Independent re-derivation** (fresh agent, own code; `recheck/RECHECK.md`). It reproduces every number
  and the decision: **replicated**. Its break attempts:
  - **The class rule does not make the excess:** Set C skip-like sites 1.06 (CGN background 0.63), random positions 0.90.
  - **Not a few genomes:** 31 of 33 genomes with ≥ 1 expected have a ratio ≥ 2; dropping the top 20 genomes leaves 3.99.
  - **Inflated by the background:** the codon before GGN is too generous a background, since CTT is commoner before any
    C-initial codon. Against the codon before CGA/CGC/CGT, CTT is **3.93** at skip-like sites (95% 3.22-4.72) and 1.06
    at in-frame-like sites.
  - **Inflated by duplicates:** orthologous sites counted once per family bring it to **2.84** (in-frame-like 1.06).
  - **TTT** is no longer depleted under the matched background (0.97); AAA, GAA and TAT stay depleted.
  - **Where:** family CAG-508 (4.02 matched); genus CAG-269 about 11× (GGN background). CAG-465 has 1 CTT site against
    3.5 expected, too few to show or exclude.
  - Not checked: whether skip-like sites are real skips (sequence classes only).
- 03:13 The post's exploratory line was rewritten from this (backup `post_draft/post_before_cuu_20261115.md`):
  - the old 153-against-26 used the same generous background, so the bullet now gives the matched estimate "about 3 to
    4 times";
  - the CUC sentence is dropped.
  - A fresh fact-check of the new bullet comes next.
- 03:27 Fact-check 6 (`post_draft/rna/factcheck6/REPORT.md`) of the rewritten post bullet found:
  - one wrong statement: "the order's other family". The 94 genomes are CAG-508 (83), CAG-465 (10) and UBA1234 (1).
  - two misleading statements:
    - "3 to 4 times": with both corrections it is 2.8;
    - CAG-465 "too few to show it either way": it shows 1 against 3.5 expected, where CAG-508's rate would give 14.1,
      Poisson P(≤ 1) = 1e-5. So it argues against an excess of that size there, unless its sequence-classed sites are
      mostly not real skips. Only 2 of its 6 genomes called in the frozen test.
  - one false reason: CUU is commoner before CG-initial codons (2.33% before CGA/CGC/CGT, against 1.51% before GGN), not
    before C-initial codons in general (1.42%).
  - UBA1234's one genome shows the excess (12 against 3.62).
  All fixed in the post (backup `post_draft/post_before_fc6_20261115.md`).
