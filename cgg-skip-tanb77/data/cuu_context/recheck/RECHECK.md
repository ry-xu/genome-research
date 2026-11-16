# Re-check: codon before a skip-like CGG in the chicken genomes (independent code)

Written 2026-11-15 03:09 by an independent re-check agent. Read: the README's Why, Sites, Expected and Predictions
sections, `results/rna_chicken/README.md` (for the site-build description), and the input tables. Not read:
`src/skip_context_chicken.py`, `result.txt`, `codon_ratios.tsv`, the README's Log, or any site-building source.
Nothing outside `recheck/` was changed.

## Method (own code)
- `recheck.py` (pass 1), `analyse.py` (rule numbers and break tests), `analyse2.py` (follow-ups). Outputs in `out/`
  (`analysis.txt`, `analysis2.txt`, `ratios_recheck.tsv`, per-site table `sites_up.tsv`). Atlas venv python `-I`,
  one process at a time; peaks 195 MB / 316 MB / 251 MB resident, 80 s / 66 s / 35 s.
- Coordinates, found by probing and then checked on every site: `cpos` is the 0-based forward index of the C;
  `+` codon = s[c:c+3], `-` codon = revcomp(s[c-2:c+1]). All 171,335 Set A + Set C codons reproduced, 0 mismatches.
- Upstream codon: the 3 nt immediately 5' of the C on the CDS strand. (On a skip walk this is always the walk's
  previous codon, because a skip consumes one base and reading restarts on the next.) All used sites have one; none is
  a stop.
- Classes: I re-derived the rule as code. Old frame = codons from the base after the CGG, no skipping, until a stop
  (capped at 60). +1 frame = codons from the first G of the CGG on a skip walk (every CGG met drops its C), until a
  stop (capped at 60). A frame that runs off the contig is undefined, so the site is unclear. This reproduces
  `site_class.tsv`'s class and both run-length columns for all 14,677 Set A sites. The main analysis uses the given
  classes. Set C and random positions use mine.
- Exclusion: 7 of the 134 chicken genomes are in `cgg_test_T.tsv`. 6 of them have Set A sites (GCA_944388605.1,
  GCA_944388625.1, GCA_944390965.1, GCA_944393615.1, GCA_946496695.1, GCA_947856825.1). That leaves **94 of 100**
  genomes with Set A sites.
- Sites used: skip-like **6,521** (in 93 genomes), in-frame-like **3,923**, unclear 3,451.
- Expected: per genome, the frequency of each codon immediately before an in-frame GGN in the pyrodigal CDS calls.
  CDSs are read from their start in frame, and the terminal stop is dropped. There are no internal stops (0 found).
  The first codon counts if codon 2 is GGN. ncRNA features are ignored. That makes 2,434,606 GGN-preceding codons in
  the 94 genomes. Expected = sum over sites of the site genome's frequency, which is the same as frequency × site
  count summed over genomes.
- Bootstrap: genomes resampled with replacement, 1,000 times, `random.Random(1)`. The units are the 93 genomes with
  ≥ 1 skip-like site (as a variant, all 94). The CTT ratio is Σobs/Σexp, with a 2.5–97.5 percentile interval.

## Rule numbers
- **P1** CTT at skip-like sites: **695 obs / 99.22 exp = 7.00** (≥ 2) → holds.
- **P2** CTT at in-frame-like sites: **98 / 59.57 = 1.65**. 7.00 > 1.65 → holds.
- **P3** skip-like ratios: CTT 7.00, TTG 2.70, TCA 2.11, CCA 2.20, GAC 2.43, so **5/5 are > 1** → holds.
- **Decision: REPLICATED.**
- The four codons depleted in the original, at skip-like sites: TTT 140/254.8 = 0.55, AAA 207/437.5 = 0.47,
  GAA 132/391.4 = 0.34, TAT 75/279.5 = 0.27. At in-frame-like sites: TTT 0.75, AAA 1.39, GAA 0.97, TAT 0.71.
- Bootstrap of CTT's skip-like ratio: **5.63–8.83** over 93 genomes (5.46–8.68 over all 94).
- For reference, the 6 excluded genomes: 445 skip-like sites, CTT 23/6.32 = 3.64.

All codons, GGN background (sorted by skip-like ratio):

| codon | aa | skip obs | skip exp | skip o/e | in-frame obs | in-frame exp | in-frame o/e |
|---|---|---:|---:|---:|---:|---:|---:|
| TGC | C | 91 | 12.2 | 7.47 | 33 | 7.0 | 4.72 |
| CTT | L | 695 | 99.2 | 7.00 | 98 | 59.6 | 1.65 |
| TCG | S | 39 | 7.5 | 5.17 | 8 | 4.4 | 1.83 |
| CCT | P | 196 | 43.8 | 4.48 | 22 | 28.0 | 0.79 |
| CCC | P | 13 | 3.0 | 4.29 | 4 | 2.1 | 1.93 |
| CTC | L | 23 | 6.1 | 3.77 | 9 | 2.8 | 3.16 |
| CAG | Q | 58 | 16.3 | 3.55 | 44 | 10.0 | 4.41 |
| CTG | L | 20 | 5.9 | 3.38 | 25 | 3.9 | 6.36 |
| GGG | G | 42 | 13.1 | 3.20 | 25 | 7.4 | 3.40 |
| TCC | S | 20 | 6.8 | 2.93 | 17 | 4.2 | 4.01 |
| TTG | L | 122 | 45.2 | 2.70 | 94 | 26.7 | 3.52 |
| GTC | V | 28 | 10.8 | 2.60 | 18 | 6.5 | 2.76 |
| ATC | I | 67 | 25.8 | 2.59 | 73 | 16.3 | 4.49 |
| CCG | P | 23 | 9.4 | 2.43 | 14 | 6.1 | 2.28 |
| GAC | D | 125 | 51.4 | 2.43 | 51 | 27.1 | 1.88 |
| CGC | R | 4 | 1.8 | 2.20 | 9 | 1.1 | 8.34 |
| CCA | P | 228 | 103.9 | 2.20 | 25 | 61.0 | 0.41 |
| TCA | S | 278 | 131.6 | 2.11 | 57 | 75.8 | 0.75 |
| AGC | S | 49 | 24.9 | 1.97 | 31 | 13.7 | 2.26 |
| TCT | S | 234 | 121.5 | 1.93 | 41 | 76.8 | 0.53 |
| GTG | V | 50 | 26.1 | 1.91 | 37 | 15.2 | 2.44 |
| AAG | K | 98 | 55.6 | 1.76 | 141 | 33.1 | 4.26 |
| AAC | N | 87 | 49.7 | 1.75 | 64 | 29.0 | 2.21 |
| AGG | R | 22 | 12.6 | 1.74 | 37 | 7.4 | 5.02 |
| GCC | A | 20 | 12.2 | 1.64 | 12 | 7.9 | 1.52 |
| ACT | T | 186 | 114.0 | 1.63 | 80 | 71.5 | 1.12 |
| CAC | H | 22 | 13.5 | 1.63 | 9 | 6.4 | 1.40 |
| GCT | A | 177 | 113.9 | 1.55 | 70 | 72.3 | 0.97 |
| GTT | V | 264 | 174.3 | 1.51 | 81 | 108.0 | 0.75 |
| GAG | E | 75 | 50.8 | 1.48 | 97 | 32.5 | 2.99 |
| ACG | T | 21 | 14.3 | 1.47 | 11 | 8.4 | 1.31 |
| TAC | Y | 24 | 20.3 | 1.18 | 47 | 12.2 | 3.86 |
| GCG | A | 19 | 16.3 | 1.17 | 14 | 9.2 | 1.52 |
| GGA | G | 190 | 165.5 | 1.15 | 97 | 98.2 | 0.99 |
| ACC | T | 12 | 10.8 | 1.11 | 14 | 6.8 | 2.05 |
| CTA | L | 66 | 61.9 | 1.07 | 27 | 35.8 | 0.75 |
| ATT | I | 293 | 275.8 | 1.06 | 212 | 167.7 | 1.26 |
| TGT | C | 95 | 89.5 | 1.06 | 42 | 53.7 | 0.78 |
| GGC | G | 28 | 30.0 | 0.93 | 21 | 18.1 | 1.16 |
| TTC | F | 20 | 21.9 | 0.92 | 62 | 11.8 | 5.24 |
| GTA | V | 162 | 186.6 | 0.87 | 73 | 112.3 | 0.65 |
| AGT | S | 116 | 139.0 | 0.83 | 46 | 83.6 | 0.55 |
| GCA | A | 184 | 239.5 | 0.77 | 84 | 137.1 | 0.61 |
| CAT | H | 58 | 76.3 | 0.76 | 43 | 48.1 | 0.89 |
| CGA | R | 5 | 6.7 | 0.75 | 6 | 3.7 | 1.60 |
| TGG | W | 35 | 47.2 | 0.74 | 37 | 28.8 | 1.28 |
| TTA | L | 235 | 344.3 | 0.68 | 141 | 207.8 | 0.68 |
| GGT | G | 117 | 174.8 | 0.67 | 46 | 106.4 | 0.43 |
| AAT | N | 275 | 422.0 | 0.65 | 164 | 258.1 | 0.64 |
| CGT | R | 11 | 17.5 | 0.63 | 15 | 10.2 | 1.47 |
| ACA | T | 180 | 299.6 | 0.60 | 58 | 175.5 | 0.33 |
| CAA | Q | 86 | 155.5 | 0.55 | 67 | 92.6 | 0.72 |
| ATG | M | 81 | 146.6 | 0.55 | 157 | 87.6 | 1.79 |
| TTT | F | 140 | 254.8 | 0.55 | 115 | 153.7 | 0.75 |
| ATA | I | 184 | 340.8 | 0.54 | 138 | 202.4 | 0.68 |
| AAA | K | 207 | 437.5 | 0.47 | 366 | 263.9 | 1.39 |
| AGA | R | 56 | 158.2 | 0.35 | 96 | 95.0 | 1.01 |
| GAA | E | 132 | 391.4 | 0.34 | 223 | 229.1 | 0.97 |
| TAT | Y | 75 | 279.5 | 0.27 | 120 | 168.1 | 0.71 |
| GAT | D | 58 | 333.1 | 0.17 | 55 | 212.7 | 0.26 |
| TAA | * | 0 | 0.0 | n/a | 0 | 0.0 | n/a |
| TAG | * | 0 | 0.0 | n/a | 0 | 0.0 | n/a |
| TGA | * | 0 | 0.0 | n/a | 0 | 0.0 | n/a |
| CGG | R | 0 | 0.8 | 0.00 | 0 | 0.5 | 0.00 |

## Break attempts

### (a) Is it carried by a few genomes, genera or genes?
- CTT skip-like sites occur in 73 of 93 genomes. Of the 33 genomes with expected ≥ 1, 31 have a ratio ≥ 2
  (median 7.09).
- The top 5 genomes by excess (obs − exp) hold 220 of the 695 CTT sites. All 5 are genus CAG-269 (ratios 15.6–20.5).
  Dropping them gives 475/86.93 = **5.46**. Dropping the top 10 gives 5.10, and the top 20 gives 3.99. The
  leave-one-genome-out range is 6.68–7.39.
- By genus (skip-like n, CTT obs/exp, ratio): CAG-269 2,511, 360/32.1 = **11.2**; Scatovivens 842, 5.78; CALXVA01
  582, 7.41; UMGS1754 383, 3.38; JAAWPK01 357, 3.17; CAG-273 278, 3.98; UMGS1663 187, 21.8; **CAG-465 184, 1/5.4 =
  0.19**; Merdicola 6.29; UMGS1864 3.48; UBA1234 4.05; CAG-793 1.94; UMGS1933 11.5; other genera pooled 4.76.
  Without genus CAG-269 the ratio is 4.99 (GGN background) or 3.06 (CGN background, see (c)).
- By GTDB family (CGN background): CAG-508 (83 genomes) 682/169.6 = 4.02. CAG-465 (9 genomes, 201 skip-like
  sites) 1/3.5 = 0.28, Poisson P(≤ 1) = 0.135. UBA1234 (1 genome) 12/3.6 = 3.31. **The excess is a result of family CAG-508
  (83 of the 94 genomes). In CAG-465 the sample is too small to show or exclude it.**
- Genes: 663 CDSs carry one CTT skip-like site and 16 carry two. With one random skip-like site per CDS
  (5,763 sites), the ratio is 7.23.
- Orthologous sites counted once per species. I grouped sites by an anchored 20-aa upstream peptide (codons −21..−2,
  leaving out the −1 codon), single linkage. A decoy test (reversed peptides) put the threshold at ≥ 10/20 identical:
  at ≥ 10/20, 24 of 200,000 random real pairs pass and 0 decoys do; at ≥ 8/20 the counts are 33 and 3. Skip-like
  sites fall into 3,791 families (largest 53); 221 CTT skip-like sites sit in 82 groups with an *identical* upstream
  peptide. Weighting each site by 1/family size gives CTT **4.86** with the GGN background and **2.84** with the CGN
  background. In-frame-like sites give 1.63 and 1.06. A family-level bootstrap (3,791 units, seed 1) gives 4.39–5.39
  with the GGN background. **Duplicated orthologous sites inflate the ratio by about 1.4× (3.93 → 2.84), but they do
  not create it.**

### (b) Does the class rule make an upstream-codon bias by itself?
- **Set C** (in-frame CGA/CGC/CGT on the same walks, classed by my re-implementation of the same rule): 998 of
  147,803 Set C sites in the 94 genomes are skip-like (0.7%), against 6,521 of 13,895 Set A sites (47%). At Set C
  skip-like sites, CTT is 15/14.21 = **1.06** with the GGN background and 15/23.69 = **0.63** with the CGN background.
  By codon (CGN background): CGA 3/10.76, CGC 4/2.94, CGT 8/9.98. Set C in-frame-like sites give 1.67 (GGN) and
  1.02 (CGN).
- **Random in-frame CDS positions** (2,000 per genome, `random.Random(7)`, classed by the same rule, background =
  all in-frame codons): skip-like 1,072 positions, CTT 11/12.24 = **0.90**; in-frame-like 1.01; positions within 30
  codons of the called stop 1.00 (n = 19,698). Random GGN positions that class as skip-like: 63, with 0 CTT against
  0.94 expected.
- **The rule does not produce a CTT excess.** Sites that pass it without a CGG show none.

### (c) Other ways it could be an artefact
- **Background context (the main caveat).** The rule's background is the codon before GGN, but the sites are the
  codon before C. CTT is more common before a C-initial codon. Under the GGN background, all Set C sites give 1.62
  and in-frame-like CGG sites give 1.65. A context-matched background (codon before in-frame CGA/CGC/CGT in the CDS
  calls, 157,058 codons) gives:
  - CTT: skip-like **695/176.8 = 3.93**, in-frame-like 98/92.8 = **1.06**, unclear 2.13.
  - Bootstrap 3.22–4.72. Dropping the top 5, 10 and 20 genomes gives 3.20, 3.11 and 2.55.
  - P3 is still 5/5 (TTG 1.25, TCA 2.72, CCA 2.18, GAC 1.93).
  - Of the depleted four, **TTT is no longer depleted (0.97)**. AAA 0.44, GAA 0.41 and TAT 0.25 stay depleted.

  So about 1.8× of the headline 7.0 (7.00/3.93) comes from G-versus-C context, not from skipping. The context-matched
  and family-collapsed ratio is about 2.8. P1 (≥ 2) still holds under every background and weighting I tried except
  the circular one below.
- A background of the codon before in-frame CGG in the CDS calls is circular. 5,711 of those 16,085 CGGs are
  skip-like sites themselves. Under it CTT gives skip-like 1.34 and in-frame-like 0.51, which still differ by about
  2.6×. The in-frame-like class is the CGG-matched control, and that is P2.
- Gene composition: using each site's own CDS codon usage as the background gives skip-like 695/83.30 = 8.34 and
  in-frame-like 1.22. Using all in-frame codons gives 9.24 and 2.20.
- Codon versus amino acid: Leu as a whole is 2.06× at skip-like sites. CTT's share of the Leu codons there is
  **0.599 against 0.176 expected**; at in-frame-like sites it is 0.249 against 0.177. Skip-like ratios by Leu codon:
  CTT 7.00, CTC 3.77, CTG 3.38, TTG 2.70, CTA 1.07, TTA 0.68. Under the CGN background, CTT ranks first among codons
  with ≥ 20 observations: CTT 3.93, TGC 2.94, CCT 2.92, CCG 2.76, TCA 2.72, CTC 2.42, CCA 2.18, TCG 2.15.
- Single nucleotide: C at position 1 is 2.43×. But ATT 1.06, GTT 1.51 and TTT 0.55 against CTT 7.00, and CTA 1.07,
  CTC 3.77 and CTG 3.38. The effect is the whole CTT codon, not one base.
- Leakage to the next codon out: CTT at codon −2 is 1.26 at skip-like sites and 1.37 at in-frame-like sites, so there
  is no excess there.
- Geometry (CTT ratio at skip-like sites, GGN background):
  - in frame inside the called CDS, ≥ 30 codons from its start: 7.12 (n = 5,529);
  - < 30 codons from its start: 10.86 (n = 182);
  - past an earlier skip or outside the call: 5.43 (n = 810);
  - +1 frame closed within 10 codons upstream, so the skip is needed: 6.79 (n = 4,631). Only 33 skip-like sites have
    the +1 frame open ≥ 60 codons upstream (0 CTT), so "the gene was in +1 all along" is not the source.
  - By the old frame's stop distance after the CGG: 0–5 codons 6.77, 6–15 codons 7.18, 16–30 codons 7.39.
- Homopolymer errors: CTT CGG has no C run (TTCGG). A C-insertion error would need a C-ending upstream codon.
  Sites whose C is a run of 1 give 7.77.

## Verdict
By the rule: **REPLICATED** (P1 7.00 ≥ 2; P2 7.00 > 1.65; P3 5/5). Bootstrap 5.63–8.83.

The break attempts did not make the CTT excess go away:
- The class rule does not create it (Set C 1.06 or 0.63; random positions 0.90).
- It is not a few genomes: 31/33 genomes are ≥ 2, and dropping the top 20 still gives 3.99.
- It is not gene composition or Leu in general.

But the headline size is inflated:
- **about 1.8×** by the GGN background, since CTT is commoner before any C-initial codon;
- **about 1.4×** by orthologous sites counted in several species.

A context-matched, family-collapsed estimate is about **2.8–3.9×**. The effect is concentrated in CAG-508,
especially genus CAG-269 at about 11×. It is not seen in CAG-465 (1 against 3.5 expected; too few sites to decide).
The TTT depletion is a background artefact (0.97 context-matched).

## What I did not check
- Whether skip-like sites are real skips. The classes come from sequence only; no reads, alignments or peptides were
  used. Any process that puts +1 shifts after CTT, biological or technical (for example an assembler or basecaller
  error specific to CTTCGG), would give the same count. I know of no such error, but sequence alone cannot exclude it.
- The other agent's numbers or code (not read by design). I did not compare.
- Families: anchored, ungapped peptide identity, not a homology search. Orthologs with an indel near the site are not
  merged, so the family collapse may under-merge.
- Context beyond the first two bases of the following codon (for example a full codon-pair model, the codon after the site, local GC).
- Phylogenetic non-independence above the family weighting (no tree-aware test).
- Random-position null only on called CDSs, not on walk stretches past a call.
- Set C and random positions are classed by my re-implementation. It matches the Set A classes exactly, but the
  producer's Set C classes do not exist to compare against.
