# Tweet check: posts/mito-uga-arginine/thread.md (21 Nov 2026)

Drafting agent's line-by-line check of the tweet draft against the outbox write-up and project data.
P = `outbox/2026-11-11-mito-uga-arginine/post.md`; W = `work/new-stops-euk/`; PC = `W/prepost/uga-arg-prepost-check.md`;
C = `work/curation/CHECK_pneumocystis.md`. Thread file: `~/.local/share/lab-dashboard-v3/posts/mito-uga-arginine/thread.md`.

## Text

| Thread claim | Source | OK? |
|---|---|---|
| P. jirovecii mitochondria "most likely" read UGA as arginine instead of tryptophan | P:1 ("most likely means arginine, not tryptophan") | yes, same strength |
| "the fungus that causes pneumocystis pneumonia" | P:1 | yes |
| No protein data yet | P:24 (no proteome in PRIDE/MassIVE/ProteomeXchange; 3 published proteomics studies have no mito-encoded protein); C:33-38 | yes |
| Main alternative: modified tryptophan tRNA reading UGA, suggested by Ma et al. 2013; would put Trp at conserved Arg sites | P:16 (Ma et al. suggested tRNA-Trp CCA is modified to read UGA); P:34 (risk: "would put tryptophan at 23 conserved arginine positions"); C:13 (Ma quote verified) | yes |
| UGA at 23 conserved arginine sites, 94–99.9% arginine | P:6 | yes |
| "across thousands of mitochondrial genomes" | P:6 ("about 4,200–6,500 genomes") | yes, rounded down to "thousands" for length |
| at none of 37 conserved tryptophan sites, which this genome encodes with UGG | P:6; PC:159 recommends "at none of the 37" wording | yes |
| eight other Pneumocystis genomes, from other mammals | P:7 (two from rats, mouse, rabbit, dog, macaque, prairie dog, squirrel) | yes |
| arginine codon 179 of 184 times at those sites | P:7; W/README.md:369-371 (AGA 113, CGA 48, CGT 9, AGG 5, CGG 4) | yes |
| closest relative, from macaques, has none of these UGAs | P:26 (macaque species, closest relative per Cissé et al. 2021, "has none"); W/results/tga_arg/fig/pneumo_cols.tsv row MT726214.1 (TGA_strictR 0, TGA_all 0); C:22 (Cissé quote verified) | yes |
| 18 more patient isolates, from one NIH group | P:8 | yes |
| UGA at its 32 arginine-position sites 564 of 576 times | P:8 (32 × 18 = 576) | yes |
| RNA from another patient reads UGA there in 7,109 of 7,298 reads | P:9; W/results/tga_arg/rna/summary.txt ("('TGA', 'R') 7298 7109") | yes. "Another patient": the RNA is a Lausanne lung-wash sample (PC:145), not one of the 18 isolates |
| Not RNA editing, which would turn UGA into CGA; only 3 reads show CGA | P:9 (3 CGA, corrected from 4 by the pre-post check) | yes. P:9 itself says "the cell does not edit UGA into an arginine codon", so "It's not RNA editing" is the source's strength |
| UGA and arginine codons swap between isolates, one way at one site and the other way at a second | P:10 (nad2:282 ref UGA, 12/18 AGA; nad1:33 ref CGA, 6/18 UGA); W/README.md:361-367 | yes |
| tRNA reading UGA unknown; genome has none with anticodon UCA | P:25; P:16 | yes |
| Prior work: proposed for a diatom mitochondrion (Yılmaz et al. 2026) and some ciliate nuclear genomes (Chen et al. 2023), not in fungi | P:17; C:14-15 (both verified; "proposed" right for Yılmaz); C:55-57 ("first case in fungi" stands) | yes. Chen 2023 covers seven spirotrich ciliates; "some ciliate nuclear genomes" fits. Lukeš 2025 calls some of Chen's cases spurious (P:17); not mentioned in the tweet |
| Fixed ending, image line | STYLE.md | yes |

Words (## 1 through the ending, without the image line): about 266.

## Figure: `figs/mito-uga-arginine/public/uga_arg_rows.png`

Script (drafting agent's scratchpad, not in the repo): `fig_uga_arg_rows.py`, adapted from `work/curation/src/fig_proline_rows.py`.
- Codons are cox1 CDS codons 430-438 (0-based) from `W/data/mito.genes.tsv.gz` (NC_020331.1) and `W/data/gbmito.genes.tsv.gz`
  (MT726214.1, *Pneumocystis* sp. 'macacae', host *Macaca mulatta*), translated with table 4: both read
  G L Q G M P R [R / UGA] I, the MPRR motif.
- Each 27-nt stretch occurs once, forward strand: NC_020331.1 (`W/data/tga_arg/NC_020331.gb`) from base 17546,
  `GGATTACAGGGTATGCCCAGATGAATT`; MT726214.1 (`W/data/gb_mito/gb_036.gb.gz`) from base 5374, `GGACTTCAAGGTATGCCTAGAAGAATT`.
- The UGA is P. jirovecii cox1 codon index 437, listed in `W/results/tga_arg/isolates.tsv` (JX499143.1 cox1 437, strict) and
  in `W/results/tga_arg/rna/summary.txt` (cox1 437 TGA, mito_Vir_cox1:457, R, strict; 253 of 253 RNA reads UGA).
- Bracket "arginine pair that binds heme": PC:140 (pairwise alignment to cow COX1 puts this UGA at R439 of the MPRR
  pair; R438/R439 bind the heme propionates); P:34.
- Caption "Arginine is at this position in 99.5% of 6,188 mitochondrial genomes": PC:140, P:34.
- Caption "tryptophan in the code these fungi are filed under": P:3 (NCBI table 4), same words as the old outbox figure legend.
- The UGA's amino-acid letter is shown as R (orange), as the proline figure shows P; the claim is "most likely".
