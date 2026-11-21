# Blinded codon-meaning inference: 8 Dipodascomycetes mitogenomes (NCBI TPA)

Generated 2026-11-11 12:13 by `/private/tmp/claude-501/-Users-ryan-lab/c532536a-a44b-4f0d-84b9-39c5fcb71590/scratchpad/blind_cun/src/report.py`. Null hypothesis: NCBI table 4 (the candidate alternative was table 3). Each codon's meaning was inferred from the data with two independent methods: (1) Codetta, which uses Pfam profiles and no relatives; (2) an ortholog-column tally against curated relatives. Organism names come from the FASTA headers and are listed for orientation only; they were not used in any inference.

## 1. Per-genome summary

**Bottom line.** Four genomes (BK074532.1, BK074262.1, BK074538.1, BK074542.1) read CTT and CTA as **proline**, not leucine. Both methods support this independently, and no other codon departs from table 4 in any genome. In these four genomes, CTC and CTG are nearly unused (0-1 codons each), so their meaning cannot be determined. The other four genomes fit table 4. BK074539.1 is a strong CUN = Leu control. BK074526.1, BK074252.1 and BK074222.1 have few CUN codons but read them as Leu wherever they can be scored, never as Pro. Every genome reads TGA as Trp and ATA as Ile, which is table 4 and rules out table 3.

| genome | modelled | deviations from table 4 | inferred code | method 1 (Codetta, CUN shown) | method 2 (ortholog columns) |
|---|---|---|---|---|---|
| BK074526.1 (Wickerhamiella azyma) | 14 genes, 3721 codons | none detected | table 4 | default: CTT/CTC/CTA/CTG = ?/?/?/?; r=0.99: ?/?/?/? | no codon deviates; CUN at conserved columns: L4 (n=4) |
| BK074532.1 (Wickerhamiella infanticola) | 14 genes, 3961 codons | CTT→P, CTA→P | table 4 except CTT and CTA = Pro; CTC/CTG not callable (CTC used 0x / 0 cons. cols, CTG used 1x / 1 cons. cols) | default: CTT/CTC/CTA/CTG = P/?/P/?; r=0.99: P/?/P/? | CTT: P at 24/24 conserved columns; CTA: P at 78/79 conserved columns (Fisher vs TTR, P vs L: CTT p=5.0e-35, CTA p=1.1e-77) |
| BK074539.1 (Wickerhamiella sp.) | 14 genes, 4105 codons | none detected | table 4 | default: CTT/CTC/CTA/CTG = L/L/L/L; r=0.99: L/L/L/L | no codon deviates; CUN at conserved columns: L171, I4, F3, Y3 ... (n=188) |
| BK074262.1 (Sugiyamaella pinicola) | 14 genes, 4080 codons | CTT→P, CTA→P | table 4 except CTT and CTA = Pro; CTC/CTG not callable (CTC used 0x / 0 cons. cols, CTG used 1x / 0 cons. cols) | default: CTT/CTC/CTA/CTG = P/?/P/?; r=0.99: P/?/P/? | CTT: P at 15/15 conserved columns; CTA: P at 41/41 conserved columns (Fisher vs TTR, P vs L: CTT p=3.6e-25, CTA p=5.8e-52) |
| BK074538.1 (Wickerhamiella siamensis) | 14 genes, 4049 codons | CTT→P, CTA→P | table 4 except CTT and CTA = Pro; CTC/CTG not callable (CTC used 0x / 0 cons. cols, CTG used 1x / 0 cons. cols) | default: CTT/CTC/CTA/CTG = P/?/P/?; r=0.99: P/?/P/? | CTT: P at 17/17 conserved columns; CTA: P at 91/92 conserved columns (Fisher vs TTR, P vs L: CTT p=2.4e-27, CTA p=5.4e-86) |
| BK074252.1 (Sugiyamaella castrensis) | 14 genes, 4114 codons | none detected | table 4 | default: CTT/CTC/CTA/CTG = ?/?/?/?; r=0.99: ?/?/?/? | no codon deviates; CUN at conserved columns: L6, I1 (n=7) |
| BK074542.1 (Wickerhamiella verensis) | 14 genes, 4003 codons | CTT→P, CTA→P | table 4 except CTT and CTA = Pro; CTC/CTG not callable (CTC used 0x / 0 cons. cols, CTG used 1x / 0 cons. cols) | default: CTT/CTC/CTA/CTG = P/?/P/?; r=0.99: P/?/P/? | CTT: P at 19/20 conserved columns; CTA: P at 88/88 conserved columns (Fisher vs TTR, P vs L: CTT p=1.2e-29, CTA p=4.1e-84) |
| BK074222.1 (Starmerella aceti) | 14 genes, 3983 codons | none detected | table 4 | default: CTT/CTC/CTA/CTG = ?/?/?/?; r=0.99: L/?/?/? | no codon deviates; CUN at conserved columns: L12, R1, F1 (n=14) |

**Codons checked in all 8 genomes.** TGA = Trp: method 2 gives BK074526 37/39, BK074532 41/41, BK074539 33/33, BK074262 42/43, BK074538 44/45, BK074252 41/42, BK074542 43/43, BK074222 40/43; Codetta calls TGA = W in every genome. That is a table-4 meaning, so it is not a deviation. ATA = Ile, not Met (table 3 would read Met): pooled conserved columns give I 297, M 4 of 536, and Codetta calls ATA = I in all 8. TAA and TAG never occur inside the modelled genes. The remaining codons are either table 4 or too rare to call (see section 6). Table 3's CUN = Thr is not supported anywhere: pooled over all genomes, CUN codons sit at 3 conserved Thr columns out of 590.

## 2. Method 2 internal controls: 2×2 tests against table-4 synonyms

Each cell counts conserved ortholog columns (≥5 relatives with a residue, ≥80% of them identical; relatives' CUN/AUA/CGA/CGC codons and stops masked) whose consensus is the named residue. Fisher's exact test is two-sided. The control row holds the codon's unambiguous table-4 synonyms, so a codon that means what table 4 says should look like its control.

**CTT vs TTA+TTG**: columns P vs L

| genome | test P / L | control P / L | odds ratio | Fisher p | test codon(s): all consensus residues |
|---|---|---|---|---|---|
| BK074526.1 | 0 / 3 | 0 / 213 | - | 1 | L3 |
| BK074532.1 | 24 / 0 | 0 / 251 | inf | 5.0e-35 | P24 |
| BK074539.1 | 0 / 59 | 0 / 83 | - | 1 | L59, I2, T1, F1 ... |
| BK074262.1 | 15 / 0 | 0 / 266 | inf | 3.6e-25 | P15 |
| BK074538.1 | 17 / 0 | 0 / 255 | inf | 2.4e-27 | P17 |
| BK074252.1 | 0 / 2 | 0 / 263 | - | 1 | L2, I1 |
| BK074542.1 | 19 / 0 | 0 / 254 | inf | 1.2e-29 | P19, Y1 |
| BK074222.1 | 0 / 8 | 0 / 229 | - | 1 | L8, R1 |

**CTA vs TTA+TTG**: columns P vs L

| genome | test P / L | control P / L | odds ratio | Fisher p | test codon(s): all consensus residues |
|---|---|---|---|---|---|
| BK074526.1 | 0 / 1 | 0 / 213 | - | 1 | L1 |
| BK074532.1 | 78 / 0 | 0 / 251 | inf | 1.1e-77 | P78, Y1 |
| BK074539.1 | 0 / 69 | 0 / 83 | - | 1 | L69, F1, M1, G1 ... |
| BK074262.1 | 41 / 0 | 0 / 266 | inf | 5.8e-52 | P41 |
| BK074538.1 | 91 / 0 | 0 / 255 | inf | 5.4e-86 | P91, Y1 |
| BK074252.1 | 0 / 4 | 0 / 263 | - | 1 | L4 |
| BK074542.1 | 88 / 0 | 0 / 254 | inf | 4.1e-84 | P88 |
| BK074222.1 | 0 / 3 | 0 / 229 | - | 1 | L3, F1 |

**CTN (all 4) vs TTA+TTG**: columns P vs L

| genome | test P / L | control P / L | odds ratio | Fisher p | test codon(s): all consensus residues |
|---|---|---|---|---|---|
| BK074526.1 | 0 / 4 | 0 / 213 | - | 1 | L4 |
| BK074532.1 | 102 / 0 | 0 / 251 | inf | 1.4e-91 | P102, Y1, T1 |
| BK074539.1 | 0 / 171 | 0 / 83 | - | 1 | L171, I4, F3, Y3 ... |
| BK074262.1 | 56 / 0 | 0 / 266 | inf | 4.2e-64 | P56 |
| BK074538.1 | 108 / 0 | 0 / 255 | inf | 2.3e-95 | P108, Y1 |
| BK074252.1 | 0 / 6 | 0 / 263 | - | 1 | L6, I1 |
| BK074542.1 | 107 / 0 | 0 / 254 | inf | 1.1e-94 | P107, Y1 |
| BK074222.1 | 0 / 12 | 0 / 229 | - | 1 | L12, R1, F1 |

**CTT vs all other table-4 Leu codons (TTR, CTC, CTA, CTG)**: columns P vs L

| genome | test P / L | control P / L | odds ratio | Fisher p | test codon(s): all consensus residues |
|---|---|---|---|---|---|
| BK074526.1 | 0 / 3 | 0 / 214 | - | 1 | L3 |
| BK074532.1 | 24 / 0 | 78 / 251 | inf | 1.3e-14 | P24 |
| BK074539.1 | 0 / 59 | 0 / 195 | - | 1 | L59, I2, T1, F1 ... |
| BK074262.1 | 15 / 0 | 41 / 266 | inf | 7.1e-13 | P15 |
| BK074538.1 | 17 / 0 | 91 / 255 | inf | 4.3e-10 | P17 |
| BK074252.1 | 0 / 2 | 0 / 267 | - | 1 | L2, I1 |
| BK074542.1 | 19 / 0 | 88 / 254 | inf | 2.7e-11 | P19, Y1 |
| BK074222.1 | 0 / 8 | 0 / 233 | - | 1 | L8, R1 |

**CTA vs all other table-4 Leu codons (TTR, CTT, CTC, CTG)**: columns P vs L

| genome | test P / L | control P / L | odds ratio | Fisher p | test codon(s): all consensus residues |
|---|---|---|---|---|---|
| BK074526.1 | 0 / 1 | 0 / 216 | - | 1 | L1 |
| BK074532.1 | 78 / 0 | 24 / 251 | inf | 2.9e-57 | P78, Y1 |
| BK074539.1 | 0 / 69 | 0 / 185 | - | 1 | L69, F1, M1, G1 ... |
| BK074262.1 | 41 / 0 | 15 / 266 | inf | 1.2e-39 | P41 |
| BK074538.1 | 91 / 0 | 17 / 255 | inf | 9.7e-69 | P91, Y1 |
| BK074252.1 | 0 / 4 | 0 / 265 | - | 1 | L4 |
| BK074542.1 | 88 / 0 | 19 / 254 | inf | 9.4e-66 | P88 |
| BK074222.1 | 0 / 3 | 0 / 238 | - | 1 | L3, F1 |

**comparator: CCN (Pro) vs TTA+TTG**: columns P vs L

| genome | test P / L | control P / L | odds ratio | Fisher p | test codon(s): all consensus residues |
|---|---|---|---|---|---|
| BK074526.1 | 102 / 0 | 0 / 213 | inf | 1.5e-85 | P102, S1, V1 |
| BK074532.1 | 1 / 0 | 0 / 251 | inf | 0.004 | I1, P1 |
| BK074539.1 | 107 / 0 | 0 / 83 | inf | 5.0e-56 | P107, Q1 |
| BK074262.1 | 49 / 0 | 0 / 266 | inf | 1.2e-58 | P49 |
| BK074538.1 | 0 / 0 | 0 / 255 | - | 1 | - |
| BK074252.1 | 106 / 0 | 0 / 263 | inf | 1.7e-95 | P106 |
| BK074542.1 | 1 / 0 | 0 / 254 | inf | 0.004 | P1 |
| BK074222.1 | 104 / 0 | 0 / 229 | inf | 3.4e-89 | P104, Y1 |

**ATA vs ATT+ATC (table 3 would make ATA Met)**: columns M vs I

| genome | test M / I | control M / I | odds ratio | Fisher p | test codon(s): all consensus residues |
|---|---|---|---|---|---|
| BK074526.1 | 3 / 40 | 1 / 40 | 3 | 0.616 | I40, L18, V14, T6 ... |
| BK074532.1 | 0 / 51 | 0 / 34 | - | 1 | I51, L18, V8, W2 ... |
| BK074539.1 | 0 / 27 | 0 / 55 | - | 1 | I27, V5, L5, W1 ... |
| BK074262.1 | 0 / 51 | 0 / 36 | - | 1 | I51, V10, L8, T2 ... |
| BK074538.1 | 0 / 26 | 0 / 61 | - | 1 | I26, L11, V9, N1 ... |
| BK074252.1 | 0 / 45 | 1 / 38 | 0 | 0.464 | I45, V11, L10, F3 ... |
| BK074542.1 | 0 / 36 | 0 / 50 | - | 1 | I36, L18, V11, Y1 ... |
| BK074222.1 | 1 / 21 | 0 / 68 | inf | 0.244 | I21, L12, V4, Y2 ... |

**TGA vs TGG (Trp vs any other residue)**: columns W vs non-W

| genome | test W / non-W | control W / non-W | odds ratio | Fisher p | test codon(s): all consensus residues |
|---|---|---|---|---|---|
| BK074526.1 | 37 / 2 | 0 / 0 | - | 1 | W37, Y1, T1 |
| BK074532.1 | 41 / 0 | 0 / 0 | - | 1 | W41 |
| BK074539.1 | 33 / 0 | 11 / 1 | inf | 0.267 | W33 |
| BK074262.1 | 42 / 1 | 1 / 0 | 0 | 1 | W42, Y1 |
| BK074538.1 | 44 / 1 | 0 / 0 | - | 1 | W44, Y1 |
| BK074252.1 | 41 / 1 | 0 / 0 | - | 1 | W41, Y1 |
| BK074542.1 | 43 / 0 | 1 / 0 | - | 1 | W43 |
| BK074222.1 | 40 / 3 | 1 / 0 | 0 | 1 | W40, F2, L1 |

## 3. Which codons fill conserved Pro and Leu columns (method 2)

In the four CUN→Pro genomes, conserved Pro columns are filled mainly by CUN and conserved Leu columns almost only by TTR. Leu is recoded onto TTA/TTG rather than lost, and Pro is carried by CUN.

| genome | Pro cols: CCN | Pro cols: CUN | Pro cols: other | Leu cols: TTR | Leu cols: CUN | Leu cols: other | CCN usage (all modelled codons) | CUN usage |
|---|---|---|---|---|---|---|---|---|
| BK074526.1 | 102 | 0 | 10 | 213 | 4 | 58 | CCT:81 CCC:0 CCA:51 CCG:0 | CTT:7 CTC:0 CTA:4 CTG:1 |
| BK074532.1 | 1 | 102 | 5 | 251 | 0 | 38 | CCT:1 CCC:0 CCA:2 CCG:0 | CTT:29 CTC:0 CTA:104 CTG:1 |
| BK074539.1 | 107 | 0 | 4 | 83 | 171 | 48 | CCT:72 CCC:8 CCA:55 CCG:5 | CTT:144 CTC:30 CTA:194 CTG:86 |
| BK074262.1 | 49 | 56 | 5 | 266 | 0 | 22 | CCT:20 CCC:2 CCA:43 CCG:2 | CTT:18 CTC:0 CTA:52 CTG:1 |
| BK074538.1 | 0 | 108 | 3 | 255 | 0 | 38 | CCT:0 CCC:0 CCA:0 CCG:0 | CTT:18 CTC:0 CTA:120 CTG:1 |
| BK074252.1 | 106 | 0 | 3 | 263 | 6 | 23 | CCT:53 CCC:1 CCA:86 CCG:2 | CTT:5 CTC:0 CTA:7 CTG:0 |
| BK074542.1 | 1 | 107 | 3 | 254 | 0 | 36 | CCT:0 CCC:0 CCA:2 CCG:0 | CTT:26 CTC:0 CTA:110 CTG:1 |
| BK074222.1 | 104 | 0 | 5 | 229 | 12 | 42 | CCT:70 CCC:0 CCA:65 CCG:1 | CTT:22 CTC:0 CTA:6 CTG:1 |

## 4. Example positions for every deviating codon

`codon #` is the codon's 1-based index in the panel gene model. The model is built from tblastn HSPs, so its first codon may not be the annotated start. `pos` is the 1-based genome coordinate span of the codon. Peptides are ±5 residues, focal residue in brackets. The panel peptide is translated with table 4, so the focal CUN prints as `l`. Lowercase marks a contested codon (CUN/AUA/CGA/CGC) shown with its table-4 letter. The relatives are the two most identical to the panel protein in that gene's alignment that have a residue in the column. `column` is the masked consensus of all relatives at that column.

### BK074532.1 (Wickerhamiella infanticola) CTT → P: 24 of 24 conserved columns (used 29× in the modelled genes; genes: atp6 1, atp8 1, cob 2, cox1 5, cox2 3, nad1 4, nad2 2, nad3 1, nad4 2, nad5 3)

| gene | codon # | pos | panel ±5 | relative 1 | relative 2 | column |
|---|---|---|---|---|---|---|
| cox1 | 128 | 23213-23215 (+) | `FWLLV[l]SLLFA` | `FWLLV[P]SMLFA` PQ436977 Entelexis sorbosivorans (83% id) | `FVVLL[P]SMLFA` NC_005972 Candida zemplinina (76% id) | P11 |
| nad1 | 186 | 13676-13678 (+) | `CIlLL[l]ITILF` | `CIPLL[P]ITLLF` PQ436977 Entelexis sorbosivorans (87% id) | `ViPLL[P]LFIiF` NC_025768 Galactomyces candidum (70% id) | P10 |
| cox2 | 10 | 10535-10537 (+) | `iSNDV[l]TlWGi` | `iFNDV[P]TPWGM` NC_005972 Candida zemplinina (81% id) | `iSNDM[P]RPWGK` PQ436977 Entelexis sorbosivorans (79% id) | P11 |
| nad5 | 9 | 36948-36950 (+) | `TLLLL[l]IIGSV` | `TLLVL[P]LiGSL` PQ436977 Entelexis sorbosivorans (81% id) | `LiiML[P]iiGSL` HG937690 Blastobotrys adeninivorans (65% id) | P10 |
| cob | 209 | 52765-52767 (+) | `RilMH[l]YFiFK` | `RiPMH[P]YFLFK` PQ436977 Entelexis sorbosivorans (76% id) | `RiPMH[P]FFLFK` NC_005972 Candida zemplinina (75% id) | P11 |
| nad4 | 212 | 45016-45018 (+) | `iLIKS[l]LFlLH` | `IMVKS[P]LFPLH` PQ436977 Entelexis sorbosivorans (86% id) | `iLVKS[P]LYPLH` HG937690 Blastobotrys adeninivorans (63% id) | P10 |
| nad2 | 408 | 35577-35579 (+) | `SFIGL[l]lLLGF` | `SLiGi[P]PLLGF` HG937690 Blastobotrys adeninivorans (52% id) | `SLIGI[P]PlTGF` NC_024095 Magnusiomyces capitatus (45% id) | P9 |
| atp6 | 91 | 56024-56026 (+) | `GTVFF[l]LIYSL` | `GSPYF[P]FIYSL` PQ436977 Entelexis sorbosivorans (76% id) | `GNKYF[P]LiYSL` NC_005972 Candida zemplinina (71% id) | P11 |

### BK074532.1 (Wickerhamiella infanticola) CTA → P: 78 of 79 conserved columns (used 104× in the modelled genes; genes: atp6 5, atp8 1, atp9 1, cob 13, cox1 11, cox2 7, cox3 7, nad1 9, nad2 2, nad3 1, nad4 11, nad5 10)

| gene | codon # | pos | panel ±5 | relative 1 | relative 2 | column |
|---|---|---|---|---|---|---|
| cob | 173 | 52657-52659 (+) | `LHYLL[l]FiiAA` | `LHYLL[P]FiiAA` PQ436977 Entelexis sorbosivorans (76% id) | `LHYLL[P]FVIAG` NC_005972 Candida zemplinina (75% id) | P11 |
| cox1 | 94 | 21681-21683 (+) | `FFMVM[l]VTMGA` | `FFLVM[P]VTMGA` PQ436977 Entelexis sorbosivorans (83% id) | `FYLVM[P]iTMGA` NC_005972 Candida zemplinina (76% id) | P11 |
| nad4 | 86 | 44638-44640 (+) | `iGiIi[l]ISiLS` | `iGiiL[P]iSILS` PQ436977 Entelexis sorbosivorans (86% id) | `iGLiL[P]LCiLG` HG937690 Blastobotrys adeninivorans (63% id) | P10 |
| nad5 | 108 | 37245-37247 (+) | `YMSHD[l]YQQRF` | `YMSHD[P]YQQRF` PQ436977 Entelexis sorbosivorans (81% id) | `YMSHD[P]HQSRF` HG937690 Blastobotrys adeninivorans (65% id) | P10 |
| nad1 | 41 | 12306-12308 (+) | `QRRLG[l]NAVGY` | `QRRLG[P]NAVGW` PQ436977 Entelexis sorbosivorans (87% id) | `QRRLG[P]NAVGY` NC_025768 Galactomyces candidum (70% id) | P10 |
| cox2 | 12 | 10541-10543 (+) | `NDVlT[l]WGiYL` | `NDVPT[P]WGMYL` NC_005972 Candida zemplinina (81% id) | `NDMPR[P]WGKHL` PQ436977 Entelexis sorbosivorans (79% id) | P11 |
| cox3 | 27 | 296-298 (+) | `TlSlW[l]iYVTF` | `KPSPW[P]iYVSF` PQ436977 Entelexis sorbosivorans (68% id) | `KPSPW[P]iYLVF` NC_005972 Candida zemplinina (60% id) | P11 |
| atp6 | 7 | 55772-55774 (+) | `YLiNS[l]LEQFE` | `YLINS[P]LEQFE` PQ436977 Entelexis sorbosivorans (76% id) | `LFiNS[P]LEQFE` NC_005972 Candida zemplinina (71% id) | P11 |

### BK074262.1 (Sugiyamaella pinicola) CTT → P: 15 of 15 conserved columns (used 18× in the modelled genes; genes: atp6 1, cob 2, cox1 1, cox3 2, nad1 2, nad2 1, nad3 1, nad4 3, nad5 2)

| gene | codon # | pos | panel ±5 | relative 1 | relative 2 | column |
|---|---|---|---|---|---|---|
| nad4 | 221 | 16954-16956 (+) | `iMiKS[l]LYlLH` | `IMVKS[P]LFPLH` PQ436977 Entelexis sorbosivorans (72% id) | `iLVKS[P]LYPLH` HG937690 Blastobotrys adeninivorans (67% id) | P10 |
| cox3 | 136 | 15576-15578 (+) | `HSSLN[l]NiEIG` | `HSSLS[P]NiEiG` NC_022155 Candida santjacobensis (69% id) | `HSSMG[P]TVEIG` NC_014613 Candida salmanticensis (67% id) | P11 |
| cob | 276 | 12405-12407 (+) | `EFYLL[l]FYAiL` | `EFYLL[P]FYAIL` PQ436977 Entelexis sorbosivorans (79% id) | `EFYLL[P]FYAiL` HG937690 Blastobotrys adeninivorans (79% id) | P11 |
| nad1 | 77 | 34778-34780 (+) | `iLVIS[l]iITLi` | `IlVIS[P]lITlI` NC_020043 Magnusiomyces magnusii (79% id) | `iLIiS[P]LiTLi` NC_024095 Magnusiomyces capitatus (78% id) | P10 |
| nad5 | 255 | 7529-7531 (+) | `LLEYT[l]TVLLi` | `LLEYS[P]iILLL` PQ436977 Entelexis sorbosivorans (68% id) | `iLEYS[P]TiLLV` NC_024095 Magnusiomyces capitatus (68% id) | P10 |
| cox1 | 72 | 26212-26214 (+) | `FFLVM[l]MTMGF` | `FFFVM[P]lSiGF` NC_022155 Candida santjacobensis (80% id) | `FFlVM[P]MAVGF` NC_020043 Magnusiomyces magnusii (79% id) | P11 |
| atp6 | 105 | 32612-32614 (+) | `LISLV[l]YNFAI` | `lFSMI[P]YNFAl` NC_060461 Saprochaete fungicola (71% id) | `LFSMI[P]YNFAi` NC_024095 Magnusiomyces capitatus (69% id) | P11 |
| nad2 | 415 | 14077-14079 (+) | `FMGVP[l]LAGFY` | `LiGiP[P]LLGFY` HG937690 Blastobotrys adeninivorans (55% id) | `LIGIP[P]lTGFY` NC_024095 Magnusiomyces capitatus (47% id) | P9 |

### BK074262.1 (Sugiyamaella pinicola) CTA → P: 41 of 41 conserved columns (used 52× in the modelled genes; genes: atp6 1, atp8 1, atp9 1, cob 7, cox1 6, cox2 6, cox3 1, nad1 5, nad2 2, nad4 6, nad5 5)

| gene | codon # | pos | panel ±5 | relative 1 | relative 2 | column |
|---|---|---|---|---|---|---|
| cob | 173 | 12096-12098 (+) | `SSVSN[l]TIQRF` | `SSVSN[P]TIQRF` PQ436977 Entelexis sorbosivorans (79% id) | `SSVAN[P]TLQRF` HG937690 Blastobotrys adeninivorans (79% id) | P11 |
| cox1 | 82 | 27624-27626 (+) | `FFYLV[l]LMLGA` | `GNYiV[P]lMLGA` NC_022155 Candida santjacobensis (80% id) | `GNYlV[P]lMlGC` NC_020043 Magnusiomyces magnusii (79% id) | P11 |
| cox2 | 95 | 1098-1100 (+) | `VLiAF[l]SFVLL` | `VLiAF[P]SFVLL` HG937690 Blastobotrys adeninivorans (82% id) | `IlIAL[P]SFVll` NC_060461 Saprochaete fungicola (79% id) | P11 |
| nad4 | 95 | 16576-16578 (+) | `iGiiL[l]iSiLS` | `iGiiL[P]iSILS` PQ436977 Entelexis sorbosivorans (72% id) | `iGLiL[P]LCiLG` HG937690 Blastobotrys adeninivorans (67% id) | P10 |
| nad1 | 39 | 34664-34666 (+) | `QRRLG[l]NAiGY` | `QRRlG[P]NAVGY` NC_020043 Magnusiomyces magnusii (79% id) | `QRRLG[P]NAiGW` NC_024095 Magnusiomyces capitatus (78% id) | P10 |
| nad5 | 224 | 7436-7438 (+) | `NSMEG[l]TlVSS` | `YSMAG[P]TPVSA` PQ436977 Entelexis sorbosivorans (68% id) | `SSMEG[P]TPVSA` NC_024095 Magnusiomyces capitatus (68% id) | P10 |
| nad2 | 228 | 13516-13518 (+) | `KIGLA[l]LHNWL` | `KIGiA[P]LHNWI` HG937690 Blastobotrys adeninivorans (55% id) | `KiGTA[P]FHNWL` NC_024095 Magnusiomyces capitatus (47% id) | P9 |
| cox3 | 146 | 15606-15608 (+) | `GAMWl[l]IGITA` | `GAMWP[P]lSiTA` NC_022155 Candida santjacobensis (69% id) | `GSMWP[P]iDIAS` NC_014613 Candida salmanticensis (67% id) | P11 |

### BK074538.1 (Wickerhamiella siamensis) CTT → P: 17 of 17 conserved columns (used 18× in the modelled genes; genes: cob 1, cox1 6, cox2 3, nad1 2, nad2 2, nad3 1, nad5 2)

| gene | codon # | pos | panel ±5 | relative 1 | relative 2 | column |
|---|---|---|---|---|---|---|
| cox1 | 95 | 5365-5367 (+) | `FFMVM[l]iTMGA` | `FFLVM[P]VTMGA` PQ436977 Entelexis sorbosivorans (84% id) | `FYLVM[P]iTMGA` NC_005972 Candida zemplinina (79% id) | P11 |
| cox2 | 10 | 31192-31194 (+) | `iSNDV[l]TlWGi` | `iFNDV[P]TPWGM` NC_005972 Candida zemplinina (82% id) | `iSNDM[P]RPWGK` PQ436977 Entelexis sorbosivorans (81% id) | P11 |
| nad1 | 185 | 33415-33417 (+) | `IWYCV[l]LLlIT` | `IWYCI[P]LLPIT` PQ436977 Entelexis sorbosivorans (83% id) | `VWlFV[P]lFPlF` NC_020043 Magnusiomyces magnusii (70% id) | P10 |
| nad5 | 9 | 9544-9546 (+) | `TLLLL[l]IiGSI` | `TLLVL[P]LiGSL` PQ436977 Entelexis sorbosivorans (81% id) | `LiiML[P]iiGSL` HG937690 Blastobotrys adeninivorans (65% id) | P10 |
| nad2 | 407 | 8096-8098 (+) | `SFMGL[l]lLLGF` | `SLiGi[P]PLLGF` HG937690 Blastobotrys adeninivorans (50% id) | `SLIGI[P]PlTGF` NC_024095 Magnusiomyces capitatus (44% id) | P9 |
| cob | 271 | 20337-20339 (+) | `llAIV[l]EFYLL` | `PQAIV[P]EFYLL` PQ436977 Entelexis sorbosivorans (87% id) | `PAAIV[P]EFYLL` NC_005972 Candida zemplinina (82% id) | P11 |
| nad3 | 54 | 8651-8653 (+) | `ISSLL[l]YSLVL` | `VTSLL[P]YSLAL` NC_025768 Galactomyces candidum (64% id) | `LSSIL[P]FSLNL` HG937690 Blastobotrys adeninivorans (64% id) | P9 |
| cox1 | 288 | 5944-5946 (+) | `AYSKK[l]IFGKI` | `VYSKK[P]iFGKi` PQ436977 Entelexis sorbosivorans (84% id) | `VYAKK[P]iFGKi` NC_005972 Candida zemplinina (79% id) | P11 |

### BK074538.1 (Wickerhamiella siamensis) CTA → P: 91 of 92 conserved columns (used 120× in the modelled genes; genes: atp6 6, atp8 3, atp9 1, cob 16, cox1 13, cox2 7, cox3 7, nad1 11, nad2 2, nad3 1, nad4 13, nad5 11)

| gene | codon # | pos | panel ±5 | relative 1 | relative 2 | column |
|---|---|---|---|---|---|---|
| cob | 155 | 19989-19991 (+) | `LVTAV[l]YLGQA` | `LVTAV[P]YLGQS` PQ436977 Entelexis sorbosivorans (87% id) | `LVTAi[P]YLGQA` NC_005972 Candida zemplinina (82% id) | P11 |
| cox1 | 107 | 5401-5403 (+) | `GNYLV[l]LMIGT` | `ANYLV[P]LMIGT` PQ436977 Entelexis sorbosivorans (84% id) | `ANYLV[P]LQiGS` NC_005972 Candida zemplinina (79% id) | P11 |
| nad4 | 95 | 16217-16219 (+) | `iGIII[l]LCILS` | `iGiiL[P]iSILS` PQ436977 Entelexis sorbosivorans (84% id) | `iGLiL[P]LCiLG` HG937690 Blastobotrys adeninivorans (63% id) | P10 |
| nad1 | 43 | 32989-32991 (+) | `QRRiG[l]NATGY` | `QRRLG[P]NAVGW` PQ436977 Entelexis sorbosivorans (83% id) | `QRRlG[P]NAVGY` NC_020043 Magnusiomyces magnusii (70% id) | P10 |
| nad5 | 108 | 9841-9843 (+) | `YMSHD[l]YQQRF` | `YMSHD[P]YQQRF` PQ436977 Entelexis sorbosivorans (81% id) | `YMSHD[P]HQSRF` HG937690 Blastobotrys adeninivorans (65% id) | P10 |
| cox2 | 12 | 31198-31200 (+) | `NDVlT[l]WGiYL` | `NDVPT[P]WGMYL` NC_005972 Candida zemplinina (82% id) | `NDMPR[P]WGKHL` PQ436977 Entelexis sorbosivorans (81% id) | P11 |
| cox3 | 27 | 22167-22169 (+) | `AlSlW[l]iYTAF` | `KPSPW[P]iYVSF` PQ436977 Entelexis sorbosivorans (69% id) | `KPSPW[P]iYLVF` NC_005972 Candida zemplinina (61% id) | P11 |
| atp6 | 7 | 21114-21116 (+) | `FLiSS[l]LEQFE` | `YLINS[P]LEQFE` PQ436977 Entelexis sorbosivorans (75% id) | `LFiNS[P]LEQFE` NC_005972 Candida zemplinina (69% id) | P11 |

### BK074542.1 (Wickerhamiella verensis) CTT → P: 19 of 20 conserved columns (used 26× in the modelled genes; genes: atp8 1, cob 2, cox1 6, cox2 2, nad1 2, nad2 2, nad3 1, nad4 2, nad5 1)

| gene | codon # | pos | panel ±5 | relative 1 | relative 2 | column |
|---|---|---|---|---|---|---|
| cox1 | 95 | 26650-26652 (+) | `FFMVM[l]iTMGA` | `FFLVM[P]VTMGA` PQ436977 Entelexis sorbosivorans (86% id) | `FYLVM[P]iTMGA` NC_005972 Candida zemplinina (80% id) | P11 |
| cox2 | 10 | 22162-22164 (+) | `iSNDi[l]TlWGi` | `iFNDV[P]TPWGM` NC_005972 Candida zemplinina (79% id) | `iSNDM[P]RPWGK` PQ436977 Entelexis sorbosivorans (78% id) | P11 |
| cob | 223 | 13311-13313 (+) | `RilMH[l]YFiFK` | `RiPMH[P]YFLFK` PQ436977 Entelexis sorbosivorans (88% id) | `RiPMH[P]FFLFK` NC_005972 Candida zemplinina (83% id) | P11 |
| nad1 | 41 | 24077-24079 (+) | `QRRLG[l]NAVGY` | `QRRLG[P]NAVGW` PQ436977 Entelexis sorbosivorans (88% id) | `QRRlG[P]NAVGY` NC_020043 Magnusiomyces magnusii (70% id) | P10 |
| nad4 | 228 | 11683-11685 (+) | `LHIWL[l]IVHSE` | `LHiWL[P]VVHSE` PQ436977 Entelexis sorbosivorans (84% id) | `LHIWL[P]VVHSE` HG937690 Blastobotrys adeninivorans (63% id) | P10 |
| nad2 | 407 | 3709-3711 (+) | `SFiGL[l]lLLGF` | `SLiGi[P]PLLGF` HG937690 Blastobotrys adeninivorans (52% id) | `SLIGI[P]PlTGF` NC_024095 Magnusiomyces capitatus (45% id) | P9 |
| atp8 | 33 | 14069-14071 (+) | `SKFVL[l]RVLSI` | `SKFVL[P]RVLSI` NC_005972 Candida zemplinina (89% id) | `SKFiL[P]RiLNi` NC_025768 Galactomyces candidum (71% id) | P10 |
| nad5 | 9 | 5075-5077 (+) | `TLLLL[l]IVGSi` | `TLLVL[P]LiGSL` PQ436977 Entelexis sorbosivorans (82% id) | `LiiML[P]iiGSL` HG937690 Blastobotrys adeninivorans (64% id) | P10 |

### BK074542.1 (Wickerhamiella verensis) CTA → P: 88 of 88 conserved columns (used 110× in the modelled genes; genes: atp6 6, atp8 1, atp9 1, cob 15, cox1 13, cox2 8, cox3 7, nad1 11, nad2 2, nad3 1, nad4 11, nad5 12)

| gene | codon # | pos | panel ±5 | relative 1 | relative 2 | column |
|---|---|---|---|---|---|---|
| cob | 155 | 13107-13109 (+) | `LVTAV[l]YLGQA` | `LVTAV[P]YLGQS` PQ436977 Entelexis sorbosivorans (88% id) | `LVTAi[P]YLGQA` NC_005972 Candida zemplinina (83% id) | P11 |
| cox1 | 107 | 28134-28136 (+) | `GNYLV[l]LMIGT` | `ANYLV[P]LMIGT` PQ436977 Entelexis sorbosivorans (86% id) | `ANYLV[P]LQiGS` NC_005972 Candida zemplinina (80% id) | P11 |
| nad5 | 108 | 5372-5374 (+) | `YMSHD[l]YQQRF` | `YMSHD[P]YQQRF` PQ436977 Entelexis sorbosivorans (82% id) | `YMSHD[P]HQSRF` HG937690 Blastobotrys adeninivorans (64% id) | P10 |
| nad1 | 67 | 24646-24648 (+) | `KEiLi[l]KEANA` | `KEIii[P]KESNG` PQ436977 Entelexis sorbosivorans (88% id) | `KEIIl[P]KESNK` NC_020043 Magnusiomyces magnusii (70% id) | P10 |
| nad4 | 93 | 11278-11280 (+) | `iGiIi[l]VSiLS` | `iGiiL[P]iSILS` PQ436977 Entelexis sorbosivorans (84% id) | `iGLiL[P]LCiLG` HG937690 Blastobotrys adeninivorans (63% id) | P10 |
| cox2 | 12 | 22168-22170 (+) | `NDilT[l]WGiYL` | `NDVPT[P]WGMYL` NC_005972 Candida zemplinina (79% id) | `NDMPR[P]WGKHL` PQ436977 Entelexis sorbosivorans (78% id) | P11 |
| cox3 | 27 | 15311-15313 (+) | `TlSlW[l]iYVTF` | `KPSPW[P]iYVSF` PQ436977 Entelexis sorbosivorans (69% id) | `KPSPW[P]iYLVF` NC_005972 Candida zemplinina (62% id) | P11 |
| atp6 | 7 | 14271-14273 (+) | `FLiSS[l]LEQFE` | `YLINS[P]LEQFE` PQ436977 Entelexis sorbosivorans (77% id) | `LFiNS[P]LEQFE` NC_005972 Candida zemplinina (71% id) | P11 |

## 5. Method 1: Codetta per-codon results

Codetta v2.0 run as cloned (`codetta.py <genome>.fna -m`, so mitochondrial Pfams are included; defaults e-value 1e-10, probability threshold 0.9999, max_fraction 0.01). A second pass (`codetta_infer.py -m -r 0.99`) reuses the same alignments with a relaxed probability threshold. `N` = Pfam consensus columns where the codon aligned. logP values are Codetta's log decoding probabilities; 0 marks the most probable meaning, and the more negative a value, the less probable that meaning. Codetta compares against the standard code, so its TGA = W 'difference' is table 4's own meaning.

| genome | total aligned codon-columns | '?' sense codons (default) | deviations from table 4 (default) | deviations from table 4 (r=0.99) | default inferred code (TCAG order) |
|---|---|---|---|---|---|
| BK074526.1 | 2339 | 31 | none | none | `FFL?S?S?YY??C?W?????P?P?H?Q?R???I?IMT?T?N?K?S?R?V?V?A?A?D?E?G?G?` |
| BK074532.1 | 2224 | 28 | CTT=P, CTA=P | CTT=P, CTA=P | `FFL?S?S?Y???C?W?P?P?????HHQ?R???IIIMT?T?NNK?S?R?V?V?A?A?D?EEG?G?` |
| BK074539.1 | 2648 | 13 | none | none | `FFLLS?S?YY??C?WWLLLLPPP?HHQQR???IIIMT?T?NNKK?SR?V?VVAAAADDEEG?GG` |
| BK074262.1 | 2345 | 26 | CTT=P, CTA=P | CTT=P, CTA=P | `FFL?S?S?YY??C?W?P?P?P?P?HHQ?R???I?IMT?T?N?K?S?R?V?V?AAA?D?E?G?GG` |
| BK074538.1 | 2249 | 27 | CTT=P, CTA=P | CTT=P, CTA=P | `FFL?S?S?YY??C?W?P?P?????HHQ?R???IIIMT?T?NNK?S?R?V?V?A?A?D?EEG?G?` |
| BK074252.1 | 2513 | 31 | none | none | `FFL?S?S?Y???C?W?????P?P?HHQ?R???I?IMT?T?N?K?S?R?V?V?A?A?D?E?G?G?` |
| BK074542.1 | 2209 | 27 | CTT=P, CTA=P | CTT=P, CTA=P | `FFL?S?S?YY??C?W?P?P?????HHQ?R???IIIMT?T?NNK?S?R?V?V?A?A?D?EEG?G?` |
| BK074222.1 | 2450 | 28 | none | none | `FFL?S?S?YY??C?W?????P?P?HHQ?R???IIIMT?T?N?K?S?R?V?V?AAA?D?E?G?G?` |

Key codons. Each cell gives `inference [relaxed inference] N; logP(L) / logP(P) / logP(T)` for CUN, `logP(I) / logP(M)` for ATA, and `logP(W)` for TGA.

| genome | CTT | CTC | CTA | CTG | ATA | TGA |
|---|---|---|---|---|---|---|
| BK074526.1 | ? [?] N=4; -0.7 / -25.0 / -5.6 | ? [?] N=0; -3.0 / -3.0 / -3.0 | ? [?] N=1; -3.2 / -4.7 / -3.8 | ? [?] N=1; -1.3 / -6.1 / -5.1 | I [I] N=205; 0.0 / -181.4 | W [W] N=37; 0.0 |
| BK074532.1 | P [P] N=17; -65.0 / 0.0 / -71.4 | ? [?] N=0; -3.0 / -3.0 / -3.0 | P [P] N=50; -248.9 / 0.0 / -203.2 | ? [?] N=0; -3.0 / -3.0 / -3.0 | I [I] N=193; 0.0 / -155.8 | W [W] N=34; 0.0 |
| BK074539.1 | L [L] N=93; 0.0 / -355.7 / -177.3 | L [L] N=13; -0.0 / -55.3 / -26.0 | L [L] N=124; 0.0 / -457.6 / -230.0 | L [L] N=65; 0.0 / -230.7 / -86.5 | I [I] N=100; 0.0 / -99.2 | W [W] N=35; 0.0 |
| BK074262.1 | P [P] N=10; -47.2 / 0.0 / -38.4 | ? [?] N=0; -3.0 / -3.0 / -3.0 | P [P] N=32; -167.8 / 0.0 / -141.9 | ? [?] N=0; -3.0 / -3.0 / -3.0 | I [I] N=203; 0.0 / -168.2 | W [W] N=36; 0.0 |
| BK074538.1 | P [P] N=8; -44.1 / 0.0 / -31.4 | ? [?] N=0; -3.0 / -3.0 / -3.0 | P [P] N=63; -312.3 / 0.0 / -274.8 | ? [?] N=0; -3.0 / -3.0 / -3.0 | I [I] N=115; 0.0 / -104.5 | W [W] N=37; 0.0 |
| BK074252.1 | ? [?] N=5; -0.6 / -18.9 / -7.2 | ? [?] N=0; -3.0 / -3.0 / -3.0 | ? [?] N=6; -0.3 / -19.2 / -11.1 | ? [?] N=0; -3.0 / -3.0 / -3.0 | I [I] N=231; 0.0 / -196.9 | W [W] N=40; 0.0 |
| BK074542.1 | P [P] N=10; -46.8 / 0.0 / -40.1 | ? [?] N=0; -3.0 / -3.0 / -3.0 | P [P] N=59; -297.1 / 0.0 / -256.4 | ? [?] N=0; -3.0 / -3.0 / -3.0 | I [I] N=177; 0.0 / -159.6 | W [W] N=37; 0.0 |
| BK074222.1 | ? [L] N=9; -0.0 / -37.4 / -20.7 | ? [?] N=0; -3.0 / -3.0 / -3.0 | ? [?] N=4; -1.5 / -11.2 / -4.0 | ? [?] N=0; -3.0 / -3.0 / -3.0 | I [I] N=127; 0.0 / -89.6 | W [W] N=42; 0.0 |

Per-codon Codetta inference with N aligned for all 64 codons (default run):

| codon | t4 | BK074526 | BK074532 | BK074539 | BK074262 | BK074538 | BK074252 | BK074542 | BK074222 |
|---|---|---|---|---|---|---|---|---|---|
| TTT | F | F (140) | F (109) | F (71) | F (126) | F (73) | F (142) | F (86) | F (108) |
| TTC | F | F (33) | F (67) | F (146) | F (31) | F (106) | F (37) | F (89) | F (83) |
| TTA | L | L (342) | L (360) | L (126) | L (381) | L (358) | L (387) | L (355) | L (381) |
| TTG | L | ? (0) | ? (3) | L (19) | ? (5) | ? (2) | ? (9) | ? (3) | ? (3) |
| TCT | S | S (37) | S (33) | S (73) | S (38) | S (18) | S (31) | S (22) | S (65) |
| TCC | S | ? (0) | ? (2) | ? (22) | ? (4) | ? (0) | ? (0) | ? (0) | ? (0) |
| TCA | S | S (86) | S (72) | S (32) | S (55) | S (93) | S (73) | S (84) | S (53) |
| TCG | S | ? (0) | ? (0) | ? (7) | ? (1) | ? (0) | ? (1) | ? (0) | ? (0) |
| TAT | Y | Y (124) | Y (122) | Y (77) | Y (122) | Y (107) | Y (140) | Y (105) | Y (135) |
| TAC | Y | Y (7) | ? (8) | Y (73) | Y (9) | Y (20) | ? (6) | Y (17) | Y (5) |
| TAA | * | ? (0) | ? (0) | ? (0) | ? (0) | ? (0) | ? (0) | ? (0) | ? (0) |
| TAG | * | ? (0) | ? (0) | ? (0) | ? (0) | ? (0) | ? (0) | ? (0) | ? (0) |
| TGT | C | C (22) | C (22) | C (28) | C (25) | C (24) | C (28) | C (24) | C (24) |
| TGC | C | ? (1) | ? (0) | ? (2) | ? (0) | ? (1) | ? (0) | ? (0) | ? (0) |
| TGA | W | W (37) | W (34) | W (35) | W (36) | W (37) | W (40) | W (37) | W (42) |
| TGG | W | ? (0) | ? (0) | W (14) | ? (2) | ? (0) | ? (0) | ? (0) | ? (0) |
| CTT | L | ? (4) | P (17) | L (93) | P (10) | P (8) | ? (5) | P (10) | ? (9) |
| CTC | L | ? (0) | ? (0) | L (13) | ? (0) | ? (0) | ? (0) | ? (0) | ? (0) |
| CTA | L | ? (1) | P (50) | L (124) | P (32) | P (63) | ? (6) | P (59) | ? (4) |
| CTG | L | ? (1) | ? (0) | L (65) | ? (0) | ? (0) | ? (0) | ? (0) | ? (0) |
| CCT | P | P (58) | ? (0) | P (53) | P (9) | ? (0) | P (31) | ? (0) | P (47) |
| CCC | P | ? (0) | ? (0) | P (6) | ? (2) | ? (0) | ? (1) | ? (0) | ? (0) |
| CCA | P | P (38) | ? (1) | P (34) | P (26) | ? (0) | P (55) | ? (0) | P (45) |
| CCG | P | ? (0) | ? (0) | ? (2) | ? (0) | ? (0) | ? (1) | ? (0) | ? (1) |
| CAT | H | H (58) | H (35) | H (21) | H (51) | H (31) | H (53) | H (38) | H (49) |
| CAC | H | ? (0) | H (11) | H (30) | H (3) | H (16) | H (3) | H (11) | H (4) |
| CAA | Q | Q (38) | Q (34) | Q (26) | Q (40) | Q (36) | Q (43) | Q (33) | Q (43) |
| CAG | Q | ? (1) | ? (0) | Q (20) | ? (3) | ? (1) | ? (2) | ? (1) | ? (0) |
| CGT | R | R (12) | R (8) | R (36) | R (28) | R (9) | R (29) | R (8) | R (11) |
| CGC | R | ? (0) | ? (0) | ? (0) | ? (0) | ? (0) | ? (0) | ? (0) | ? (0) |
| CGA | R | ? (0) | ? (0) | ? (0) | ? (0) | ? (0) | ? (0) | ? (0) | ? (0) |
| CGG | R | ? (0) | ? (0) | ? (1) | ? (0) | ? (0) | ? (0) | ? (0) | ? (0) |
| ATT | I | I (121) | I (90) | I (69) | I (106) | I (122) | I (105) | I (86) | I (161) |
| ATC | I | ? (10) | I (21) | I (100) | ? (8) | I (46) | ? (8) | I (23) | I (35) |
| ATA | I | I (205) | I (193) | I (100) | I (203) | I (115) | I (231) | I (177) | I (127) |
| ATG | M | M (72) | M (57) | M (69) | M (61) | M (62) | M (69) | M (58) | M (62) |
| ACT | T | T (41) | T (37) | T (64) | T (41) | T (27) | T (45) | T (23) | T (63) |
| ACC | T | ? (0) | ? (1) | ? (8) | ? (0) | ? (0) | ? (4) | ? (0) | ? (0) |
| ACA | T | T (60) | T (74) | T (59) | T (76) | T (90) | T (79) | T (87) | T (60) |
| ACG | T | ? (0) | ? (1) | ? (5) | ? (4) | ? (0) | ? (1) | ? (0) | ? (0) |
| AAT | N | N (87) | N (65) | N (41) | N (71) | N (55) | N (81) | N (57) | N (83) |
| AAC | N | ? (4) | N (12) | N (50) | ? (4) | N (17) | ? (1) | N (17) | ? (4) |
| AAA | K | K (39) | K (41) | K (16) | K (38) | K (37) | K (41) | K (39) | K (49) |
| AAG | K | ? (0) | ? (2) | K (31) | ? (1) | ? (2) | ? (0) | ? (0) | ? (0) |
| AGT | S | S (61) | S (72) | ? (64) | S (100) | S (55) | S (102) | S (64) | S (81) |
| AGC | S | ? (0) | ? (2) | S (15) | ? (2) | ? (6) | ? (3) | ? (2) | ? (0) |
| AGA | R | R (39) | R (32) | R (19) | R (19) | R (36) | R (22) | R (34) | R (36) |
| AGG | R | ? (0) | ? (0) | ? (0) | ? (0) | ? (0) | ? (0) | ? (0) | ? (0) |
| GTT | V | V (49) | V (39) | V (87) | V (74) | V (54) | V (63) | V (45) | V (84) |
| GTC | V | ? (0) | ? (1) | ? (22) | ? (1) | ? (0) | ? (1) | ? (0) | ? (1) |
| GTA | V | V (89) | V (91) | V (74) | V (75) | V (86) | V (92) | V (91) | V (68) |
| GTG | V | ? (1) | ? (2) | V (24) | ? (3) | ? (3) | ? (2) | ? (3) | ? (2) |
| GCT | A | A (103) | A (58) | A (104) | A (86) | A (71) | A (92) | A (89) | A (97) |
| GCC | A | ? (3) | ? (9) | A (31) | A (8) | ? (3) | ? (6) | ? (7) | A (9) |
| GCA | A | A (59) | A (69) | A (42) | A (56) | A (86) | A (62) | A (57) | A (45) |
| GCG | A | ? (0) | ? (3) | A (9) | ? (5) | ? (0) | ? (4) | ? (0) | ? (1) |
| GAT | D | D (50) | D (37) | D (26) | D (44) | D (40) | D (42) | D (44) | D (38) |
| GAC | D | ? (0) | ? (1) | D (23) | ? (1) | ? (3) | ? (1) | ? (1) | ? (4) |
| GAA | E | E (47) | E (48) | E (24) | E (53) | E (43) | E (54) | E (43) | E (52) |
| GAG | E | ? (3) | E (6) | E (36) | ? (1) | E (11) | ? (1) | E (10) | ? (1) |
| GGT | G | G (91) | G (116) | G (155) | G (112) | G (126) | G (125) | G (103) | G (116) |
| GGC | G | ? (0) | ? (0) | ? (3) | ? (2) | ? (0) | ? (1) | ? (0) | ? (1) |
| GGA | G | G (65) | G (53) | G (15) | G (43) | G (50) | G (50) | G (64) | G (56) |
| GGG | G | ? (0) | ? (3) | G (14) | G (8) | ? (0) | ? (2) | ? (3) | ? (2) |

## 6. Method 2 per-codon consensus at conserved columns (all 64 codons)

Each cell gives `usage / N conserved columns: top residues`. A codon is called when N ≥ 5 and its top residue holds ≥ 80% of the columns. A **bold** cell is a codon whose call differs from table 4.

| codon | t4 | BK074526 | BK074532 | BK074539 | BK074262 | BK074538 | BK074252 | BK074542 | BK074222 |
|---|---|---|---|---|---|---|---|---|---|
| TTT | F | 218 / 120: F99, Y7, L6 ... | 202 / 100: F88, Y6, L3 ... | 111 / 56: F46, Y5, L4 ... | 218 / 126: F119, Y3, L3 ... | 150 / 64: F57, L4, Y3 | 216 / 124: F114, Y5, L4 ... | 186 / 84: F75, Y4, L3 ... | 183 / 81: F70, L7, Y2 ... |
| TTC | F | 41 / 27: F25, W1, I1 | 109 / 56: F54, Y2 | 187 / 104: F96, L4, Y3 ... | 42 / 26: F25, Y1 | 162 / 97: F90, Y4, L2 ... | 48 / 33: F33 | 127 / 73: F71, Y2 | 132 / 70: F67, Y3 |
| TTA | L | 548 / 254: L213, F9, V8 ... | 602 / 267: L249, I5, V3 ... | 199 / 80: L70, I6, F2 ... | 603 / 274: L261, I7, Y2 ... | 614 / 270: L254, I5, V3 ... | 636 / 273: L260, I7, V2 ... | 605 / 265: L251, I5, Y3 ... | 579 / 247: L228, V5, F3 ... |
| TTG | L | 1 / 0 | 4 / 2: L2 | 29 / 15: L13, V1, Y1 | 7 / 5: L5 | 8 / 1: L1 | 7 / 4: L3, I1 | 6 / 4: L3, V1 | 4 / 2: L1, G1 |
| TCT | S | 61 / 30: S19, G3, T3 ... | 59 / 24: S23, A1 | 127 / 51: S43, A3, P2 ... | 78 / 31: S29, A1, T1 | 39 / 15: S12, A3 | 61 / 29: S26, A1, C1 ... | 43 / 17: S16, G1 | 96 / 40: S40 |
| TCC | S | 0 / 0 | 4 / 2: S2 | 38 / 16: S13, K1, L1 ... | 7 / 4: S3, C1 | 1 / 0 | 0 / 0 | 1 / 1: S1 | 1 / 0 |
| TCA | S | 134 / 66: S52, G6, T3 ... | 131 / 60: S50, A3, P2 ... | 60 / 25: S20, A4, H1 | 116 / 51: S47, T1, P1 ... | 176 / 76: S67, P2, R1 ... | 136 / 58: S52, T3, A2 ... | 161 / 72: S61, A3, P2 ... | 98 / 38: S31, P2, A2 ... |
| TCG | S | 0 / 0 | 1 / 0 | 12 / 7: S6, T1 | 4 / 3: S3 | 0 / 0 | 3 / 1: C1 | 1 / 0 | 0 / 0 |
| TAT | Y | 227 / 123: Y104, F11, W3 ... | 237 / 123: Y114, F4, N2 ... | 135 / 65: Y60, N2, F1 ... | 244 / 130: Y122, F5, L2 ... | 226 / 119: Y110, N3, F3 ... | 254 / 129: Y124, F4, L1 | 222 / 114: Y108, F2, N2 ... | 253 / 129: Y116, F4, L3 ... |
| TAC | Y | 11 / 8: Y8 | 17 / 11: Y10, F1 | 128 / 70: Y64, F3, K1 ... | 15 / 8: Y8 | 31 / 16: Y16 | 9 / 4: Y4 | 29 / 19: Y18, F1 | 11 / 9: Y6, R1, I1 ... |
| TAA | * | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 |
| TAG | * | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 |
| TGT | C | 27 / 10: C6, T1, I1 ... | 32 / 15: C13, A1, F1 | 39 / 15: C9, M2, F1 ... | 36 / 14: C12, V1, T1 | 33 / 16: C14, A1, F1 | 37 / 15: C12, V2, T1 | 32 / 15: C13, A1, F1 | 36 / 16: C13, A1, F1 ... |
| TGC | C | 1 / 1: C1 | 0 / 0 | 5 / 3: C2, T1 | 0 / 0 | 2 / 0 | 0 / 0 | 0 / 0 | 0 / 0 |
| TGA | W | 54 / 39: W37, Y1, T1 | 60 / 41: W41 | 49 / 33: W33 | 56 / 43: W42, Y1 | 64 / 45: W44, Y1 | 59 / 42: W41, Y1 | 62 / 43: W43 | 63 / 43: W40, F2, L1 |
| TGG | W | 1 / 0 | 0 / 0 | 17 / 12: W11, V1 | 4 / 1: W1 | 0 / 0 | 0 / 0 | 1 / 1: W1 | 1 / 1: W1 |
| CTT | L | 7 / 3: L3 | **29 / 24: P24** | 144 / 64: L59, I2, T1 ... | **18 / 15: P15** | **18 / 17: P17** | 5 / 3: L2, I1 | **26 / 20: P19, Y1** | 22 / 9: L8, R1 |
| CTC | L | 0 / 0 | 0 / 0 | 30 / 12: L12 | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 |
| CTA | L | 4 / 1: L1 | **104 / 79: P78, Y1** | 194 / 76: L69, F1, M1 ... | **52 / 41: P41** | **120 / 92: P91, Y1** | 7 / 4: L4 | **110 / 88: P88** | 6 / 4: L3, F1 |
| CTG | L | 1 / 0 | 1 / 1: T1 | 86 / 36: L31, Y1, N1 ... | 1 / 0 | 1 / 0 | 0 / 0 | 1 / 0 | 1 / 1: L1 |
| CCT | P | 81 / 65: P64, S1 | 1 / 1: I1 | 72 / 55: P54, Q1 | 20 / 13: P13 | 0 / 0 | 53 / 40: P40 | 0 / 0 | 70 / 54: P54 |
| CCC | P | 0 / 0 | 0 / 0 | 8 / 6: P6 | 2 / 2: P2 | 0 / 0 | 1 / 1: P1 | 0 / 0 | 0 / 0 |
| CCA | P | 51 / 39: P38, V1 | 2 / 1: P1 | 55 / 45: P45 | 43 / 32: P32 | 0 / 0 | 86 / 64: P64 | 2 / 1: P1 | 65 / 50: P49, Y1 |
| CCG | P | 0 / 0 | 0 / 0 | 5 / 2: P2 | 2 / 2: P2 | 0 / 0 | 2 / 1: P1 | 0 / 0 | 1 / 1: P1 |
| CAT | H | 71 / 53: H47, G1, A1 ... | 53 / 41: H38, T1, R1 ... | 29 / 19: H19 | 69 / 51: H49, N1, F1 | 46 / 36: H34, R1, F1 | 70 / 53: H51, N1, F1 | 56 / 44: H42, R1, F1 | 63 / 47: H43, Y1, N1 ... |
| CAC | H | 0 / 0 | 12 / 12: H12 | 40 / 28: H28 | 5 / 4: H4 | 19 / 17: H17 | 3 / 3: H3 | 10 / 9: H9 | 4 / 3: H3 |
| CAA | Q | 62 / 42: Q37, N2, K1 ... | 61 / 36: Q34, D1, H1 | 45 / 23: Q23 | 69 / 40: Q38, G1, N1 | 66 / 37: Q35, D1, H1 | 74 / 41: Q38, G1, S1 ... | 64 / 38: Q36, D1, H1 | 63 / 37: Q35, S1, E1 |
| CAG | Q | 1 / 1: Q1 | 3 / 1: Q1 | 30 / 14: Q12, Y1, N1 | 4 / 1: Q1 | 2 / 2: Q2 | 1 / 0 | 1 / 1: Q1 | 1 / 0 |
| CGT | R | 17 / 10: R10 | 14 / 12: R12 | 54 / 29: R28, K1 | 37 / 29: R27, K1, N1 | 14 / 12: R12 | 34 / 26: R25, N1 | 13 / 12: R12 | 14 / 12: R12 |
| CGC | R | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 | 1 / 1: L1 |
| CGA | R | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 |
| CGG | R | 0 / 0 | 0 / 0 | 2 / 0 | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 |
| ATT | I | 222 / 58: I36, L11, V6 ... | 179 / 48: I31, L8, V3 ... | 100 / 34: I27, V5, L1 ... | 207 / 36: I30, L2, F2 ... | 204 / 62: I44, L9, V5 ... | 188 / 45: I36, V4, L3 ... | 160 / 55: I42, L7, V4 ... | 251 / 83: I63, L8, V8 ... |
| ATC | I | 16 / 7: I4, M1, T1 ... | 32 / 3: I3 | 141 / 39: I28, L5, V5 ... | 16 / 8: I6, V1, L1 | 84 / 23: I17, L5, V1 | 9 / 2: I2 | 33 / 9: I8, F1 | 50 / 11: I5, V3, L2 ... |
| ATA | I | 362 / 93: I40, L18, V14 ... | 341 / 86: I51, L18, V8 ... | 147 / 39: I27, V5, L5 ... | 371 / 74: I51, V10, L8 ... | 245 / 52: I26, L11, V9 ... | 397 / 77: I45, V11, L10 ... | 353 / 71: I36, L18, V11 ... | 265 / 44: I21, L12, V4 ... |
| ATG | M | 85 / 46: M32, L7, I3 ... | 93 / 42: M40, F1, V1 | 103 / 48: M39, L6, V2 ... | 95 / 45: M40, L2, C1 ... | 97 / 43: M42, I1 | 92 / 45: M41, I2, L2 | 93 / 46: M43, F2, L1 | 93 / 47: M39, L2, I2 ... |
| ACT | T | 69 / 31: T22, S4, V1 ... | 66 / 27: T25, N1, D1 | 93 / 49: T40, L4, V2 ... | 89 / 42: T40, P1, Y1 | 51 / 21: T21 | 84 / 40: T39, S1 | 50 / 21: T21 | 101 / 49: T41, S4, L2 ... |
| ACC | T | 0 / 0 | 2 / 1: T1 | 17 / 5: T4, S1 | 2 / 0 | 0 / 0 | 5 / 3: T3 | 0 / 0 | 0 / 0 |
| ACA | T | 91 / 46: T37, S4, A3 ... | 140 / 58: T51, S2, G1 ... | 99 / 41: T35, L2, G1 ... | 133 / 45: T43, A1, N1 | 158 / 73: T63, S3, V2 ... | 133 / 43: T41, N2 | 162 / 70: T63, S3, N2 ... | 102 / 39: T38, S1 |
| ACG | T | 0 / 0 | 1 / 1: T1 | 10 / 4: T3, D1 | 9 / 3: T3 | 0 / 0 | 0 / 0 | 1 / 1: I1 | 1 / 0 |
| AAT | N | 239 / 89: N66, S4, T3 ... | 216 / 58: N51, Y2, S2 ... | 87 / 33: N29, H1, R1 ... | 203 / 65: N61, Y2, T1 ... | 199 / 55: N51, H1, Y1 ... | 223 / 73: N66, K2, Y2 ... | 205 / 58: N53, F1, H1 ... | 224 / 78: N66, S5, G2 ... |
| AAC | N | 9 / 4: N2, D1, S1 | 24 / 13: N12, T1 | 101 / 39: N34, D1, H1 ... | 9 / 6: N6 | 34 / 15: N15 | 5 / 2: N2 | 31 / 12: N12 | 10 / 4: N4 |
| AAA | K | 96 / 44: K38, R3, L1 ... | 98 / 44: K37, G2, N2 ... | 29 / 15: K12, R1, N1 ... | 87 / 43: K40, R1, S1 ... | 92 / 45: K41, N2, I1 ... | 85 / 41: K39, R1, N1 | 94 / 44: K41, N2, I1 | 111 / 46: K41, Y1, S1 ... |
| AAG | K | 2 / 0 | 3 / 3: K2, A1 | 53 / 29: K27, R1, N1 | 3 / 1: K1 | 7 / 1: K1 | 2 / 1: K1 | 0 / 0 | 0 / 0 |
| AGT | S | 88 / 46: S32, G5, A4 ... | 130 / 49: S44, G3, A1 ... | 100 / 44: S37, G3, A2 ... | 168 / 61: S56, A2, G2 ... | 117 / 49: S42, G4, N1 ... | 169 / 63: S58, A2, G2 ... | 132 / 52: S46, G3, N1 ... | 120 / 55: S50, T2, R1 ... |
| AGC | S | 0 / 0 | 7 / 5: S4, L1 | 35 / 10: S6, G2, T1 ... | 2 / 1: S1 | 16 / 6: S6 | 2 / 2: S2 | 6 / 3: S3 | 2 / 1: S1 |
| AGA | R | 50 / 38: R34, K3, T1 | 51 / 36: R34, L1, K1 | 32 / 19: R18, N1 | 33 / 22: R21, K1 | 51 / 36: R35, K1 | 34 / 20: R19, K1 | 50 / 36: R35, K1 | 52 / 35: R32, L2, K1 |
| AGG | R | 0 / 0 | 0 / 0 | 0 / 0 | 3 / 0 | 0 / 0 | 1 / 0 | 0 / 0 | 1 / 1: P1 |
| GTT | V | 65 / 29: V20, I5, L2 ... | 58 / 23: V18, I2, A1 ... | 128 / 40: V33, I4, T1 ... | 106 / 44: V40, I2, A1 ... | 82 / 35: V30, I3, A1 ... | 96 / 38: V35, I2, A1 | 62 / 28: V22, T2, I2 ... | 103 / 49: V44, I2, G1 ... |
| GTC | V | 0 / 0 | 2 / 1: I1 | 33 / 5: V3, I1, S1 | 2 / 0 | 1 / 0 | 0 / 0 | 0 / 0 | 1 / 0 |
| GTA | V | 116 / 48: V40, L4, I1 ... | 154 / 65: V52, T4, I4 ... | 131 / 37: V24, L8, I2 ... | 128 / 49: V39, I4, L3 ... | 162 / 54: V42, I4, T3 ... | 136 / 48: V37, I6, F2 ... | 164 / 59: V45, I8, L3 ... | 89 / 32: V25, I2, P1 ... |
| GTG | V | 1 / 0 | 4 / 1: V1 | 42 / 8: V7, L1 | 5 / 1: V1 | 4 / 2: V2 | 2 / 1: V1 | 5 / 4: V4 | 3 / 1: V1 |
| GCT | A | 120 / 64: A44, G8, S6 ... | 84 / 40: A33, S3, L1 ... | 144 / 70: A49, S9, F2 ... | 117 / 55: A50, S3, P1 ... | 91 / 45: A43, S1, L1 | 116 / 61: A57, G2, S1 ... | 121 / 65: A53, S5, L1 ... | 124 / 62: A56, S2, K1 ... |
| GCC | A | 3 / 2: A2 | 11 / 7: A6, S1 | 47 / 20: A12, G3, L2 ... | 9 / 6: A6 | 6 / 3: A3 | 6 / 4: A4 | 10 / 5: A5 | 12 / 9: A9 |
| GCA | A | 71 / 37: A28, S2, T2 ... | 101 / 50: A42, S3, G1 ... | 51 / 18: A14, I2, D1 ... | 81 / 35: A34, P1 | 114 / 56: A42, S7, G2 ... | 81 / 33: A29, P1, V1 ... | 75 / 35: A30, S2, G1 ... | 68 / 25: A19, G3, V1 ... |
| GCG | A | 0 / 0 | 6 / 2: A2 | 13 / 6: A6 | 5 / 3: A3 | 0 / 0 | 4 / 1: A1 | 0 / 0 | 1 / 1: A1 |
| GAT | D | 84 / 55: D49, E4, N1 ... | 76 / 49: D47, N2 | 48 / 30: D28, N1, E1 | 86 / 49: D47, E1, Q1 | 77 / 48: D46, N2 | 96 / 49: D47, E1, N1 | 80 / 51: D49, N2 | 82 / 45: D45 |
| GAC | D | 0 / 0 | 4 / 1: D1 | 44 / 22: D20, N1, H1 | 4 / 2: D2 | 6 / 3: D3 | 2 / 1: D1 | 4 / 0 | 6 / 4: D4 |
| GAA | E | 73 / 51: E50, K1 | 73 / 54: E52, Q1, D1 | 41 / 25: E23, H1, D1 | 85 / 57: E55, D1, L1 | 69 / 46: E44, Q1, D1 | 87 / 59: E55, D3, L1 | 68 / 49: E47, Q1, D1 | 91 / 57: E54, D2, Q1 |
| GAG | E | 5 / 3: E3 | 15 / 7: E6, D1 | 56 / 35: E34, Q1 | 10 / 3: E2, D1 | 25 / 14: E14 | 3 / 1: E1 | 21 / 11: E11 | 3 / 1: E1 |
| GGT | G | 148 / 107: G95, A6, S3 ... | 174 / 126: G120, A3, S2 ... | 208 / 160: G151, A5, S2 ... | 180 / 131: G130, A1 | 186 / 138: G131, A3, S2 ... | 184 / 129: G128, A1 | 161 / 120: G116, A3, S1 | 164 / 132: G122, A4, S3 ... |
| GGC | G | 0 / 0 | 0 / 0 | 4 / 2: G2 | 2 / 1: G1 | 0 / 0 | 0 / 0 | 0 / 0 | 1 / 0 |
| GGA | G | 87 / 67: G64, S2, F1 | 70 / 58: G55, S2, A1 | 21 / 18: G18 | 68 / 48: G48 | 78 / 57: G56, A1 | 68 / 53: G53 | 99 / 73: G69, S2, V1 ... | 79 / 55: G53, S1, A1 |
| GGG | G | 0 / 0 | 8 / 4: G4 | 17 / 11: G11 | 10 / 8: G8 | 1 / 0 | 2 / 2: G2 | 4 / 3: G3 | 5 / 4: G4 |

## 7. Method 2 sensitivity: CTT and CTA under other settings

Each cell gives `Pro / Leu / other (N)` at conserved columns. The settings are: primary (masked panel protein in MAFFT, ≥5 relatives, ≥80% identity); `t4` (panel protein aligned unmasked, as table 4); `hsp-agree` (only codons where the MAFFT column equals the tblastn HSP column); and alternative conservation thresholds.

| genome | codon | primary | t4 panel | hsp-agree | n≥7, ≥90% | n≥5, ≥60% | n≥3, ≥80% |
|---|---|---|---|---|---|---|---|
| BK074526.1 | CTT | 0/3/0 (3) | 0/3/0 (3) | 0/3/0 (3) | 0/3/0 (3) | 0/3/2 (5) | 0/3/0 (3) |
| BK074526.1 | CTA | 0/1/0 (1) | 0/1/0 (1) | 0/1/0 (1) | 0/0/0 (0) | 0/1/1 (2) | 0/1/0 (1) |
| BK074532.1 | CTT | 24/0/0 (24) | 24/0/0 (24) | 24/0/0 (24) | 21/0/0 (21) | 24/0/0 (24) | 24/0/0 (24) |
| BK074532.1 | CTA | 78/0/1 (79) | 78/0/1 (79) | 78/0/1 (79) | 72/0/0 (72) | 88/0/4 (92) | 78/0/1 (79) |
| BK074539.1 | CTT | 0/59/5 (64) | 0/59/5 (64) | 0/58/5 (63) | 0/39/1 (40) | 0/75/16 (91) | 0/59/5 (64) |
| BK074539.1 | CTA | 0/69/7 (76) | 0/69/8 (77) | 0/69/7 (76) | 0/44/3 (47) | 0/91/25 (116) | 0/69/7 (76) |
| BK074262.1 | CTT | 15/0/0 (15) | 15/0/0 (15) | 15/0/0 (15) | 15/0/0 (15) | 16/0/1 (17) | 15/0/0 (15) |
| BK074262.1 | CTA | 41/0/0 (41) | 41/0/0 (41) | 41/0/0 (41) | 34/0/0 (34) | 44/0/2 (46) | 41/0/0 (41) |
| BK074538.1 | CTT | 17/0/0 (17) | 17/0/0 (17) | 17/0/0 (17) | 14/0/0 (14) | 17/0/0 (17) | 17/0/0 (17) |
| BK074538.1 | CTA | 91/0/1 (92) | 91/0/1 (92) | 91/0/1 (92) | 84/0/0 (84) | 101/0/3 (104) | 91/0/1 (92) |
| BK074252.1 | CTT | 0/2/1 (3) | 0/2/1 (3) | 0/2/1 (3) | 0/0/0 (0) | 0/2/1 (3) | 0/2/1 (3) |
| BK074252.1 | CTA | 0/4/0 (4) | 0/4/0 (4) | 0/4/0 (4) | 0/1/0 (1) | 0/5/0 (5) | 0/4/0 (4) |
| BK074542.1 | CTT | 19/0/1 (20) | 19/0/1 (20) | 19/0/1 (20) | 16/0/0 (16) | 20/0/1 (21) | 19/0/1 (20) |
| BK074542.1 | CTA | 88/0/0 (88) | 88/0/0 (88) | 88/0/0 (88) | 81/0/0 (81) | 96/0/3 (99) | 88/0/0 (88) |
| BK074222.1 | CTT | 0/8/1 (9) | 0/8/1 (9) | 0/8/1 (9) | 0/5/1 (6) | 0/13/1 (14) | 0/8/1 (9) |
| BK074222.1 | CTA | 0/3/1 (4) | 0/3/1 (4) | 0/3/1 (4) | 0/2/1 (3) | 0/3/1 (4) | 0/3/1 (4) |

## 8. Codon usage per genome (all 64 codons, modelled genes only)

These are counts over the codons in the gene models (core genes cox1-3, cob, atp6/8/9, nad1-6, nad4L; intron ORFs excluded; segment junction codons trimmed).

| codon | t4 | BK074526 | BK074532 | BK074539 | BK074262 | BK074538 | BK074252 | BK074542 | BK074222 |
|---|---|---|---|---|---|---|---|---|---|
| TTT | F | 218 | 202 | 111 | 218 | 150 | 216 | 186 | 183 |
| TTC | F | 41 | 109 | 187 | 42 | 162 | 48 | 127 | 132 |
| TTA | L | 548 | 602 | 199 | 603 | 614 | 636 | 605 | 579 |
| TTG | L | 1 | 4 | 29 | 7 | 8 | 7 | 6 | 4 |
| TCT | S | 61 | 59 | 127 | 78 | 39 | 61 | 43 | 96 |
| TCC | S | 0 | 4 | 38 | 7 | 1 | 0 | 1 | 1 |
| TCA | S | 134 | 131 | 60 | 116 | 176 | 136 | 161 | 98 |
| TCG | S | 0 | 1 | 12 | 4 | 0 | 3 | 1 | 0 |
| TAT | Y | 227 | 237 | 135 | 244 | 226 | 254 | 222 | 253 |
| TAC | Y | 11 | 17 | 128 | 15 | 31 | 9 | 29 | 11 |
| TAA | * | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| TAG | * | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| TGT | C | 27 | 32 | 39 | 36 | 33 | 37 | 32 | 36 |
| TGC | C | 1 | 0 | 5 | 0 | 2 | 0 | 0 | 0 |
| TGA | W | 54 | 60 | 49 | 56 | 64 | 59 | 62 | 63 |
| TGG | W | 1 | 0 | 17 | 4 | 0 | 0 | 1 | 1 |
| CTT | L | 7 | 29 | 144 | 18 | 18 | 5 | 26 | 22 |
| CTC | L | 0 | 0 | 30 | 0 | 0 | 0 | 0 | 0 |
| CTA | L | 4 | 104 | 194 | 52 | 120 | 7 | 110 | 6 |
| CTG | L | 1 | 1 | 86 | 1 | 1 | 0 | 1 | 1 |
| CCT | P | 81 | 1 | 72 | 20 | 0 | 53 | 0 | 70 |
| CCC | P | 0 | 0 | 8 | 2 | 0 | 1 | 0 | 0 |
| CCA | P | 51 | 2 | 55 | 43 | 0 | 86 | 2 | 65 |
| CCG | P | 0 | 0 | 5 | 2 | 0 | 2 | 0 | 1 |
| CAT | H | 71 | 53 | 29 | 69 | 46 | 70 | 56 | 63 |
| CAC | H | 0 | 12 | 40 | 5 | 19 | 3 | 10 | 4 |
| CAA | Q | 62 | 61 | 45 | 69 | 66 | 74 | 64 | 63 |
| CAG | Q | 1 | 3 | 30 | 4 | 2 | 1 | 1 | 1 |
| CGT | R | 17 | 14 | 54 | 37 | 14 | 34 | 13 | 14 |
| CGC | R | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 |
| CGA | R | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| CGG | R | 0 | 0 | 2 | 0 | 0 | 0 | 0 | 0 |
| ATT | I | 222 | 179 | 100 | 207 | 204 | 188 | 160 | 251 |
| ATC | I | 16 | 32 | 141 | 16 | 84 | 9 | 33 | 50 |
| ATA | I | 362 | 341 | 147 | 371 | 245 | 397 | 353 | 265 |
| ATG | M | 85 | 93 | 103 | 95 | 97 | 92 | 93 | 93 |
| ACT | T | 69 | 66 | 93 | 89 | 51 | 84 | 50 | 101 |
| ACC | T | 0 | 2 | 17 | 2 | 0 | 5 | 0 | 0 |
| ACA | T | 91 | 140 | 99 | 133 | 158 | 133 | 162 | 102 |
| ACG | T | 0 | 1 | 10 | 9 | 0 | 0 | 1 | 1 |
| AAT | N | 239 | 216 | 87 | 203 | 199 | 223 | 205 | 224 |
| AAC | N | 9 | 24 | 101 | 9 | 34 | 5 | 31 | 10 |
| AAA | K | 96 | 98 | 29 | 87 | 92 | 85 | 94 | 111 |
| AAG | K | 2 | 3 | 53 | 3 | 7 | 2 | 0 | 0 |
| AGT | S | 88 | 130 | 100 | 168 | 117 | 169 | 132 | 120 |
| AGC | S | 0 | 7 | 35 | 2 | 16 | 2 | 6 | 2 |
| AGA | R | 50 | 51 | 32 | 33 | 51 | 34 | 50 | 52 |
| AGG | R | 0 | 0 | 0 | 3 | 0 | 1 | 0 | 1 |
| GTT | V | 65 | 58 | 128 | 106 | 82 | 96 | 62 | 103 |
| GTC | V | 0 | 2 | 33 | 2 | 1 | 0 | 0 | 1 |
| GTA | V | 116 | 154 | 131 | 128 | 162 | 136 | 164 | 89 |
| GTG | V | 1 | 4 | 42 | 5 | 4 | 2 | 5 | 3 |
| GCT | A | 120 | 84 | 144 | 117 | 91 | 116 | 121 | 124 |
| GCC | A | 3 | 11 | 47 | 9 | 6 | 6 | 10 | 12 |
| GCA | A | 71 | 101 | 51 | 81 | 114 | 81 | 75 | 68 |
| GCG | A | 0 | 6 | 13 | 5 | 0 | 4 | 0 | 1 |
| GAT | D | 84 | 76 | 48 | 86 | 77 | 96 | 80 | 82 |
| GAC | D | 0 | 4 | 44 | 4 | 6 | 2 | 4 | 6 |
| GAA | E | 73 | 73 | 41 | 85 | 69 | 87 | 68 | 91 |
| GAG | E | 5 | 15 | 56 | 10 | 25 | 3 | 21 | 3 |
| GGT | G | 148 | 174 | 208 | 180 | 186 | 184 | 161 | 164 |
| GGC | G | 0 | 0 | 4 | 2 | 0 | 0 | 0 | 1 |
| GGA | G | 87 | 70 | 21 | 68 | 78 | 68 | 99 | 79 |
| GGG | G | 0 | 8 | 17 | 10 | 1 | 2 | 4 | 5 |
| total |  | 3721 | 3961 | 4105 | 4080 | 4049 | 4114 | 4003 | 3983 |

Gene models (codons per gene):

| genome | cox1 | cox2 | cox3 | cob | atp6 | atp8 | atp9 | nad1 | nad2 | nad3 | nad4 | nad4L | nad5 | nad6 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| BK074526.1 | 540 | 236 | 267 | 380 | 252 | 48 | 76 | 343 | 254 | 95 | 468 | 65 | 624 | 73 |
| BK074532.1 | 528 | 243 | 272 | 370 | 244 | 44 | 76 | 337 | 450 | 101 | 473 | 79 | 650 | 94 |
| BK074539.1 | 528 | 247 | 280 | 384 | 254 | 48 | 76 | 342 | 495 | 125 | 416 | 85 | 644 | 181 |
| BK074262.1 | 523 | 247 | 289 | 383 | 248 | 48 | 76 | 331 | 475 | 112 | 482 | 88 | 620 | 158 |
| BK074538.1 | 557 | 243 | 272 | 384 | 253 | 48 | 76 | 341 | 450 | 100 | 482 | 79 | 627 | 137 |
| BK074252.1 | 527 | 235 | 280 | 380 | 248 | 48 | 76 | 311 | 504 | 148 | 477 | 88 | 602 | 190 |
| BK074542.1 | 556 | 243 | 272 | 384 | 253 | 44 | 70 | 336 | 449 | 101 | 480 | 79 | 642 | 94 |
| BK074222.1 | 459 | 236 | 273 | 381 | 254 | 48 | 76 | 341 | 480 | 111 | 490 | 72 | 641 | 121 |

## 9. Reference set (method 2)

These are curated relatives with GenBank CDS annotation. Core-gene CDSs were selected by gene name; intron ORFs, maturases, endonucleases and partial CDSs were dropped, and the longest copy was kept. They were translated with table 4 to build the MAFFT L-INS-i alignments. For the tally, each relative's CUN/AUA/CGA/CGC codons and stops were masked to X, so a relative's own code (table 3 vs 4, or a CUN reassignment) cannot leak into the consensus.

| accession | organism | genes | which |
|---|---|---|---|
| HG937690 | Blastobotrys adeninivorans | 14 | atp6,atp8,atp9,cob,cox1,cox2,cox3,nad1,nad2,nad3,nad4,nad4L,nad5,nad6 |
| NC_002659 | Yarrowia lipolytica | 14 | atp6,atp8,atp9,cob,cox1,cox2,cox3,nad1,nad2,nad3,nad4,nad4L,nad5,nad6 |
| NC_005972 | Candida zemplinina | 7 | atp6,atp8,atp9,cob,cox1,cox2,cox3 |
| NC_014613 | Candida salmanticensis | 14 | atp6,atp8,atp9,cob,cox1,cox2,cox3,nad1,nad2,nad3,nad4,nad4L,nad5,nad6 |
| NC_016116 | Candida galli | 14 | atp6,atp8,atp9,cob,cox1,cox2,cox3,nad1,nad2,nad3,nad4,nad4L,nad5,nad6 |
| NC_020043 | Magnusiomyces magnusii | 14 | atp6,atp8,atp9,cob,cox1,cox2,cox3,nad1,nad2,nad3,nad4,nad4L,nad5,nad6 |
| NC_022155 | Candida santjacobensis | 14 | atp6,atp8,atp9,cob,cox1,cox2,cox3,nad1,nad2,nad3,nad4,nad4L,nad5,nad6 |
| NC_024095 | Magnusiomyces capitatus | 14 | atp6,atp8,atp9,cob,cox1,cox2,cox3,nad1,nad2,nad3,nad4,nad4L,nad5,nad6 |
| NC_025768 | Galactomyces candidum | 13 | atp6,atp8,atp9,cob,cox1,cox2,cox3,nad1,nad2,nad3,nad4,nad5,nad6 |
| NC_060461 | Saprochaete fungicola | 14 | atp6,atp8,atp9,cob,cox1,cox2,cox3,nad1,nad2,nad3,nad4,nad4L,nad5,nad6 |
| PQ436977 | Entelexis sorbosivorans | 10 | atp6,atp9,cob,cox1,cox2,cox3,nad1,nad4,nad4L,nad5 |

## 10. Scripts, data and commands

All scripts are in `/private/tmp/claude-501/-Users-ryan-lab/c532536a-a44b-4f0d-84b9-39c5fcb71590/scratchpad/blind_cun/src/`. They run as `uv run --no-project [--with scipy] python -I <script> ...`, and every data path is passed as an argument. Working data live in `/private/tmp/claude-501/-Users-ryan-lab/c532536a-a44b-4f0d-84b9-39c5fcb71590/scratchpad/blind_cun/work/`.

| step | script / command | output |
|---|---|---|
| download panel and references | NCBI E-utilities `efetch` (rettype fasta / fasta_cds_na) and `esearch`/`esummary`, each into its own `mktemp -d blind_cun/dl.XXXX` directory, ≤3 requests/s | `dl.*/panel.fna`, `dl.*/refs_cds_na.fna`, `dl.*/*.json` |
| select reference core-gene CDSs | `prep_refs.py refs_cds_na.fna work/refs ACC,...` | `work/refs/ref_<gene>.{t4,mask}.faa`, `refs_selected.tsv` |
| reference MSAs | `build_ref_msa.py work/refs work/msa 4` (MAFFT L-INS-i) | `work/msa/ref_<gene>.msa` |
| panel gene models | `find_genes.py work/refs work/msa work/genomes/<acc>.fna work/genes` (tblastn -db_gencode 4 -seg no, masked queries) | `work/genes/<acc>.{genes.json,models.tsv,hsps.tsv}` |
| align panel to MSAs, record columns | `align_tally.py work/refs work/msa work/genes work/tally ACC,... 4` (mafft --add --localpair --maxiterate 1000; masked and t4 modes) | `work/tally/records.tsv`, `work/tally/aln/` |
| per-codon tallies | `tally_codons.py records.tsv ref_order.json outdir [min_nref] [min_frac]` (5 0.8 primary; 7 0.9, 5 0.6, 3 0.8 sensitivity) | `work/tally*/codon_table.tsv`, `conserved_records.tsv` |
| Codetta (method 1) | `HMMER_NCPU=4 uv run --no-project --python 3.12 --with 'numpy<2' --with 'scipy<1.14' python -E -s dl.bllj/codetta/codetta.py work/codetta/<acc>.fna -m`, then `codetta_infer.py <acc>.fna -m -r 0.99 --inference_output <acc>.relaxed_r0.99.genetic_code.out` | `work/codetta/<acc>.fna.Pfam-A_enone.hmm.1e-10_0.9999_0.01_excl-tvuy.genetic_code.out`, `*.relaxed_r0.99.genetic_code.out` |
| this report | `report.py work refsum_dir src report.md` (needs scipy) | `report.md` |

## 11. Caveats

- **Codetta setup differed from stock.** The bundled HMMER 3.3.2 does not configure on Apple-silicon arm64 (it needs SSE). Codetta therefore ran with Homebrew HMMER 3.4 (`hmmscan`/`hmmpress`/`hmmsearch`) symlinked into its expected `hmmer-3.3.2/bin`. The easel miniapps (`esl-sfetch`, `esl-translate`) were built from the bundled 3.3.2 easel source after running `autoconf`. Codetta needed numpy<2 and ran on Python 3.12. Pfam-A_enone was downloaded from eddylab.org/publications/Shulgina21/Pfam-A_enone.tar.gz. Codetta ran as cloned with `-m`, which includes the mitochondrial Pfams that it excludes by default; that suits mitogenomes. Setup took well under the 45-minute budget.
- **Codetta is underpowered for rare codons in mitogenomes.** Each genome has only about 2,200-2,650 aligned codon-columns. Rare codons (CTC, CTG, CCC/CCG, CGN, AGG and others) therefore come out '?', and the per-codon logP and N are given in section 5. Where CUN is rare (BK074526.1, BK074252.1, BK074222.1), Codetta leaves CUN as '?'. For CTT, logP(L) is near 0 in all three while logP(P) is -19 to -37. CTA has logP(P) of -11 to -19, except BK074526.1, where CTA aligned in only 1 column and is uninformative. The relaxed run calls CTT = L in BK074222.1 only.
- **Low CUN counts in three genomes.** BK074526.1 has 12 CUN codons and 4 scorable columns, BK074252.1 has 12 and 7, and BK074222.1 has 29 and 14. All scorable columns are Leu (12/14 in BK074222.1, with 0 Pro). This rules out a wholesale CUN→Pro reading in these genomes. It does not prove every CUN codon is Leu, and CTC/CTG are essentially unobserved there too.
- **CTC and CTG are uncalled in the four CUN→Pro genomes**, where each is used 0-1 times. The data show only that CTT and CTA are Pro. Whether CTC/CTG would also be Pro, or are unassigned, cannot be inferred from these sequences.
- **How Pro is encoded varies.** BK074538.1 uses no CCN codons at all in the modelled genes, BK074532.1 uses 3 and BK074542.1 uses 2, so CUN is effectively their only Pro codon family. BK074262.1 keeps CCN as well: its conserved Pro columns are 49 CCN and 56 CUN. Leu in all four is encoded almost only by TTR. Genome-wide CC/GG dinucleotides are not depleted in these genomes (checked separately), so a sequence-assembly artifact that turns CC into CT is unlikely.
- **Reference choice and masking.** There are 11 relatives with GenBank CDS annotation (section 9). The closest is 80-87% identical to the panel proteins in the example genes, and some lack genes (section 9). Relatives were translated with table 4 only to build the alignments. Their own CUN/AUA/CGA/CGC codons and stops were masked out of every consensus. A Pro consensus therefore comes only from the relatives' CCN codons, and the relatives' own code cannot create the signal. The panel protein was also masked at those codons before alignment (primary). An unmasked table-4 alignment ('t4 panel'), a check against tblastn HSP columns, and three conservation thresholds give the same calls (section 7).
- **HSP column cross-check.** MAFFT columns agree with the tblastn HSP columns for 30680/31709 (96.8%) of panel codons in masked mode. A key-mismatch bug at first left this cross-check empty. It was fixed, the alignments were reused unchanged, and the calls did not change.
- **Gene models are approximate.** They are chained tblastn HSPs, not annotated CDSs. Starts and stops are not modelled, intron ORFs are excluded, and a few codons at exon junctions are dropped where segment ranges overlap. All 32016 model codons were checked against the genome sequence at their reported coordinates and all match. Example flanking peptides follow the spliced model. For 4 of the 64 example rows, the ±5 window crosses an exon junction, so it differs at its edge from the contiguous genome translation; the focal codon is correct in all 64.
- **ATA.** At conserved columns the Ile share for ATA (pooled ~0.55) is lower than for ATT/ATC (~0.72). ATA sits at variable hydrophobic positions more often, where the consensus is L or V. Met appears at only 4 of 536 pooled columns, and Codetta gives logP(M) ≤ -89 in every genome. ATA is Ile (table 4), not Met (table 3).
- **Other low-consensus cells**, such as TCT in BK074526.1, have their top residue equal to the table-4 meaning and reflect variable columns. In section 6, no codon other than CTT/CTA has a non-table-4 top residue at 2 or more conserved columns.
- **No tRNA analysis.** aragorn, tRNAscan-SE and Infernal are not installed and were not installed for this run, so anticodon evidence (for example a tRNA-Pro(UAG)) was not examined. The inference rests on protein-level evidence only.
- **Blinding.** Organism names were visible in the panel FASTA headers (genus names Wickerhamiella, Sugiyamaella, Starmerella) and in the reference esummary records. No external information on the genetic codes of these taxa was consulted. Nothing under the excluded lab directories was read. Early on, one `ls` of the scratchpad's top level showed other files' names; none were opened.

