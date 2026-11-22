<!-- Copied from the lab record (work/silent-operon/, 22 Nov 2026). Paths such as out/, src/, recheck/ or frozen-core/ refer to the original working folder; the small tables from out/ are in this folder (see ../README.md), the rest is not included. -->
# A_processing: can the divergent operon mature, and how is it expressed? (22 Nov 2026)

Rules: `RULES.md` (frozen before results; A1 added after the control rows, A2 post hoc). Scripts: `src/`; tables: `out/`.
Run: `uv run python src/locate.py && uv run python src/operons.py`, then `src/t1_stems.py`, `src/t1b_fold.py 16S|23S`,
`src/t1c_helix.py`, `src/t2_upstream.py`, `src/t2_genes.py`, `src/t2_cons.py`, `src/t2_spacer_boxa.py`, `src/t3_content.py`.
Loci: 15 carriers, 6 genera; mature ends from full-CM cmsearch (RF00177/RF02541/RF00001) in `out/operons.tsv`.
The 8 Aquicultor_A genomes (*Ca.* A. secundus, the RNA species) carry identical operons, so they count once.
Controls: E. coli K-12 (7 operons) and *Thermobispora bispora* (4 operons; two ~6%-divergent types, both transcribed).

## T1 RNase III processing stems

**Finding.** In A. secundus and DASXVE01, the divergent 16S has no leader-trailer stem above background, while the
normal 16S in the same genomes has one. The divergent 23S keeps a short leader-trailer core (14-18 bp) at the mature
ends in all 6 genera, conserved by compensatory changes. Its overall pairing is weaker than the normal stem and is
only marginally above shuffles (duplex dG p 0.015-0.15 against 0.005 for normal).

**Calibration.** Neither 16S method passed the frozen bar on *T. bispora*. The 200x200 duplex gave p 0.075-0.19; the
whole-precursor fold called stems in 3 of 4 copies. So the 16S calls are descriptive. The 23S whole-precursor fold is
calibrated: stems in E. coli 7/7 and *T. bispora* 3/3.

| locus (one per distinct sequence) | 16S helix bp (p_helix; duplex dG p) | 23S helix bp (p_helix) | 23S stem in whole fold |
|---|---|---|---|
| E. coli (7) | 35 (0.005; 0.005) | 17-18 (0.10-0.13) | 7/7 |
| *T. bispora* (4 / 3) | 14-22 (0.17-0.85; 0.075-0.19) | 27-37 (0.005-0.01) | 3/3 |
| A. secundus normal | 23 (0.005-0.035; 0.005-0.01) | 31 (0.005-0.01) | yes (60 pairs) |
| A. secundus **div** | **13 (0.52-0.60; 0.89-0.91)** | 15 (0.33-0.45) | **no (0 pairs)**; 15-bp duplex core present |
| DASXVE01 normal / **div** | 21 (0.04; 0.005) / **11 (0.81; 0.36)** | 30 / 17 (0.005 / 0.10) | yes / yes |
| JASEIR01 normal / div | 29 (0.005) / n.a. | n.a. | n.a. |
| CADDYS01 normal / div (2 genomes) | 9 (0.98) / 13-20 (0.10-0.67) | normal n.a. (23S IVS: my 3' end was wrong, see Recheck) / 18-23 (0.03-0.13) | n.a. / yes |
| DASQQM01 div | 23 (0.025; dG 0.58) | 18 (0.13) | yes |
| Kmv33 div | 19 (0.09; dG 0.34) | 20 (0.045) | yes |

- **Divergent 23S stem.** Weaker than normal, but with a conserved core. A core of 14-18 bp pairs the leader block AATAGCGTGCAT(C/A/T)GA with the trailer block
  TC(G/T/A)ATGCACGCTATT. It sits 1-3 nt from both mature ends, the classic 23S RNase III geometry, in all 6 genera
  (`out/t1_stems.tsv`). Between genera the sequence changes while the pairing holds: A-U becomes C-G, and U-G becomes U-A.
  That is selection to keep the stem. It is about as long as E. coli's 23S stem (17-18 bp) and about half the normal
  copy's (30-31 bp). The A. secundus version, 15 bp with one G-U, is the weakest: the whole-precursor MFE does not
  close it.
- **Divergent 16S.** There is no long leader-trailer helix. Instead, right after the 16S 3' end, the trailer holds a
  GC-rich hairpin (CTCCGGGTTCCCCTGwGGGCAACCC...). The hairpin is in all 14 carrier genomes where that region is
  assembled, once per genome, and nowhere else in any genome. Its stem keeps pairing through a T-A to A-T swap
  (Aquicultor_A, DAHYRC vs CADDYS, DASQQM, Kmv33).
- **Mechanism.** The 23S probably can still be cut out, through a shorter, E. coli-length stem that is kept by
  selection. In A. secundus the 16S probably cannot be cut out by the normal RNase III route. That fits Deneb's RNA reads: divergent reads run from the 16S 3' end into the spacer, while normal
  reads stop at the 16S end. Non-accumulation of the 23S is not explained by a missing excision stem.

## T2 Promoter and antitermination

**Finding.** The divergent operon lacks the normal operon's conserved expression signals. It has no leader boxA. Its
promoter is consensus-weak in A. secundus and DASXVE01, and its upstream region is not conserved across genera.
Instead, it carries a strong boxA in the 16S-23S spacer, at a fixed position about 40 nt before the 23S, in 4 of 6
genera.

- **Normal operons.**
  - Promoter: TTGACA-17-TATAnT with an extended -10 TG, scoring 11/12 in every genus (shuffle p 0.005-0.025), 400-520
    nt upstream of the 16S. The discriminator is AT-rich (0-1 GC of 7), not E. coli's GC-rich one (6/7).
  - Leader boxA: AGCTCTTTAAGA (2 mm, p 0.03-0.065) in Aquicultor_A and CADDYS01; TGGTCTTTGAAA (3 mm) in DASXVE01,
    JASEIR01 and Kmv33.
  - Same-strand gene about 650 nt upstream.
  - The 300 nt upstream are conserved across 5 genera: 66% identity against 49% for shuffled pairs.
- **Divergent operons.**
  - Leader boxA: best 3-4 mm (p 0.18-0.89). The exception is the DAHYRC copy (2 mm, p 0.03), and its sister CADDYS copy
    has 3 mm.
  - Promoter: 9/12 (p 0.62-0.68) in A. secundus and DASXVE01; 11-12/12 in CADDYS01, DASQQM01, Kmv33 and JASEIR01, but
    not at a shared position.
  - The upstream gene is on the opposite strand in 4 genera, so the operon needs its own promoter and does not get
    read-through.
  - Upstream identity across 6 genera is 51%, against 48% for shuffled pairs (95th percentile 51%), so it is essentially
    unconserved.
  - Spacer boxA: TGCTCTTTAA(A/C)A at 0-2 mm (p 0.005) in A. secundus, CADDYS01, DASQQM01 and Kmv33; 3 mm at the same
    spot in DAHYRC; absent in DASXVE01. It sits inside the 23S leader stem. Among normal spacers only DASXVE01 has one.
- **Mechanism.** Selection keeps the divergent 23S's own antitermination and excision signals, but not its promoter or
  leader. That suggests weak or conditional transcription initiation, and a 16S stretch that is transcribed without an
  antitermination complex, which normally also chaperones rRNA folding.

## T3 Intergenic content

**Finding.** Neither operon type has spacer tRNAs; that is family-wide (as in *T. bispora*), unlike E. coli. The
divergent operon has lost its 5S in 3 of 5 assessable genera, A. secundus among them. In the other 2 the 5S is intact.

- **5S.** The normal operon has a complete 5S in all 5 genera with downstream sequence (score 100-104, 229-1,107 nt after
  the 23S). In the divergent operon:
  - DASQQM01 and Kmv33: complete 5S (score 95 / 91, which is 88-92% of normal; 88-89% identical to the normal 5S;
    31 nt after the 23S).
  - Aquicultor_A (all 8 genomes), CADDYS01 (both genomes) and DASXVE01: no 5S. There is no cmsearch hit at E <= 10 and
    no blastn remnant within 1,500 nt, all of it on the contig. A same-strand protein gene starts about 210-310 nt
    after the 23S.
- **Mechanism.** In A. secundus a divergent 50S would need the normal 5S in trans. Losing the cis 5S while keeping the
  23S excision stem fits a 23S whose role does not depend on making its own complete 50S unit.

## Overall mechanism reading
The operon keeps signals for transcribing and cutting out its 23S: the spacer boxA, and a short 23S stem core with
compensatory changes. Both are weaker than in the normal operon. It is not built to make mature 16S in A. secundus: there is no 16S
leader-trailer stem, and the 16S runs unprocessed into the spacer in the RNA. Its expression signals (promoter, leader
boxA) are not under the selection the normal operon's are, and in A. secundus it has no 5S. Non-accumulation of the
23S is therefore not explained by a broken excision stem. It points downstream of RNase III: turnover of a
cut-but-unassembled 23S (changed PTC, no cis 5S), or a precursor that is cut and degraded as such.

## Recheck (fresh in-session subagent; own coordinates and scripts in `recheck/`; it did not see out/ or this file)
- **5S: confirmed.** 5/7 divergent operons have none: Aquicultor_A x2, CADDYS, DAHYRC, DASXVE; the best cmsearch score
  was 4-9.5 even with `--max -T 0`, and blastn found nothing. DASQQM01 and Kmv33 carry a divergent-type 5S 31 nt after
  the 23S: the two are 96% identical to each other and 85-89% identical to the normal 5S. All normal operons have a 5S.
- **boxA: confirmed.** The divergent spacer boxA has 0-2 mm (P < 0.01) in A. secundus, Kmv33, DASQQM and CADDYS.
  Upstream of the divergent 16S there is none, except DAHYRC at 2 mm (P 0.036). Normal operons show the reverse.
- **Stems: same direction, and the recheck is more sceptical about the divergent 23S.**
  - Divergent 16S: background by all 3 statistics (0-1 of 7). Normal 16S: 3 of 4 significant by dG.
  - Divergent 23S: 2 of 7 clearly significant (DASQQM, DAHYRC), against 5 of 5 normal.
  - Calibration: the 23S passes by dG (10/10 controls). No 16S statistic passes on *T. bispora*.
  - The recheck also found that my CADDYS01 and JASEIR01 normal 23S 3' ends were wrong, because an intervening sequence
    misplaces the cmsearch end. With the corrected ends, the normal 23S stem is present in all 5 normal operons. My
    "CADDYS01 normal 23S no stem" was an artefact.
- **Net.** The divergent excision stems are weaker than normal: the 16S looks like background and the 23S is short or
  marginal. The lead "16S cannot be excised" is supported in direction but stays uncalibrated.
