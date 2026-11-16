# CGG re-check: 9 genomes, judged blind (written 2026-11-11 16:12 EST)

## Short answer

| label | verdict |
|---|---|
| G1 | **Not Arg. +1 frameshift at CGG is obligatory.** The 4 nt CGG+N are read as one **Gly** (the C is effectively skipped). |
| G2 | **Read as Arg.** Low confidence: only 26 sites. |
| G3 | **+1 frameshift at CGG is obligatory.** The residue is probably Gly, but there are too few sites to call it. |
| G4 | **Not Arg. +1 frameshift at CGG is obligatory, and Gly is incorporated.** |
| G5 | **Read as Arg.** |
| G6 | **Read as Arg.** |
| G7 | **Not Arg. +1 frameshift at CGG is obligatory, and Gly is incorporated.** |
| G8 | **Read as Arg.** The residue evidence is weaker than for G5 and G6. |
| G9 | **Not Arg. +1 frameshift at CGG is obligatory, and Gly is incorporated.** |

No genome shows CGG as a stop. No genome shows CGG reassigned while staying in frame. No genome shows a −1 frameshift at CGG.

## 1. Method

**Reference proteins.** I took 8 RefSeq isolates that use the standard code (table 11): Bsub GCF_000009045.1, Cdif GCF_000009205.2, Chon GCF_000981035.1, Fbut GCF_023347275.1, Aint GCF_034334985.1, Bhan GCF_002222595.2, Alac GCF_018118585.1 and Leli GCF_000146185.1. Genes were called with pyrodigal. Every residue encoded by CGG, and every stop, was written as X, so the references never vote on CGG. That leaves 25,806 proteins of at least 60 aa.

**Alignment.** All reference proteins were aligned to each genome with `miniprot -S` (frameshift-aware, no splicing). I did not use DIAMOND or `sprot_p50`.

**Per-codon statistics** (`03_analyze.py`):
- Best alignment per reference.
- In-frame locus: the codon sits in a codon-to-residue column and the ±7 aa identity is at least 0.3.
- Consensus residue: at least 2 references have a clean column and at least 60% of them agree.
- Frameshift near the codon: a frameshift within 9 nt, supported by at least half of the covering references.
- Junction count: the number of merged frameshift events whose old-frame 5-codon window contains the codon.
- Repeated copies are collapsed by their ±30 nt context.
- These statistics are computed the same way for CGG and for every null codon.

**Residue at a CGG followed by a frameshift** (`04_fs_detail.py`). For each such event I built a hypothesis translation: the old frame, then `X` for the CGG plus the next nucleotide (4 nt), then the new frame. I aligned it locally to the reference window (BLOSUM62) and read the reference residue facing `X`.

**Null for residue identity** (`11_bestref.py`). For each locus I took the residue of the single most similar reference, split into identity bins (below 0.6, 0.6 to 0.8, 0.8 and above).

**tRNAs.** I used the codon-atlas `trna.parquet` (`06_trna.py`). As a cross-check, I ran blastn with all 9 genomes' tRNAs against each genome (`10_trna_blast.py`).

**Controls.**
- RF2/prfB +1 frameshift (`05_rf2.py`).
- Planted edits in G6 (`07_plant.py`, `08_score_plant.py`).

**Scripts.** All are in `src/recheck/`:
- Pipeline: `01_refs.py`, `02_run_miniprot.sh`, `03_analyze.py`, `04_fs_detail.py`, `05_rf2.py`, `06_trna.py`, `07_plant.py`, `08_score_plant.py`, `09_table.py`, `10_trna_blast.py`, `11_bestref.py`.
- Helpers: `mp.py` (PAF parser), `show.py`, `pairs.py`, `ggn.py`, `contig_gc.py`, `same_frame_list.py`, `trna_all.py`.

**Outputs.** Everything is in `out/` and `trna/`. The main files are:
- `out/codon_stats.txt`
- `out/table09.txt`
- `out/bestref.txt`
- `out/fs_detail_summary.json`
- `out/*.fs_cgg_sites.tsv`
- `out/rf2_control.tsv`
- `out/plant_control.txt`
- `trna/summary.txt`

## 2. Per-genome table

Column definitions:
- **Aligned CGG**: unique CGG loci that sit in-frame at an aligned position, or in the old frame at a frameshift. The number in brackets is how many have a consensus residue.
- **Frame**: what happens after the CGG, as same / +1 / −1. Same means no frameshift within 9 nt.
- **Residue opposite**: consensus counts (all clean pairs in brackets). At shifted sites this is the residue facing the CGG+N unit. At same-frame sites it is the residue facing the CGG.
- **Events with CGG**: frameshift events with CGG in the old-frame window, out of all frameshift events.
- **Null codons**: in-frame loci for CGA / CGC / AGG, given as n: same/+1/−1 and the top consensus residue. Fuller detail is in table 2b.
- **tRNA-Arg**: anticodons, written 5'→3' as DNA.

| label | Aligned CGG | Frame (same/+1/−1) | Residue opposite | Events with CGG | Null codons | tRNA-Arg | Verdict |
|---|---|---|---|---|---|---|---|
| G1 | 181 (134) | 0 / 178 / 3 | CGG+N: G127 A3 K1 L1 N1 S1 (pairs G657 A36 S22 … R1, of 774) | 181 / 212 | CGA 91: 90/1/0 R38 · CGC 23: 23/0/0 R15 · AGG 590: 584/5/1 R263 K27 | ACG, CCT, TCT | +1 frameshift, Gly |
| G2 | 26 (2) | 26 / 0 / 0 | same: R2 (pairs R15 K3 Q3 L3 T2 G1) | 0 / 41 | CGA 327: 326/0/1 R70 · CGC 98: 97/1/0 R25 · AGG 1677: 1676/1/0 R452 K51 | ACG only (incomplete, see caveats) | Arg (low confidence) |
| G3 | 14 (1) | 1 / 13 / 0 | CGG+N: A1 (pairs G6 A3 E2 C1 R1 H1) | 14 / 18 | CGA 130: 130/0/0 R44 · CGC 44: 44/0/0 R17 · AGG 576: 574/2/0 R252 | ACG, CCT, TCT | +1 frameshift; residue undetermined (Gly-leaning) |
| G4 | 46 (26) | 6 / 40 / 0 | CGG+N: G25 S1 (pairs G121 S5 R5 A4 …, of 145) | 40 / 72 | CGA 117: 117/0/0 R29 · CGC 44: 43/1/0 R23 · AGG 464: 459/5/0 R184 | ACG, CCT, TCT | +1 frameshift, Gly |
| G5 | 308 (75) | 308 / 0 / 0 | same: R64 K4 E2 Q2 P1 M1 S1 (pairs R315 K72 Q29 E20 …, R 0.56) | 0 / 35 | CGA 994: 994/0/0 R253 · CGC 704: 703/1/0 R221 · AGG 2755: 2753/2/0 R827 | ACG, CCG (2 loci), CCT, TCT | Arg |
| G6 | 87 (37) | 87 / 0 / 0 | same: R31 K5 P1 (pairs R158 K44 E14 …, R 0.54) | 0 / 17 | CGA 331: 331/0/0 R122 · CGC 149: 149/0/0 R58 · AGG 1363: 1362/1/0 R538 | ACG, CCG, CCT, TCT | Arg |
| G7 | 299 (148) | 13 / 281 / 5 | CGG+N: G134 A3 E2 T2 Q2 K1 S1 (pairs G675 A40 S37 … R7, of 916); same: M1 K1 R1 | 288 / 352 | CGA 205: 196/8/1 R63 · CGC 89: 86/3/0 R25 · AGG 916: 897/18/1 R332 | ACG, CCT, TCT | +1 frameshift, Gly |
| G8 | 45 (5) | 44 / 1 / 0 | same: R4 K1 (pairs R30 K16 N9 S6 …, R 0.36) | 2 / 38 | CGA 252: 252/0/0 R56 · CGC 85: 85/0/0 R12 · AGG 1779: 1777/2/0 R582 | ACG, CCT, TCG, TCT×3 | Arg (moderate confidence) |
| G9 | 35 (16) | 6 / 28 / 1 | CGG+N: G14 A1 Y1 (pairs G61 A3 Y3 …, of 73) | 29 / 48 | CGA 311: 311/0/0 R81 · CGC 164: 164/0/0 R48 · AGG 509: 508/1/0 R158 | ACG, CCT, TCT | +1 frameshift, Gly |

### 2b. CGG compared with null codons on identical statistics

This uses only the `03` in-frame-locus statistics, computed the same way for every codon. Each cell gives: in-frame loci, frameshift within 9 nt as +1/−1, and junction events.

| label | CGG | CGA | CGC | AGG | AGA | GCG | CCG |
|---|---|---|---|---|---|---|---|
| G1 | 30, 29/1, 181 | 91, 1/0, 2 | 23, 0/0, 0 | 590, 5/1, 9 | 4596, 29/2, 38 | 228, 4/0, 3 | 115, 1/1, 0 |
| G2 | 26, 0/0, 0 | 327, 0/1, 0 | 98, 1/0, 2 | 1677, 1/0, 0 | 15805, 5/0, 3 | 1344, 1/0, 0 | 735, 0/1, 1 |
| G3 | 5, 4/0, 14 | 130, 0/0, 0 | 44, 0/0, 0 | 576, 2/0, 2 | 4780, 1/1, 1 | 496, 2/0, 0 | 252, 0/0, 0 |
| G4 | 16, 10/0, 40 | 117, 0/0, 0 | 44, 1/0, 0 | 464, 5/0, 5 | 5730, 12/2, 16 | 210, 0/0, 0 | 130, 2/1, 3 |
| G5 | 308, 0/0, 0 | 994, 0/0, 0 | 704, 1/0, 1 | 2755, 2/0, 4 | 10677, 3/1, 4 | 2111, 1/0, 1 | 1338, 1/2, 4 |
| G6 | 87, 0/0, 0 | 331, 0/0, 0 | 149, 0/0, 0 | 1363, 1/0, 2 | 4420, 1/3, 2 | 1345, 0/0, 0 | 1869, 0/0, 1 |
| G7 | 93, 76/5, 288 | 205, 8/1, 7 | 89, 3/0, 3 | 916, 18/1, 29 | 5275, 68/6, 74 | 1195, 21/1, 22 | 586, 8/0, 2 |
| G8 | 45, 1/0, 2 | 252, 0/0, 0 | 85, 0/0, 0 | 1779, 2/0, 1 | 10576, 1/7, 4 | 906, 0/0, 0 | 339, 0/1, 3 |
| G9 | 11, 4/1, 29 | 311, 0/0, 0 | 164, 0/0, 1 | 509, 1/0, 3 | 5824, 3/1, 6 | 1334, 0/0, 0 | 457, 0/1, 0 |

**Stop test.** I compared how often each codon is in frame versus out of frame at the same aligned loci:
- CGG is 0.009–0.11 in every genome (lowest in G3 and G9, where most CGG sites shift).
- CGA/CGC/AGG are 0.13–0.50.
- TAA/TAG/TGA are at most 0.0006.

So in G2, G5, G6 and G8 CGG is rare but not stop-like. In the frameshift genomes, the in-frame count is low because most CGG sites shift.

**Identity-stratified best-reference residue** (`out/bestref.txt`). This gives the R fraction in three identity bins: below 0.6, 0.6 to 0.8, and 0.8 and above.

| label | CGG | CGA | AGG |
|---|---|---|---|
| G5 | 0.38 / 0.64 / 0.73 (n 137/105/37) | 0.38 / 0.64 / 0.79 | 0.38 / 0.62 / 0.81 |
| G6 | 0.37 / 0.50 / 0.70 (n 35/40/10) | 0.38 / 0.70 / 0.84 | 0.45 / 0.64 / 0.85 |
| G2 | 0.17 / 0.38 / 1.00 (n 6/8/4) | 0.44 / 0.69 / 0.77 | 0.46 / 0.67 / 0.84 |
| G8 | 0.19 / 0.43 / 1.00 (n 21/14/1) | 0.42 / 0.71 / 0.82 | 0.45 / 0.63 / 0.83 |

- In G8, below 0.8 identity CGG faces R 10/35 and K 12/35, against R 0.53 and K 0.16 for AGG.
- A third (4 of 12) of G8's K-facing CGGs lie in one ~10 kb cluster (NZ_LT906477.1, 887.8–897.2 kb). Each has one or two references at identity 0.36–0.71.

## 3. Positive controls and null results

**prfB (RF2) programmed +1 frameshift.** This is the natural positive control.
- Full-length RF2 (B. subtilis P28367, C. perfringens Q0TNE1, E. coli P07012) was aligned with `miniprot -S`.
- The expected F1 (+1) was called at RF2 position 24, at the canonical `CTT|T|GAC` slip site, in all 8 genomes where RF2 aligns (G2–G9):
  - With the BACSU or CLOP1 query, in the right direction.
  - In G8, only by CLOP1.
  - In G2, CLOP1 placed it 4 codons upstream, while BACSU was exact.
- No spurious frameshift was called anywhere else in RF2.
- G1 has no RF2 alignment at all, and tblastn found no RF1/RF2 homolog. This is probably a gap in the MAG.

**Planted frameshifts in G6** (a genome with no natural CGG events). The edits went through the full production pipeline: miniprot, then 03/04, then scoring.

| planted | detected within 9 nt, right direction | wrong direction | residue facing the CGG+N unit |
|---|---|---|---|
| +1 (random base inserted) | 149/150 | 0 | – |
| −1 (base deleted) | 149/150 | 0 | – |
| cgg_gly (C inserted before a Gly codon, giving CGG+N then +1) | 146/150 | 0 | consensus G128, no consensus 16, not found 6; pairs G694 S11 A6 … |
| cgg_arg (Arg codon replaced by CGG, then one base inserted) | 99/100 | 0 | consensus R84 K2, no consensus 13; pairs R506 K35 … |

The called event lies 0 codons from the edit in 58–75% of cases and within 2 codons in 92–96%. The residue readout therefore separates "the CGG+N unit is Gly" from "the CGG+N unit is Arg" at the level of single sites.

**Null codons** (tables 2 and 2b):
- **Frameshifts.** In the five frameshift genomes, CGG is followed by a +1 frameshift at 80–98% of aligned sites. CGA, CGC and AGG are followed by one at 0–2.3% of their in-frame loci (G7, the noisiest MAG, is highest: CGA 4.4%, CGC 3.4%, AGG 2.1%). In G2, G5, G6 and G8, CGG and the nulls alike are at 0–1% (G8: one CGG, in a 7-copy element).
- **Junction counts.** CGG sits in 56–85% of all frameshift events in G1, G3, G4, G7 and G9, against at most 2.1% for CGA and CGC and at most 8.2% for AGG.
- **Residues.** The consensus residue at CGA, CGC and AGG is R in every genome.
- **The nucleotide after CGG** follows Gly-codon third-base usage (`ggn.py`). Fraction A / T / C / G after CGG versus at the third position of Gly codons:
  - G1: .58/.30/.05/.07 versus .61/.26/.05/.08.
  - G4: .70/.20/.03/.07 versus .63/.28/.05/.04.
  - G7: .70/.11/.06/.14 versus .57/.25/.06/.12.
  - This fits "C, then a Gly codon".

## 4. Surprises and caveats

**Surprises**
- **What the +1 genomes share.** None of G1, G3, G4, G7 or G9 has a tRNA-Arg CCG or TCG, in either the parquet tables or the blastn cross-check. Among the Arg genomes, G5 (2 loci) and G6 have CCG and G8 has TCG. G2 has neither, but its tRNA detection is clearly incomplete (next point). This fits the residue and frame evidence: with no tRNA to decode CGG, the 4 nt CGGN end up encoding Gly.
- **G2's tRNA set is incomplete.** Its tRNA set has only tRNA-Arg ACG. blastn with every tRNA from the 9 genomes also finds no tRNA for AGA or AGG, yet AGA is G2's main Arg codon (consensus R5406). So G2's tRNA complement is incompletely detected, and its tRNA data cannot be used.
- **G2's CGGs are mostly in near-identical DNA.** G2 avoids CGG more strongly than the other three Arg genomes (in/out 0.019). Many of its CGGs sit in windows 0.86–1.0 identical to the Cdif, Chon or Bsub reference, where the reference also has CGG. Only 14 G2 CGGs lie in more diverged sequence, and there the residues are mixed (R4 L3 K2 Q2 T1 N1, plus one more).
- **G2's Arg verdict therefore rests on three things:**
  - no frameshift at any of 26 CGGs;
  - no stop-like behaviour;
  - R at every near-identical locus and in both consensus calls.

  Read it as "standard Arg, thin evidence".
- **G8's residue support is weaker.** The R share of clean pairs at CGG is 0.36, against 0.58 for CGA. The best-reference residue shows a K excess, a third of which comes from one ~10 kb cluster. The consensus calls (R4 K1), the tRNA-Arg TCG, the lack of frameshifts and the in-frame occurrence all point to Arg. I would not call a reassignment to Lys on 5 consensus sites.
- **G8 has one 7-copy element with CGG then +1.** It is collapsed to 1 site, faces S, and is preceded by CCT/ACT. It is a repeat-specific case (possibly a pseudogene family), not general CGG behaviour.
- **G4 has TGA in frame at 3 unique loci** (consensus R1). This is an unrelated curiosity, and too few to interpret.
- **G9's tRNA table has two odd entries**: a "???" with a 4-nt anticodon CACC and one with "GA". Neither has the form NCCG, which a tRNA reading CGGN as a quadruplet would need.

**Same-frame CGG in the +1 genomes.** These are a small minority:
- G1: 0 of 181.
- G3: 1 of 14.
- G4: 6 of 46.
- G7: 13 of 299.
- G9: 6 of 35.

They are mostly single-reference loci at low identity, where the frameshift may simply be misplaced, plus short near-identical segments that look transferred:
- G7 contig WQVL01000047.1: 2.9 kb, GC 0.36 against 0.33 for the genome. Three CGGs there face R (best identity 0.50–0.93).
- G9 contig CASDSU010000015.1: GC 0.344 against 0.311.
- G4 1429–1430 kb: 0.64–1.0 identity to the Chon reference, which itself has CGG there.
- G4 1151 kb: 0.8 identity to the Cdif reference (residues W and Q).

**Caveats**
1. **DNA alone cannot tell two mechanisms apart.** One is skipping the C and reading GGN with the normal tRNA-Gly; the other is decoding the 4 nt CGGN as a quadruplet. Both put Gly at that position, and that is what I call "+1 frameshift, Gly incorporated".
2. **The exact position of the shift is uncertain.** It is fixed by alignment only to within about ±1–2 codons (planted control: 58–75% exact). The residue readout handles this with a local re-alignment of the explicit hypothesis translation, and the cgg_arg and cgg_gly controls show it returns the right residue.
3. **Background frameshifts.** MAG sequencing errors and pseudogenes give non-CGG frameshifts in every genome (the null rows). G7 is noisiest, but CGG still dominates there: 288 of 352 events.
4. **Possible bias of the scoring matrix.** miniprot scores CGG against reference X/R with the standard code. If anything this favours aligning CGG as Arg in frame, not the Gly-plus-frameshift pattern that was found.
5. **G3 has few sites.** It has 14 aligned CGG. The frameshift is clear (13/14, and 14/18 events), but the residue is only Gly-leaning (6 of 14 clean pairs; one consensus call, A). Its residue is undetermined.
6. **MAG tRNA tables are small** (25–38 tRNAs in G1, G3, G4, G7 and G9), so a missing CCG or TCG tRNA is supporting evidence only. G6 is also a small MAG (28 tRNAs) and does show CCG.
7. **The reference set is limited.** It is 8 Firmicutes, which are distant from these mostly uncultured Clostridia. Many loci have only 1–2 references at identity 0.3–0.6, which is why consensus counts are well below locus counts.
8. **Process notes.**
   - One command accidentally invoked system `python3` with an empty script. It did nothing; everything else ran through `uv`.
   - After all analyses and these verdicts were done, I briefly saw the last lines of a sibling agent's task output in this session's task directory. The visible part showed only script names (a KACI/GTDB join for TANB77 and a genetic-code-table check), and no genome labels or verdicts. Nothing here was changed because of it.
   - I did not read anything else outside `recheck_cgg/`, apart from the codon-atlas genomes, the trna.parquet files and Swiss-Prot (for the three RF2 sequences).
