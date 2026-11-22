# Audit: "mitochondrial AGG = methionine in *Porphyridium aerugineum*" (post draft 2026-11-11)

Independent auditor, 2026-11-11. Draft: `lab/outbox/2026-11-11-mito-agg-methionine/post.md` + `image.png`.

**Verdict: blocking errors: YES. Four statements must be fixed before posting, and all four are text fixes. The finding
itself holds:** every independent recomputation reproduces the AGG↔ATG pattern (AGG never pairs with an arginine codon
in either direction or under any of three alignment methods), the reads, and the blinded panel. The errors are:

1. **Wrong numbers:** "319 red algal mitochondrial genomes ... 1,169". The same census gives **390 red algal genomes (389
   with aligned sites) and 1,418**. The author's own working table (`scratchpad/rhodo_agg.tsv`, 14:48) also gives
   390 / 7 / 1,418.
2. **Wrong citation:** "Lee et al. 2022, *BMC Biol.*" should be **Kim et al. 2022, *BMC Biol.* 20:2** (authors Kim D, Lee J,
   Cho CH, Kim EJ, Bhattacharya D, Yoon HS).
3. **The Fisher p does not match the test as worded.** p 7e-14 is the table [[17,0],[0,34]] (ATG vs Arg, with the 7
   "other" sites dropped). With the fractions printed ("17 of 21 ... 34 of 37") the test gives p 3e-11 (ATG vs not) or
   2.6e-8.
4. **Overclaim about *P. purpureum*:** "reads AGG normally" / "as arginine". MT483996.1 has AGG at only **2** aligned
   sites, and the blinded panel says only "not Met".

Should-fix items (non-blocking) are listed in section 5.

Scripts: `scripts/` (pair_nt.py, fisher.py, census_check.py, red_by_file.py, nt_identity.py). Outputs: `fwd/`,
`rev/`, `census_check.txt`, `red_per_genome.tsv`, `nt_identity.txt`. GenBank records were fetched by efetch into
separate `mktemp -d` dirs under the scratchpad (gb.k797/NC_062302.1.gb, gb.055I/MT483995.1.gb, gb.T410/MT483996.1.gb).

---

## 1. Code audit: `src/newleads/codon_pair.py`

**Result: no bug found.** Frame, direction, CDS extraction and gene matching are correct.
- **CDS extraction and frame.**
  - Every CDS in both records has codon_start 1 and a length that is a multiple of 3.
  - Neither record has duplicate gene names.
  - My own table-4 translation of `f.extract()` equals the record's `/translation` for 17/17 (SAG) and 18/18 (CCMP1948)
    CDSs, so the extraction handles introns and joins correctly.
- **Gene matching.** The 14 shared genes are atp6, atp9, cob, cox1–3, nad1–6, nad4l and rps3. The header of
  `isolate_pair.txt` says "13 shared CDS + orfs", which is a cosmetic inconsistency. The figure footer ("14") is right.
- **Codon pairing.** Codon i of A is paired with codon j of B through `aln.aligned` blocks of a global BLOSUM62 protein
  alignment, and the codons are sliced as `[3i:3i+3]`. The direction is correct: A = SAG 1380-2, B = CCMP1948.
- **Caveats (minor).**
  - The protein alignment uses table-4 translations, in which AGG = R. This slightly favours R–R pairings over
    R(AGG)–M, which biases *against* the claim.
  - Only the first optimal alignment is used.

**Independent recomputation.** I aligned at the nucleotide level with MAFFT `--globalpair` and with a Biopython global
nt aligner. A codon is paired only if all 3 bases align gap-free to one in-frame codon. As a third route I used a MAFFT
protein alignment with AGG masked as X in both isolates.

| Quantity | Post / codon_pair | MAFFT nt | Biopython nt | AGG-masked protein |
|---|---|---|---|---|
| Codon identity | 2850/4320 (66%) | 2847/4315 (66.0%) | 2847/4308 (66.1%) | 2848/4320 (65.9%) |
| CCMP at SAG AGG sites: AGG / ATG / other / Arg | 37 / 17 / 4 / **0** (58) | 37 / 16 / 4 / **0** (57) | 37 / 16 / 4 / **0** (57) | 37 / 16 / 5 / **0** (58) |
| CCMP at SAG AGA+CGN (78): same / other Arg / ATG or AGG / else | 41 / 34 / 0 / 3 | 41 / 34 / 0 / 3 | 41 / 34 / 0 / 3 | 41 / 34 / 0 / 3 |
| CCMP at SAG ATG: ATG / AGG | 42 / 16 (71) | 42 / 16 (70) | 42 / 16 (70) | 42 / 16 (71) |
| Reverse direction (CCMP's AGG → SAG): AGG / ATG / other / Arg | 37 / 16 / 4 / 0 (57) | 37 / 16 / 4 / 0 | 37 / 16 / 4 / 0 | 37 / 16 / 4 / 0 |

- **ATG 17 vs 16.** The one difference is SAG rps3 codon 131. codon_pair pairs it with ATG, but it sits in an indel region
  (`PQTr---LNV` / `PQTmFCLLNS`) that the nt aligners leave unpaired, so the count is 16–17 depending on alignment.
- **A conservative miss in cox1.** The CCMP1948 cox1 annotation drops one codon at three exon junctions (CCMP codons 72,
  149, 293). At the junction "IVWAHH|YTVGRD" the genomic sequence after the CCMP exon starts "AGG…", with the same
  intron start as SAG. So CCMP1948 very likely also has AGG at SAG cox1 codon 305, the strict-Met column cox1:310.
  codon_pair misses this site, which makes the true tally more likely 38 AGG of 59.
- **The isolates are far apart.** Over the 14 genes the two isolates share 87% nt identity and 88% aa identity
  (`nt_identity.txt`). That is deep for "two isolates of one species", and they may be distinct species. This is
  informational only: it makes the shared AGG pattern stronger evidence, not weaker.

**Fisher test** (`scripts/fisher.py`, one-sided):

| 2x2 table | p |
|---|---|
| [[17 ATG, 0 Arg],[0 ATG, 34 Arg]], changed sites with "other" dropped | **6.8e-14**: the post's number |
| [[17 ATG, 4 not],[0 ATG, 37 not]], ATG among all changed sites | 3.0e-11 |
| [[17, 4],[3, 34]], the literal "17 of 21 vs 34 of 37" | 2.6e-8 |
| Audit nt alignment [[16, 4],[0, 37]] | 8.4e-11 |

The post prints the fractions of one table and the p of another. The conclusion is unchanged in every version.

## 2. Census numbers (`scripts/census_check.py`; my own strict-column recount, streaming the 3 census files)

My recount re-implements the strict-column definition: each genome is weighted 1/(genomes sharing its taxonomy minus
the last rank), and a column is strict when its top residue has ≥ 0.9 of the weight. It runs on mito + gbmito +
gbmitoorf (18,015 + 1,275 + 3,408 = 22,698 genomes, which matches "about 22,700").

- **SAG 1380-2.** Reproduced exactly:
  - strict Met: ATG 13, AGG 7 (plus AGT 1, ATA 1);
  - strict Arg: AGG 0, AGA 13;
  - 52 AGG sites in total; the consensus at them is M 19, L 12, I 4, ….
- **Relatives at SAG's 52 AGG columns.** Reproduced exactly: 233 codons, ATG 84, Arg codon 1. The single Arg codon is CGA
  in NC_062303.1 rps3, at an R-consensus column that is not strict.
- **All red algae.** Not reproduced.
  - The census has 390 red algal genomes, 389 with sites. Their declared tables are 4 (372), 1 (17) and 3 (1).
  - AGG sits at strict Met 7 times (all in NC_062302.1) and at strict Arg **1,418** times.
  - No subset I tried gives 319 / 1,169: table 4 only gives 372 / 1,231; mito + gbmito gives 363 / 1,336; mito +
    gbmito with table 4 gives 348 / 1,161; Florideophyceae only gives 334 / 1,134.
  - No file under the project holds 319 or 1,169.
- **The relatives' own AGG use** (`red_per_genome.tsv`; consensus counts are from `census_check.txt`):

  | genome | AGG sites | at R-consensus | at strict Arg | at strict Met |
  |---|---|---|---|---|
  | *P. purpureum* MT483996.1 | 2 | 2 | 1 | 0 |
  | NC_062401.1 | 8 | 4 | 0 | 0 |
  | NC_062400.1 | 6 | 5 | 2 | 0 |
  | NC_062303.1 | 11 | 7 | 2 | 0 |
  | NC_062397.1 | 9 | 8 | 4 | 0 |

  AGG is rare in every relative. In *P. purpureum*, the evidence for "reads AGG as arginine" is two codons.
- **CCMP1948 (MT483995.1) is not in the census.** Its "AGG at Met positions" rests on the blinded panel (Codetta M,
  25 columns; Method 2 n 50) and on the pairing above.

## 3. Reads (`results/agg_met/reads/`, made by `src/codes_census/reads_codons.py`)

**Script.**
- It runs minimap2 `-ax sr` against the gene nt sequences.
- It drops flag 0x904 reads (unmapped, secondary, supplementary).
- MAPQ ≥ `MIN_MAPQ`: the default is 20, and the README says 0 was used, which no saved file records.
- With `MAX_DIV`, reads with NM above that fraction of read length are dropped.
- It counts only reads covering all 3 bases of a codon.
- Two minor issues:
  - Each mate counts as a read, so overlapping mates count one fragment twice. This inflates n but not shares.
  - The per-site "other" column keeps only the top 2 non-reference codons (`most_common(3)` minus the reference).

**Run metadata (SRA).**
- DNA, SRR15114507: *P. aerugineum*, library "Paer", WGS, HiSeq 2500, BioProject PRJNA744153. That is the van Beveren
  et al. 2022 project, from the ESE lab in Orsay.
- RNA, SRR1296983: MMETSP0313, "Whole transcriptome sequencing of *Porphyridium aerugineum* SAG 1380-2", TruSeq mRNA
  (poly-A), biome freshwater, 0 psu.
- Both runs are strain SAG 1380-2.

| Claim | File | Verdict |
|---|---|---|
| 7 strict sites, 1,253–1,837 reads each | cob67 1253, cox1-78 1734, cox1-289 1634, cox1-304 1275, nad4-313 1787, nad5-103 1837, nad5-123 1809 | supported |
| ≥ 99.7% AGG at each | 99.68% (cob67 1249/1253), 99.71, 99.72, 99.78, 99.82, 99.89, 100% | true only after rounding: the minimum is 99.68% |
| 30,855 of 30,909 (99.8%) at the 19 Met-consensus AGG sites | recomputed from tsv: 30,909 / 30,855 | supported |
| "4 say ATG" | ATG itemized 4 times; 54 non-AGG reads, of which 35 are itemized and 19 are not (top-2 truncation) | lower bound: ATG is 4–23 |
| RNA 201 reads on the genes | only in README (stderr was not saved) | unverifiable from files |
| RNA 13 of 13 at 5 strict sites | 7 + 2 + 2 + 1 + 1, all AGG; nad5 sites 0 reads; Met-consensus total 33 / 33 | supported |

## 4. Figure (`image.png`, `src/agg_met/fig_agg.py`) vs `isolate_pair.txt`

| Row | Bars | Total |
|---|---|---|
| AGG | 37 / 17 / 4 | 58 |
| ATG | 42 / 16 / 13 other (ATA 5, GTA 3, 5 singles) | 71 |
| AGA, CGN | 41 / 34 / 3 | 78 |

- All rows match `isolate_pair.txt` exactly.
- The footer "14 shared mitochondrial genes" is correct.
- The 3 "other" in the AGA/CGN row has no number label (cosmetic).
- The figure inherits the 17-vs-16 alignment dependence noted in section 1.

## 5. Per-sentence fact-check of post.md

Grades are SUPPORTED, WRONG, UNSUPPORTED or OVERCLAIM. "Scratch" means
`scratchpad/`.

| # | Sentence (abridged) | Verdict | Source / note |
|---|---|---|---|
| 1 | In *P. aerugineum* mitochondria AGG most likely means Met, not Arg. | SUPPORTED | §1–3; blinded panel `blind/report.md` (R06, R02 = Met). "Most likely" is appropriately hedged. |
| 2 | No red alga was known to read AGG that way. | SUPPORTED | `notes/porphyridium_aggmet_lit.md` searches; my searches (2026-11-11) found no red algal AGG change. |
| 3 | Genes are read three letters (one codon) at a time. | SUPPORTED | Textbook. |
| 4 | Standard code: 61 sense codons for 20 aa, 3 stops. | SUPPORTED | Textbook / NCBI table 1. |
| 5 | Some lineages, many of them mitochondria, read codons differently; NCBI lists 26 alternatives. | SUPPORTED (nit) | gc.prt / Biopython list 27 tables: 1–6, 9–16, 21–33, so 26 non-standard. Nit: table 11 has standard assignments and differs only in start codons. |
| 6 | Red algal mitochondria use a code in which AGG is Arg. | SUPPORTED | NCBI table 4 (AGG = R). In the census, 372/390 red algal records declare table 4, 17 table 1 and 1 table 3, all AGG = R. Kim et al. 2022 annotate with table 4. |
| 7 | *P. aerugineum* is a single-celled freshwater red alga. | SUPPORTED | Kim et al. 2022 ("unicellular red alga *Porphyridium*"); MMETSP0313 biosample: "Aquatic: freshwater", salinity 0 psu. |
| 8 | Two isolates sequenced by two labs: SAG 1380-2 (NC_062302.1), CCMP1948 (MT483995.1). | SUPPORTED | GenBank submitters: van Beveren et al. (Orsay) and Kim et al. (Sungkyunkwan). NCBI nuccore holds only these two (MZ748310 = NC_062302.1). Note: 87% nt identity, so possibly two species (§1). |
| 9 | In both, AGG sits at Met positions; where genes differ AGG swaps with ATG, never an Arg codon. | SUPPORTED | SAG: census. CCMP1948 is not in the census, so its support is the blinded panel R02 plus the pairing. Arg = 0 in all routes and both directions. |
| 10 | At strict Met positions SAG has ATG 13, AGG 7. | SUPPORTED | `census_check.txt`; also AGT 1, ATA 1. |
| 11 | At strict Arg: AGG 0; AGA 13. | SUPPORTED | `census_check.txt`. |
| 12 | Across **319** red algal mitogenomes, AGG at strict Met 7 (all this genome), at strict Arg **1,169**. | **WRONG** | Recount: **390** red algal genomes (389 with sites), 7 (all NC_062302.1), **1,418**. Author's own scratch `rhodo_agg.tsv` agrees (390 / 7 / 1,418). No file holds 319 / 1,169. The qualitative point stands. |
| 13 | Relatives at SAG's AGG positions: ATG 84 of 233, Arg once. | SUPPORTED | `census_check.txt`: 37/14, 47/19, 50/20, 49/18 (+1 CGA), 50/13. |
| 14 | Same codon at 66% of aligned positions. | SUPPORTED | 2850/4320; audit 2847/4315. |
| 15 | At SAG AGG: CCMP AGG 37, ATG 17, other 4, Arg 0 (58). | SUPPORTED (caveat) | `isolate_pair.txt`. Independent nt alignments give ATG 16 of 57 (rps3 131 is ambiguous). A cox1 annotation gap hides one more AGG→AGG site. Arg 0 is robust. |
| 16 | At SAG AGA/CGN (78): different Arg codon 34, ATG/AGG never. | SUPPORTED | All three audit routes. |
| 17 | Changed AGG → ATG 17 of 21; changed Arg → Arg 34 of 37; one-sided Fisher p 7e-14. | **WRONG** (statistic/wording mismatch) | 7e-14 is [[17,0],[0,34]], with "other" dropped. The fractions as written give 3e-11 (ATG vs not) or 2.6e-8. See §1. |
| 18 | Two genomes from two labs agreeing also rules out an assembly error in one. | SUPPORTED (soften) | An error in one genome cannot make a pattern shared by both. CCMP1948 has no public reads (SRA lists only the two SAG 1380-2 runs), so its codons are not read-checked. Better to credit the SAG reads too. |
| 19 | DNA: each of 7 strict sites reads AGG in at least 99.7% of 1,253–1,837 reads. | SUPPORTED after rounding | Minimum is 99.68% (cob67, 1,249/1,253). Strictly, "at least 99.6%". |
| 20 | All Met-consensus AGG: 30,855 of 30,909 (99.8%) say AGG and 4 say ATG. | SUPPORTED / "4 ATG" is a lower bound | The tsv lists only the top-2 other codons per site; 19 of 54 non-AGG reads are not itemized. |
| 21 | RNA (poly-A transcriptome, only 201 reads on these genes): 13 of 13 AGG at 5 covered strict sites. | Poly-A and 13/13: SUPPORTED. 201: UNSUPPORTED by files | SRR1296983 = MMETSP0313, TruSeq mRNA, SAG 1380-2. 201 appears only in the README; mapping stderr was not saved. |
| 22 | So AGG is in both DNA and mRNA; this is not RNA editing. | First clause SUPPORTED; second OVERCLAIM (mild) | Excludes mRNA editing of these codons, in only 33 RNA reads at Met-consensus sites. Says nothing about tRNA editing or modification, which is a candidate decoding route. → "not mRNA editing". |
| 23 | Blinded re-check (Codetta + own Swiss-Prot likelihood, 9 coded genomes): AGG = Met in both isolates, none of 5 relatives, and the *Chromochloris* control. | SUPPORTED | `blind/report.md`, `calls.md`, `blind_prereg.md` (rule fixed 15:00 before launch), `blind_key.tsv`. The post-hoc v3 filter removed only R05 TAA→Met, which does not touch AGG. |
| 24 | The tRNA is unresolved. | SUPPORTED | |
| 25 | SAG has a CCU-anticodon tRNA and no UCU; its relatives have UCU. | SUPPORTED (nit) | Scratch `porph_trnas.tsv`. The CCU tRNA at 12459–12535 is GenBank's tRNA-Arg; none of SAG's 9 annotated tRNAs is UCU. UCU appears in **4 of 5** relatives (none found in NC_062397.1). The tRNA tables live only in scratch; copy them into `results/`. |
| 26 | CCMP1948 has neither, but carries a second Met-like tRNA (anticodon AUG). | SUPPORTED but misleading | Scratch `pa2_trnas.tsv` and GenBank: tRNA at 11348–11422 with tRNA-Met stems (TGATGAG/CTCATTA). Its anticodon 5′-AUG-3′ pairs with the **His codon CAU**, not with AGG. GenBank calls it `trnH-AUG` (tRNA-His). Readers will take "AUG" as the Met codon. |
| 27 | The decoder may be imported from the cytoplasm. | SUPPORTED (as speculation) | The SAG mitogenome encodes only 9 tRNAs (GenBank), so most tRNAs must be imported anyway. |
| 28 | No NCBI code (tables 1–33) has AGG = Met. | SUPPORTED | Biopython / gc.prt: AGG is R, S, G, K or stop in every table. |
| 29 | Known in *Chromochloris zofingiensis* (Noutahi 2019 MBE; Žihala & Eliáš 2019 GBE) and in uncultivated Bacilli and their phages (Shulgina & Eddy 2021 eLife). | SUPPORTED (nit) | Noutahi et al. 2019 MBE 36:766: "In Chromochloris, we identify tRNA variants decoding AGG as methionine". Žihala & Eliáš 2019 GBE 11:2992: "Arg-to-Met reassignment of the AGG codon". Shulgina & Eddy 2021 eLife 10:e71402: 8 Bacilli genomes; "Two phage assemblies were predicted to translate AGG as methionine" (same baboon-faecal samples, "may have adapted"). Nit: "their phages" → "two phage assemblies from the same samples". |
| 30 | I found no report of a changed AGG in any red alga. | SUPPORTED | As #2. |
| 31 | Neither genome paper (van Beveren et al. 2022 GBE; **Lee** et al. 2022 BMC Biol.) mentions a code change. | **WRONG author**; code claim SUPPORTED | BMC Biol 20:2 is **Kim** D, Lee J, Cho CH, Kim EJ, Bhattacharya D, Yoon HS 2022. Its full text uses "translation_table 4" and mentions no code change. van Beveren et al. 2022 GBE 14:evac037 says only "using the applicable genetic code"; *P. aerugineum* is not named in its main text (it is in Table S1), but SRR15114507 belongs to its BioProject PRJNA744153. |
| 32 | A Codetta survey of mitochondrial genomes (Shulgina & Eddy 2023) did not report one. | OVERCLAIM (weak evidence) | Shulgina & Eddy 2023 *Bioinformatics* 39:btac802 benchmarked **82 mitogenomes with known alternative codes**. Its unexpected calls were AUA→Met in 2 green algae and AAG→Met in a chytrid; it names no red alga. Not shown that *P. aerugineum* was in the set (likely not; NC_062302 dates from 2022). "Survey" suggests coverage it may not have had. |
| 33 | Census: ~22,700 mitogenomes (all RefSeq + one per extra species from GenBank), declared codes, profile alignment, dominant-residue positions. | SUPPORTED | 18,015 + 1,275 + 3,408 = 22,698 in the census files. "All of RefSeq" was not checked against a release count. |
| 34 | This genome was flagged by that scan. | SUPPORTED | `porph_agg_columns.tsv`, census. |
| 35 | Second isolate compared codon by codon. | SUPPORTED | codon_pair.py. |
| 36 | Raw reads (DNA and RNA) mapped, codon read at each site. | SUPPORTED | reads_codons.py. |
| 37 | Blinded re-check with two methods and a pre-written pass rule. | SUPPORTED | `blind_prereg.md`. |
| 38 | Not known: whether the proteins carry Met. | SUPPORTED | |
| 39 | No mass-spec data exist for *P. aerugineum*; the public sets are for *P. purpureum*, **which reads AGG normally**. | First part SUPPORTED (PRIDE only; soften to "I found none"); last clause OVERCLAIM | PRIDE "Porphyridium": PXD052401 and PXD011275, both *P. purpureum*. MassIVE and jPOST were not checked. *P. purpureum* has AGG at 2 aligned sites; the panel says only "not Met". |
| 40 | Which tRNA reads AGG, and whether AGG is read only as Met or sometimes Arg. | SUPPORTED | |
| 41 | *P. purpureum* and the other relatives sampled **read AGG as arginine**. | OVERCLAIM | AGG is rare in all five: 2, 8, 6, 11 and 9 aligned sites. It sits at R-consensus columns 26/36 times, strict Arg 9, strict Met 0; the panel says "not Met". The data say "rarely use AGG; where they do, mostly at Arg positions, never at conserved Met". |
| 42 | Other *P. aerugineum* strains are unsequenced. | SUPPORTED | NCBI nuccore (mitochondrion) holds only SAG 1380-2 and CCMP1948. SRA holds only the 2 SAG 1380-2 runs. |
| — | Figure | SUPPORTED | §4. |
| — | "Note for Ryan" block | not checked | Internal, and still PENDING. |

## 6. Suggested wording fixes

**Must fix:**
- **#12.** "Across 390 red algal mitochondrial genomes, AGG sits at conserved methionine positions 7 times (all 7 in this
  genome) and at conserved arginine positions 1,418 times."
- **#17**, either of:
  - "Changed AGG sites go to ATG 17 times and to an arginine codon never; changed arginine sites go to another arginine
    codon 34 times and to ATG never (one-sided Fisher p 7e-14)", or
  - "ATG at 17 of 21 changed AGG sites vs 0 of 37 changed arginine sites: one-sided Fisher p 3e-11."
- **#31.** "Kim et al. 2022, *BMC Biol.*"
- **#39 and #41.**
  - #39: "...the public sets (PRIDE) are for *P. purpureum*, which rarely uses AGG."
  - #41: "*P. purpureum* and the other relatives sampled rarely use AGG, and where they do it sits mostly at arginine
    positions, never at conserved methionine; the blinded check did not call it methionine in any of them."

**Should fix:**
- **#19.** "...reads AGG in 99.7% or more (lowest: 1,249 of 1,253) of 1,253 to 1,837 reads", or "at least 99.6%".
- **#20.** "...30,855 of 30,909 reads (99.8%) say AGG; of the other 54, at least 4 say ATG". Alternatively, rerun the
  tally without the top-2 truncation.
- **#22.** "So AGG is in both the DNA and the messenger RNA: it is not created by mRNA editing."
- **#26.** "...a second tRNA with a methionine-tRNA body but anticodon AUG (which would pair with the histidine codon
  CAU, not AGG; GenBank calls it tRNA-His)."
- **#32.** "...and Codetta's 82-mitogenome benchmark (Shulgina & Eddy 2023, *Bioinformatics*) reported no such change".
  Or drop it, since *P. aerugineum* was probably not in that set.
- **#15.** Optionally: "ATG 17 (16 with nucleotide alignment)". The cox1 annotation gap in CCMP1948 hides one more
  shared AGG.
- **#29.** "...and two phage assemblies from the same samples".
- **#18.** Optionally: "Two independently assembled genomes agree, and SAG 1380-2's raw reads confirm its codons."
- **Provenance.**
  - Save the mapping stderr: the 201 and 213,142 read totals are not in any file.
  - Copy the tRNA tables (`porph_trnas.tsv`, `pa2_trnas.tsv`) and `rhodo_agg.tsv` from scratch into `results/agg_met/`.
  - Fix the "13 shared CDS" header in `isolate_pair.txt`.

## Sources (web)
- Kim et al. 2022 BMC Biol 20:2: https://pmc.ncbi.nlm.nih.gov/articles/PMC8742464/ (full text via Europe PMC); GenBank MT483995.1 REFERENCE.
- van Beveren et al. 2022 GBE 14:evac037: https://pmc.ncbi.nlm.nih.gov/articles/PMC8995046 (full text via Europe PMC).
- Noutahi et al. 2019 MBE 36:766 (abstract via search); Žihala & Eliáš 2019 GBE 11:2992: https://pmc.ncbi.nlm.nih.gov/articles/PMC6821328
- Shulgina & Eddy 2021 eLife 10:e71402: PMC8629427 full text via Europe PMC (https://elifesciences.org/articles/71402)
- Shulgina & Eddy 2023 Bioinformatics btac802: https://pmc.ncbi.nlm.nih.gov/articles/PMC9825746/
- SRA runinfo for SRR15114507 and SRR1296983 (NCBI eutils); PRIDE API search "Porphyridium"; NCBI nuccore esearch for
  *P. aerugineum* mitochondrion.
