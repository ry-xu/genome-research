# Pre-post check: CUU and CUA read as proline in four yeast mitochondria (cun-pro) (independent check, 2026-11-11)

Claim checked, from the outbox draft `outbox/2026-11-11-mito-proline-codons/post.md` (Note dated 12:44):
- In the mitochondria of four Trichomonascaceae yeasts (*Wickerhamiella infanticola* BK074532, *W. siamensis* BK074538, *W. verensis* BK074542, *Sugiyamaella pinicola* BK074262), CUU and CUA mean proline, not leucine.
- "No known genetic code reads any codon other than CCN as proline."
- The change arose twice: once in the three *Wickerhamiella*, once in *S. pinicola*.

I read the post and its image, README lines 37–154, `results/cun_pro/` (audit report, fact-check, blinded re-check report, read tables, ARAGORN calls, grid TSV, the Wolters table and OR703884 alignment in `lit/`), the scripts in `src/codes_census/` that build the numbers (sitescan_org, codon_contrast, pro_leu, pro_gene, pro_grid), the tree files in `results/wickerhamiella_tree/`, and ENA/NCBI metadata for every run. I re-derived the numbers with my own loader and strict-column definition, which is separate from the census code (scripts and outputs in `prepost/cun-pro_work/`). A separate agent ran the literature search. I then checked its key citations myself: Ling 2014 full text, Vargas-Rodriguez 2021, Schuntermann 2023, NCBI Taxonomy, and Wolters Table S1. I also ran four Europe PMC searches of my own. I did not re-map raw reads or rebuild the tree. I edited nothing in `/Users/ryan/lab` outside `prepost/` and messaged no agent.

## Verdict: FIX, THEN POST.

- **The finding survives.** Every headline number I re-derived matches (section 2), and none is wrong in a way that changes the conclusion. The signal does not depend on the census settings:
  - It is unchanged when the four genomes, or all Dipodascomycetes, are dropped from the consensus.
  - It holds at 95% and 99% thresholds.
  - Of the CUU/CUA codons the four have at strongly conserved positions of any amino acid, 286 of 291 sit at proline positions.
- **Novelty holds for the headline** (section 1). No NCBI table or reported genetic code reads anything but CCN as proline, and no one has reported this in these yeasts.
- **One "Known before" sentence is false** as written: "found no report of any other codon read as proline". *Streptomyces* bacteria carry proline tRNAs that put proline at some alanine, asparagine or threonine codons. That is partial mistranslation, not a code change, so the headline stands, but the sentence must change.
- **The Ling 2014 summary is wrong** on where the *Ashbya* tRNA came from. Ling says it most likely came from the threonine tRNA, not directly from tRNA-His.
- **Three number wordings overstate or mislabel** (section 5):
  - *S. pinicola*'s "36 of its 80" uses the wrong denominator;
  - the DS02 "1 of 1,942 reads" is a lower bound from a truncated list;
  - "never at conserved leucine positions" is true, but it is not evidence, because the relatives barely use CUN for leucine either.
- **The main risk is still RNA editing** (section 3). The census itself contains the one known lookalike: hornwort mitochondria, which edit U to C. The post handles editing honestly; it could say more.
- **Fixes before posting** (section 5): required fixes 1–6. All are wording changes; nothing needs re-running.

## 1. Novelty

The question has two parts: has anyone reported this for these yeasts, and has any codon other than CCN ever been reported as proline?

**These yeasts.**
- **Sarhan et al.** (*Sci. Data* 13:28; doi 10.1038/s41597-025-06447-x, online 18 Dec 2025, volume year 2026):
  - "annotated using 'MFannot' v1.37, using their corresponding genetic code";
  - the full text never mentions codons, proline, CUN, *Wickerhamiella*, *Sugiyamaella* or Trichomonascaceae;
  - the four BK records carry only a source feature, with no genes and no translation table.

  The post's "deposited without gene annotation" and "annotated each genome with its taxon's code" are both right.
- **NCBI Taxonomy:** I re-fetched all four species. Each has mitochondrial code 4 (MGCId 4; taxids 1116880, 1571656, 2714595, 378059). The post is right.
- **Wolters et al. 2023** (*Front. Microbiol.* 14:1268944):
  - "The translation table for each species was estimated using codetta v2.0". So the claim that they ran Codetta holds.
  - But their decision rule only chose between table 3 and table 4: "we used the Saccharomyces code when the CUN reassignment was supported and the fungal code for all others". A proline call for CUU/CUA could not have changed their table. Their raw Codetta output is not in the paper or the supplement. It may be on figshare (10.6084/m9.figshare.23799789), but the listed zips are not obviously Codetta output, and neither agent downloaded them.
  - So "Codetta missed it" is not established. The likelier reading is that Codetta flagged it and their rule ignored it. Our blinded re-check's Codetta called P for both codons in all four genomes. See wording fix 6.
  - Their Table S1 (in `results/cun_pro/lit/`): *W. infanticola* (OR703884) is the only one of the four included, and it is assigned table 4. It was assembled from SRR5938406, the same run as BK074532. So "identical" is expected; it is not a second assembly from independent data. *W. siamensis*, *W. verensis* and *S. pinicola* are absent.
  - The dN/dS quote checks out: "the highest values were found in the order Dipodascales for yeasts in the Wickerhamiella / Starmerella clade and the grouping of yeasts most closely related to that clade". But only one of their at least seven clade members is affected, so CUN-as-proline can explain little of the clade-level result. See optional fixes.
- **Christinaki et al. 2022** (*GBE* 14:evac073) took codes from NCBI Taxonomy ("…Trichomonascaceae… employ genetic code NCBI 4"). It has only two Trichomonascaceae and no mention of proline.
- **Other sources, no hits:**
  - Ahrendt et al. 2026 (*NAR*, about 10,000 JGI mitogenomes) has nothing on CUN or these genera.
  - Groenewald 2023 (*Stud. Mycol.* 105) confirms that *Wickerhamiella*, *Starmerella* and *Sugiyamaella* are in Trichomonascaceae, so the family name is right.
  - My own Europe PMC searches for these genera with "mitochondrial" plus "codon" or "genetic code" returned 50 hits. The only relevant one was Wolters.

**Any codon read as proline.**
- **NCBI genetic codes:** 27 tables, last updated 23 Sep 2024. Proline is CCN only in all of them (the subagent parsed the NCBI page; I checked Biopython's 27 tables). "26 alternatives" and "all 27 NCBI codes (tables 1–33, with gaps)" both hold.
- **New codes from 2021–2026** are none of them proline:
  - Shulgina & Eddy 2021;
  - Codetta (Bioinformatics 39:btac802, 2023);
  - Fakih 2026 (ACA→Asp);
  - Melnykov 2026;
  - Maturová et al. 2026 (ciliates);
  - Parks 2026 (UGA→Trp);
  - Barcytė et al. 2025;
  - Žihala & Eliáš 2019.

  My searches for "codon reassignment" with proline and for CUN with proline found only engineering papers.
- **Natural mistranslation to proline exists.** This contradicts post line 15. I confirmed both abstracts via Europe PMC.
  - Vargas-Rodriguez et al. 2021, *PNAS* 118:e2110797118: Streptomyces ProRSx with tRNA-ProA "can translate GCU alanine codons as both alanine and proline … the initial example of a natural tRNA synthetase/tRNA pair for dedicated mistranslation of sense codons".
  - Schuntermann et al. 2023, *J. Biol. Chem.* 299:104852: tRNAs with Asn or Thr anticodons in a proline-tRNA body: "these tRNAs translate asparagine and threonine codons with proline".

  Both are ambiguous decoding of a few codons by an extra tRNA, mostly shown in *E. coli*. Neither is a code change, and neither is in an organelle. The headline ("No known genetic code…") is safe; line 15's "no report of any other codon read as proline" is not. See wording fix 1.

**Other fates of CUN in yeast mitochondria.**
- **Ling et al. 2014:** *NAR* 42(1):499–508, doi 10.1093/nar/gkt842. The citation is correct.
  - "CUU and CUA codons correspond to alanine in A. gossypii", confirmed by mass spectrometry.
  - "the strict avoidance of CUG and CUC codons seems noteworthy".
  - **Error in the post (line 14).** On where the tRNA came from, I read the full text myself (PMC3874161; the tRNA names are inline images, so I identified them from context). The Ala tRNA "clearly clusters with yeast mitochondrial [tRNA-Thr], from which it most likely also derived via gene duplication". It is "an orthogonal Ala tRNA … either via duplication of [tRNA-Thr(UAG)] or its evolutionary precursor tRNA-His". tRNA-His is the ancestor of the Thr tRNA. So "most likely arose from tRNA-His" misstates it. See wording fix 2.
- **Other fates, none another amino acid:**
  - absence (*Kluyveromyces*, which "avoids CUN codons altogether"; Ling 2014; CoreTracker, Noutahi 2017);
  - a +1 frameshift at CUGY (*Kazachstania*);
  - "uncertain" Thr in *Hansenula saturnus* and *K. thermotolerans* (NCBI table 3 notes).

  So threonine and alanine are the only reported meanings, and proline would be the third. The post's framing ("changed meaning twice", "a third meaning") holds.
- **RNA editing in fungal mitochondria:** chytrids edit tRNAs only. Laforest, Roewer & Lang 1997 (*NAR* 25:626): "Editing was not found at other positions of the tRNAs, or in the mitochondrial mRNAs". No yeast mitochondrial mRNA editing is reported. "I know of no such editing in fungal mitochondria" is accurate.

**Bottom line:** the headline claims hold today: no genetic code reads a non-CCN codon as proline, nobody has reported this in these yeasts, and it would be the third meaning of CUN in yeast mitochondria. One supporting sentence (line 15) is false and one citation detail (line 14) is wrong.

## 2. Correctness: the headline numbers

Sources:
- "census" means my own loader over `data/codes/{mito,gbmito,gbmitoorf}.census.npz` (`prepost/cun-pro_work/census_check.py` → `census_check.out`).
- Reads come from `results/cun_pro/reads*/` via `reads_check.py` and `ds02_other.py`.

| post claim | source | check |
|---|---|---|
| About 22,700 genomes; 18,015 of 18,031 RefSeq plus one per additional species | census npz; ENA/NCBI counts | holds (22,698; 18,015 RefSeq). "One per additional species" is approximate. 4,683 GenBank records are unique names but about 3,980 binomials, 61 of which also have a RefSeq genome. Harmless. |
| Per genome, 44–82 CUU/CUA at strict-proline positions; 0 at strict-leucine positions; UUA/UUG at leucine 140–166 | census, threshold 0.9, family weight (each of the four weighs 1/105) | holds: CTN@P 78/82/82/44 (CTT 16/12/13/11, CTA 62/70/69/33, CTC and CTG 0); CTN@L 0/0/0/0; TTR@L 140/149/147/166. Lowest proline share at any CTN@P column is 0.903. At 0.95: 71/75/75/41; at 0.99: 37/37/37/19; CTN@L is 0 at every threshold. Unchanged when the four, or all 151 Dipodascomycetes, are dropped from the consensus. |
| *S. pinicola* "uses CCN at 36 of its 80 proline positions" | census | **the 36 holds; the denominator doesn't.** *S. pinicola* has 91 sites in strict-proline columns: CTA 33, CCA 20, CCT 15, CTT 11, CCC 1, and 11 other codons. 80 is only CCN plus CUN. See wording fix 4. The four genomes have 91, 91, 93 and 91 such positions. The post never gives the scale for "44 to 82". |
| Other 17 *Wickerhamiella*: CCN 80–85, CUN 0–1 | census | holds (only BK074539 has 1) |
| CUA Fisher test: the four at p 2.1e-28 to 1.9e-46; 147 other Dipodascomycetes, lowest 0.019 from a single CUA | census; one-sided Fisher [[CUA@P, CUA@L], [syn@P, syn@L]] | holds: 2.1e-28 (*S. pinicola*), 2.8e-39, 4.4e-45, 1.9e-46. Lowest other is *Geotrichum* BK072907 at 0.0186, [1, 2, 0, 158]; next is *Tardiomyces* BK074691 at 0.0205. 2 of 147 are below 0.05. |
| CUU: the four at 2.5e-6 to 1.8e-8; others lowest 0.049 | census | holds (1.8e-8, 2.8e-8, 9.2e-7, 2.5e-6; *Middelhovenomyces* BK073468 0.0486; 1 of 147 below 0.05) |
| Reads: 286 sites, all supported; 3,782–46,781 reads per genome; 98.9–99.6% agreeing | `results/cun_pro/reads/*.tsv` | holds: 78+82+82+44 = 286, all covered; reads 3,782 / 46,781 / 13,970 / 16,746; agreeing 98.86 / 99.60 / 99.57 / 98.92%. No site has a non-majority. *W. infanticola* has 2 sites with only 1 read. The choice of ORF copy doesn't change it. |
| PYCC 7069 (2024): CUU/CUA at all 82 sites, 1,432 of 1,450 reads | `reads_isolates/`; ENA ERR11532760 | holds (98.76%; lowest site 8 of 9). Thin at 3 sites, which have 1, 2 and 2 reads. ENA: *W. siamensis* PYCC 7069, Yeast Genomics Lab@Nova (Portugal), public 2024-05-07. |
| DS02 PacBio (2017): CUU/CUA commonest at all 65 sites with ≥ 3 reads | `reads_isolates/`; ENA SRR4787334 | holds: 67 of 78 sites have reads, 1,942 reads, 1,540 matching the assembly codon (79.3%); every one of the 65 sites with ≥ 3 reads has CTT or CTA as its commonest triplet. 61 of the 65 have the assembly's own codon as majority. |
| DS02: "a proline CCN codon shows up in 1 of 1,942 reads" | same | **lower bound only.** The "other" column keeps only the top 3 triplets per site, minus the codon itself. 119 of the 402 non-matching reads are not listed, so CCN could be among them. See wording fix 5. |
| CUC and CUG used 0–1 times per genome | census, all aligned sites (`cuc_cug.py`) | holds: CTC 0/0/0/0; CTG 0/0/1/1 (out of about 3,800–4,400 aligned codons per genome). Same as the fact-check. |
| 27 NCBI tables; Pro = CCN in all | Biopython; NCBI page (subagent) | holds |
| Blinded re-check: CUU at proline 75 of 76; CUA 298 of 300; Codetta calls P in all four | `results/cun_pro/blind_recheck/report.md` | holds (24+15+17+19 = 75; 78+41+91+88 = 298). The post discloses that the agent was steered toward CUN by the table 3 vs 4 framing and could see organism names. Good. |
| Tree: 49 genomes, 4,403 sites, CUN masked; *Wickerhamiella* trio at 100; *S. pinicola* inside *Sugiyamaella*, five nodes at SH-aLRT 98–100 and UFBoot 96–100; other 15 *Sugiyamaella* use CCN | `results/wickerhamiella_tree/`; `pro_grid_tree.tsv` | holds: SH-aLRT 98.1/99.9/100/99.5/100; UFBoot 96/96/100/99/100. The full run finished (4 "BETTER TREE FOUND" lines; end 15:48). The 15 other *Sugiyamaella* have CUN 0 at the figure's proline columns. The most recent common ancestor of all four spans 45 of 49 tips, so one origin would need many reversals; two origins is the parsimonious reading. |
| All four keep tRNA-Pro(UGG), tRNA-Leu(UAA) and a UAG tRNA; relatives that read CUN as Leu have a similar one | `results/cun_pro/dipo.aragorn.tsv` | holds. The three *Wickerhamiella* share an identical UAG tRNA; *S. pinicola*'s differs. Controls BK074539, BK074529 and BK074251 also have a UAG tRNA. "No more like tRNA-Pro than like other tRNAs": not re-checked. |
| OR703884 identical to BK074532 | `results/cun_pro/lit/OR703884_vs_BK074532.paf` | holds: both are 56,486 bp, the same circle rotated, with blocks of 36,482 and 20,004 bp matching exactly. Both come from run SRR5938406 (Wolters Table S1, sheet 3). |
| NCBI Taxonomy gives table 4 | efetch, all four taxids | holds |
| Figure: 81 of 126 columns, 12 proteins, 49 genomes | `pro_grid_tree.tsv`; `src/codes_census/pro_grid.py` | holds; every row sums to 81. The figure uses its own column filter, so *S. pinicola* shows 34 CCN / 42 CUN rather than the census's 36/44. That is consistent. |

**Extra check, not in the post:** I counted the four genomes' CUN codons at strict columns of any amino acid, not just proline and leucine (`cun_by_res.py` → `cun_by_res.out`). 286 of 291 are at proline columns; the other five are 4 at phenylalanine and 1 at alanine. This answers "maybe CUN just wanders" better than the leucine statement does (section 3a).

**README only, not in the post:** README line 80 gives the *S. pinicola* CCN-at-proline read tally as "12,123 of 12,348". That matches the read table's own strict-proline flags (36 rows). On my census site keys it is 11,888 of 12,122. One site differs between the two sets, the same kind of 81-vs-82 drift the audit found. Both are about 98%, so this is not a post issue.

## 3. Weak parts

**a. "Never at conserved leucine positions" is true but not evidence.**
- **Not circular.** The consensus is census-wide. Each of the four weighs 1/105 of a family. The four are read with table 4, which counts their CUN as leucine and so biases the census *against* the finding. The result is also unchanged when they, or all Dipodascomycetes, are excluded.
- **But not distinctive.** Relatives barely use CUN for leucine either:
  - among the 147 other Dipodascomycetes, 35 also have 0 CUN at strict-leucine positions, and 108 of all 151 have 7 or fewer;
  - close *Wickerhamiella*, *Starmerella* and *Sugiyamaella* relatives have 0–7, except BK074539 (*Wickerhamiella* sp.) with 93. Farther away, *Magnusiomyces* has 142–151.

  These yeasts spell leucine UUA/UUG. So the line 5 sentence ("never at conserved leucine positions. Leucine is spelled UUA/UUG instead.") reads as if CUN moved from leucine to proline in these four. In fact CUN was already nearly unused, which fits codon capture.
- **Where the evidence is.** It is all on the proline side: 78–82 CUN at proline positions, against 0–1 in 17 sister *Wickerhamiella* and no more than 3 in any of the other 147 Dipodascomycetes (the whole-census top list in `census_check.out` holds every genome with 4 or more, and none of them is a Dipodascomycete apart from the four). The new 286-of-291 count says it more directly. See wording fix 3.

**b. RNA editing is the main alternative, and the census contains a lookalike.**
- **The only other genomes with this signal edit their RNA.** In the whole census, the next genomes with CUU/CUA at proline positions are hornwort mitochondria: *Phaeoceros* 23, *Anthoceros* 22, 21 and 9, *Nothoceros* 6. Hornworts are known to edit U to C. So U-to-C editing does produce this exact DNA signature, and the census can see it.
- **What argues against editing here:**
  - Hornworts also keep 41–55 CUN at leucine positions; editing is site-by-site.
  - The three *Wickerhamiella* have lost CCN almost entirely: 0–2 among about 4,000 aligned codons per genome. That looks like a code in which CUU/CUA have taken over proline, not like scattered repair of individual sites.
  - No mitochondrial mRNA editing is known in any fungus (section 1).
- **Where editing fits best.** *S. pinicola* is mixed: 36 CCN and 44 CUN at proline positions, with both spellings within the same genes (cob 5/6, nad5 5/6, cox1 7/4, atp6 5/1, cox2 2/7). That fits an ambiguous intermediate stage, but it would also fit editing best. The second origin is therefore the weaker one.
- The post's "Not known" bullet is honest. One sentence on hornworts would preempt the obvious reply (optional fix).

**c. Strain mix-up or contamination: unlikely to matter.**
- **Independent cultures.** Two of the three *Wickerhamiella* have a second culture sequenced by another group: *W. siamensis* PYCC 7069 (Portugal, 2024) and *W. infanticola* DS02 (PacBio, 2017). Both carry the codons.
- **A mix-up can't produce the pattern.** With CUN masked, the tree places each genome where its name says, with strong support. A contaminant or mislabelled culture would land in the wrong place and would not show 99% read agreement.
- **But the second origin rests on one run.** *S. pinicola* (NRRL YB-2263, SRR16988755) and *W. verensis* (PYCC 9704, ERR12422161) each have a single sequencing run. There is no RNA-seq for any of the four species. The post says nothing about *S. pinicola*'s single source (optional fix).

**d. Did Codetta miss it?** There is no evidence that it did. Wolters's decision rule could only pick table 3 or table 4, and our own Codetta run calls proline on these genomes (section 1). The Note's "someone may have seen it and not published" fits that reading. Line 16's "assigned it table 4 after running Codetta" reads as if Codetta saw nothing (wording fix 6).

**e. Two origins.**
- **Robust:** *S. pinicola* sits deep inside 15 *Sugiyamaella* that all use CCN, at 96–100 on both runs, and the trio is at 100. One origin would need about 40 reversions.
- **Not robust, but it doesn't affect the claim:**
  - the backbone: the node joining the trio's region to the rest of *Wickerhamiella*/*Starmerella* has support 10.4 / 54;
  - *Wickerhamiella* is not monophyletic: *W. azyma* sits on a long branch near *Blastobotrys*, and *Starmerella* nests inside;
  - the figure shows that topology without supports (section 4).

**f. A mechanism puzzle worth one line (optional).**
- The *Wickerhamiella* have almost no CCN left, yet all keep tRNA-Pro(UGG).
- Either that tRNA reads CUU/CUA despite a mismatch at the middle base of the codon, or the UAG tRNA is charged with proline.
- The post's "How" bullet already says it is open.

**g. Wolters dN/dS sentence (line 16).**
- The quote is accurate.
- Only *W. infanticola* among their at least seven *Wickerhamiella*/*Starmerella* genomes is affected, so this can explain at most one branch of a clade-level signal.
- "I have not tested how much…" is honest but invites over-reading (optional fix).

## 4. Image

`outbox/2026-11-11-mito-proline-codons/image.png` (1643 × 1597; identical to `results/cun_pro/pro_grid_tree.png`). It shows a tree of 49 mitochondrial genomes and, beside it, a grid of 81 strict-proline positions coloured by codon: CCN blue, CUN orange, other codon grey, no site light grey.

- **Correct.** Rows and counts match `pro_grid_tree.tsv`: *S. pinicola* 34 CCN / 42 CUN / 1 other / 4 no site; trio CCN 0–2, CUN 67–78; all others CUN 0. The tree is the CUN-masked IQ-TREE tree, rooted on *Geotrichum*. It says what the post says: one solid orange band across three sister tips, and one striped row deep inside a blue *Sugiyamaella* block.
- **Honest.**
  - Columns are ordered by gene, not sorted to make bands.
  - The column filter keeps only columns where at least half of the present sites are CCN or CUN. It cannot favour CUN, since 4 rows of 49 can't pass it alone.
  - But the caption omits that filter. It gives only "≥ 60% of these 49 genomes have an aligned codon".
  - The legend says "CUN", but only CUU and CUA occur.
- **Does it read at a glance?** Mostly. The orange band is immediate. Missing:
  - **No title.** The only explanation is a 2-line, 9-pt caption at the bottom.
  - **Labels too small.** At phone width (about 600 px), the 49 species labels are about 4–5 px high and unreadable, and the four genomes aren't marked. Bold or orange labels for the four would let a reader find them.
  - **No support values or scale bar.** That is acceptable for a codon figure. But the long *W. azyma* branch and the *Starmerella*-inside-*Wickerhamiella* nesting draw the eye and raise questions the figure can't answer. Showing supports at just the two key nodes (trio 100; *S. pinicola* in *Sugiyamaella* 96–100) would carry the "two origins" claim.
- **Suggested title:** "CUU and CUA have replaced CCN at proline positions in two branches" or "Proline positions in 49 yeast mitochondria".
- **Suggested legend:** "CUU or CUA (leucine in the standard code)".
- **Suggested caption addition:** "and where at least half of present codons are CCN or CUN".

## 5. Wording fixes

Line numbers are from `cat -n post.md`.

Required:
1. **Line 15, proline novelty (false as written).**
   - Current: "My search and a separate subagent's found no report of any other codon read as proline."
   - Replace with: "No genetic code in any organism or organelle is known to read another codon as proline. The closest cases are partial mistranslation: some *Streptomyces* bacteria carry extra proline tRNAs that sometimes put proline at alanine, asparagine or threonine codons (Vargas-Rodriguez et al. 2021, *PNAS*; Schuntermann et al. 2023, *J. Biol. Chem.*)."
2. **Line 14, Ling 2014 tRNA origin.**
   - Current: "read by a tRNA that most likely arose from tRNA-His"
   - Replace with: "read by a tRNA that most likely arose by duplication of the threonine tRNA (itself descended from tRNA-His)".
3. **Line 5, leucine side.**
   - Current: "CUU and CUA sit where related proteins have proline and never at conserved leucine positions. Leucine is spelled UUA/UUG instead."
   - Replace with: "CUU and CUA sit where related proteins have proline: of the four genomes' CUU/CUA codons at strongly conserved positions, 286 of 291 are at proline positions. Leucine is spelled UUA/UUG, as in their relatives, which barely use CUN at all."
   - The 286/291 figure is from `prepost/cun-pro_work/cun_by_res.out`; re-count it on the census code before using it.
   - Minimal alternative if the new number isn't wanted: keep the sentence and add "(their relatives barely use CUN for leucine either)".
4. **Line 6, scale and the *S. pinicola* denominator.**
   - Current: "44 to 82 CUU/CUA codons at positions where…" → Replace with: "44 to 82 CUU/CUA codons among their 91–93 positions where…".
   - Current: "(*S. pinicola* also still uses CCN at 36 of its 80 proline positions.)" → Replace with: "(*S. pinicola* still uses both: at its 91 such positions, 44 CUU/CUA, 36 CCN and 11 other codons.)"
5. **Line 9, DS02 CCN count (lower bound).**
   - Current: "; a proline CCN codon shows up in 1 of 1,942 reads."
   - Replace with: "." (cut the clause). The per-site list kept only the 3 commonest triplets, so 119 of the 402 mismatching reads were never looked at for CCN. The 65-of-65 result carries the point on its own.
6. **Line 16, Wolters and Codetta.**
   - Current: "and assigned it table 4 after running Codetta, a code-prediction tool."
   - Replace with: "and assigned it table 4. They ran Codetta, a code-prediction tool, but used it only to choose between table 4 and baker's yeast's table 3, so a proline call would not have changed their choice."

Optional:
- **Line 1 (headline):** define CCN at first use, e.g. "…any codon other than CCN (CCU, CCC, CCA, CCG) as proline." CCN is only defined implicitly later.
- **Line 7, plain language:** "A one-sided Fisher test of CUA at proline against its leucine synonyms gives the four p 2.1e-28 to 1.9e-46" → "Per genome, CUA sits at proline positions far more often than the other leucine codons do (one-sided Fisher test, p 1.9e-46 to 2.1e-28 in the four)…".
- **Acronyms and jargon:**
  - "SH-aLRT support 100" → "support 100%". Line 22's "fast run, SH-aLRT" and "1,000 ultrafast bootstraps" can go to the Note.
  - "ENA" → "the European Nucleotide Archive".
  - "PYCC 7069" → "PYCC 7069, Portuguese Yeast Culture Collection".
- **Line 24 (editing), add:** "The only other genomes in the census with CUU/CUA at proline positions are hornwort mitochondria, which are known to edit U to C; but hornworts keep CCN and use CUN for leucine too, while the three *Wickerhamiella* have almost no CCN left."
- **"Not known", add:** "*S. pinicola*, the second origin, rests on one sequencing run of one strain (NRRL YB-2263)."
- **Line 16, Wolters dN/dS:** either cut the last two sentences, or add "only *W. infanticola* among their *Wickerhamiella/Starmerella* genomes is affected". That caveat makes clear the effect is probably small.
- **Line 16:** add "(both from the same sequencing run)" after "identical to BK074532", so no one reads it as an independent confirmation.
- **Line 16:** Sarhan is fine as "December 2025" (online date). If you give a full reference, use *Sci. Data* 13:28 (2026).
- **If Codetta is cited:** cite v2.0 as Shulgina & Eddy 2023, *Bioinformatics* 39:btac802. Wolters cited the 2021 *eLife* paper for v2.0.
- **Image:** add a title; change the legend to "CUU or CUA"; add the "≥ half CCN or CUN" filter to the caption; bold or colour the four species names; show supports at the two key nodes (section 4).
- **Length:** the body is about 28 dense lines plus nested bullets. Candidates to trim:
  - the Wolters dN/dS sentences;
  - the method detail in "A blinded re-check" (keep "a separate agent with its own code called CUU/CUA as proline in exactly these four of eight genomes, and was steered to CUN but not told proline");
  - "Done → Census" second sentence.

Not errors:
- **The "Note for Ryan" is accurate as far as I checked it:** DS02 79.3%, the bootstrap run finishing at 15:48, 18,015 of 18,031, and the PYCC 7069 provenance caveat. I did not check the 83.8% leucine-site figure for DS02.
- **README line 80's *S. pinicola* read tally** is self-consistent (section 2).
