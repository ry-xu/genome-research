# Lab notes (full agent write-ups)

This file joins the four underlying lab reports behind `README.md`:
1. Asgard *Ca.* Lokiarchaeum ossiferum B35 (`outbox/2026-11-14-asgard-protein-intron/post.md`)
2. B35 artefact confirmation tests (`work/archaeal-introns/confirm/RESULT.md`)
3. 3,666-genome archaeal DNA census (`work/archaeal-intron-census/RESULT.md`)
4. TACK RNA scan and *Methanosuratincola* `c0` (`work/tack-rna/RESULT.md`)

---

## 1. Asgard *Ca.* Lokiarchaeum ossiferum B35 (`L02` and `S859`)

**In the cultured Asgard archaeon *Candidatus* Lokiarchaeum ossiferum B35, a 31-base stretch inside a RuBisCO-family protein gene is missing from most of the gene's RNA (385 of 411 molecules). It is cut out at the place and length I predicted from the DNA before opening any read. None of 46 DNA fragments from the culture lacks it. With the stretch removed, two annotated genes join into one protein that about 66 other Asgard genomes encode without a break. Its closest known relative, *L. ossiferum* B-36, found separately and not among those 66, has the gene in one piece and lacks the stretch entirely. A second gene, with no homolog in the other 106 Asgard genomes searched, loses 46 bases the same way. The only endonuclease-cut (BHB) intron in an archaeal protein-coding gene shown by experiment so far is in cbf5 of Crenarchaeota. Asgard genomes carry self-splicing group II introns (a 2024 preprint), but I found no report of an endonuclease-cut intron in an Asgard protein gene.**

Archaea remove introns from their tRNAs with an endonuclease that cuts at a bulge-helix-bulge (BHB) fold. The same machinery removes a short intron from the mRNA of one protein gene, cbf5, in two orders of Crenarchaeota (Watanabe et al. 2002; Yoshinari et al. 2006; Yokobori et al. 2009). Outside cbf5, the reported BHB cases rest on RNA-seq read patterns or sit outside the coding sequence. I was not looking for introns. The lab's reading ledger compares every gene in a genome with known proteins and lists genes that are cut by an extra stretch of DNA. In Asgard genomes it found many such stretches of about 30 bases, in several protein families.

- **The test, fixed before the reads were opened.** B35 has a closed genome (CP104013.1) and public RNA-seq of its enrichment culture. The ledger reported two cuts in it, and I wrote a prediction for each before opening any reads (`work/archaeal-introns/culture/PREREG_CULTURE.md`):
  - L02 (RuBisCO-family gene, BHB motif in its relaxed form at the edges): removed in RNA, at 4487605 + 31 (0-based).
  - L01 (HdrA-like gene, no BHB): not removed. The ledger also calls this one a cut, so it is the control.
- **What the reads show.** The first 10 million read pairs of two RNA runs (SRR27313660, SRR32148620) and 8 million of the culture's DNA (SRR22273278) were used.
  - **L02:** 385 RNA molecules lack the stretch; 383 of them lack exactly the predicted 31 bases. 26 read through it. In DNA: 0 fragments lack it and 46 read through.
  - **L01:** 0 RNA molecules lack it and 24 read through. Both predictions came out as written. (Neither library is stranded, so by the rule every molecule counts. The script I registered handles only the stranded case, so I applied its thresholds to all molecules with a second script, written after the count. On message-strand reads alone L01 has too few reads for a word.)
  - **Controls:** at places 200 bases to either side, 0 molecules in RNA or DNA lack a stretch. The RNA reads match the closed genome at 99.6–99.9 % identity.
  - **A second reader:** a separate agent recounted all three places from the raw reads with its own code, without seeing my counts. At L02 it found 393 molecules without the stretch and 24 through (392 at the predicted place), and 0 without it in DNA.
- **What the splice makes.** The stretch holds the stop codon of the first annotated gene (UYP47557, annotated as a RuBisCO-like enolase). Removing it joins that gene in frame to the next one, UYP47558 (annotated ribulose bisphosphate carboxylase), 22 codons before that gene's annotated start, giving one 441-amino-acid protein. The cut falls exactly between two codons.
- **The relatives have no break.** I searched 107 Asgard genomes with the spliced protein. In about 66 of them, a copy of this gene spans the junction in one piece, with no stop and no stretch. (In my first count of 88, an independent check found that about 22 were closer to a distant B35 paralog than to this gene.) Four close relatives, extracted and translated independently, have the same residues across the junction as B35's spliced protein (…RPLLNNMI|KPCTGH…). The 31-base stretch is found in no listed genome except B35.
  - The closest is a genome NCBI names as another strain of the species, *Candidatus* Lokiarchaeum ossiferum B-36 (GCA_041410205.1), which is not among the 107. Whole-genome identity is only about 88 %, so it may be a sister species. Its copy (XEO79081.1) is one uninterrupted gene, 95 % identical to B35's spliced protein, and the stretch is nowhere in its genome. So the stretch most likely arrived in B35's line after the two split.
- **A second gene.** Outside the registered test, the commonest removal in the whole RNA set is 46 bases at 859933 (0-based). Removing them joins two annotated "hypothetical protein" genes (UYP44411 and UYP44410) into one 235-amino-acid frame. I counted 856 molecules without it; the second reader counted 822 without it and 74 through. DNA: 0 without it. This gene has no homolog in the other 106 Asgard genomes at the threshold used (DIAMOND, e ≤ 1e-5), so conservation cannot test it. The ledger missed it, because it only sees genes with a known relative. Of the removals checked in this culture's RNA, this one has the strongest motif-and-fold support for the archaeal endonuclease (another agent's check; only this site, L02 and the tRNA-SeC intron were folded). It passes the BHB motif in its strict form, and its predicted lowest-energy fold has all four central-helix pairs and both 3-base bulges. The central helix forms with probability 0.71–0.76 at 37 °C, a level random 46-base stretches reach with p about 0.002 once the choice between two possible start positions is allowed for. At 20 °C it is 0.91. Both are more than the culture's tRNA-SeC intron reaches (0.43 at 37 °C, 0.34 at 20 °C).
- **Circles.** Spliced-out archaeal BHB introns are typically ligated into circles. Another agent (Altair) looked for reads across a circle junction in the first 10 million spots of two runs and found none at either gene. On 19 Nov I searched all 65.7 million spots of two libraries rich in L02 (469A, 469B; the A202 pair has as many L02 reads) with the same probes. L02 has 1,917 joined reads and 0 circle reads; the 46-base site has 5,044 joined reads and 0 circle reads. The only intron circles the probes found are 3 reads of the 54-base tRNA-Tyr intron, as tandem repeats. The culture's 30-base tRNA-SeC intron, removed in 52 reads, also shows no circle. So no circle of a ~30-base intron was seen at any site probed, and the missing L02 circle counts neither way. Why the 46-base site shows none is not known. (My rule fixed beforehand called the 0 "counts against". Its expectation divides by mature tRNA reads, of which the probes found only 3 because ARAGORN's intron edges are approximate. With the 52 exact-edge tRNA-SeC reads added it still says "counts against", and it gives the same verdict for the tRNA-SeC intron itself.)
- **Is it made while the library was made? (19 Nov, rules fixed before counting: `work/archaeal-introns/confirm/`.)** All eight B35 RNA libraries are Illumina from one lab; none is RT-free, and no later B35 DNA is public. Inside two of these libraries (first 10 million spots of A202_A and 469B):
  - **Not by PCR.** 100 groups of duplicate reads at L02 (same run, same ends, so mostly copies of one molecule) hold 2 or more reads each. If the stretch were skipped while copying, about 13.5 groups would hold both forms; none does (P 1.4 × 10⁻⁶). This does not test a skip in the first, RT, copy.
  - **Not a general skip over hairpins.** As DNA, L02's stretch folds weakly (−3.6 kcal/mol). Of 3,735 places in these reads that fold at least as well, none loses 10 % of its reads (3 removal reads against 1.2 million read-throughs, summed over places). As RNA, Altair's scan gave the same (below).
  - **Protein.** In a B35 proteome from the Schleper group (PXD074767, measured May 2024, released 2026; 15 runs, searched against the annotated genes), UYP47557 is seen through 5 peptides in 9–15 runs each, and its own unspliced C-terminal peptide in none. That is weak support for the spliced form; the junction peptide was not in their database.
  - **No second copy** of the stretch or its hairpin exists in B35 or B-36, so whether other copies are removed cannot be tested.
  - **An independent attempt to kill it failed.** A fresh agent, given only the raw data, counted the junction directly in the reads without the aligner (0 of 51 DNA reads lack the stretch; 285 of 316 and 342 of 349 RNA reads do). It also found the share the same whichever strand read 1 came from, no B-36 reads and no strain difference fixed in both RNA cultures, and a third library (469A) at 99 %. Its verdict: not killed. The routes it could not exclude are a B35 line that lost the 31 bases after the 2022 DNA sample and took over the culture, and a reverse-transcriptase skip caused by something unique to this sequence. Experiments (a) and (e) below would settle both.

![](image.png)

- **Known before:**
  - **cbf5.** Watanabe et al. 2002 (FEBS Lett) showed splicing in vivo by RT-PCR. Yoshinari et al. 2006 (BBRC) detected the protein and showed the cut made in vitro by a lab-made endonuclease. Yokobori et al. 2009 (BMC Evol Biol) mapped the intron across Desulfurococcales and Sulfolobales.
  - **Candidates from RNA-seq.** Berkemer et al. 2015 report circles from introns inside reading frames in *Methanopyrus kandleri*. Danan et al. 2012 found circles inside genes of *Sulfolobus solfataricus* but no spliced mRNA. In *Haloferax volcanii* HVO_1309 the intron is in the 5′ UTR (Schwarz et al. 2020).
  - **Asgard.** Neither the B35 paper (Rodrigues-Oliveira et al. 2023) nor the MK-D1 paper (Imachi et al. 2020) mentions introns or splicing in its main text. Mattick et al. 2024 (bioRxiv) report group II introns in Asgard genomes, from sequence analysis. Group II introns are a different class, which splice themselves out as RNA rather than being cut by the tRNA endonuclease. A 2026 preprint (Goswami et al.) finds many hidden tRNA introns in Asgard archaea; it does not cover protein genes.
- **Not known:**
  - **What cuts it.** L02 passes the BHB motif only in its relaxed form, and its predicted fold does not form the bulge-helix-bulge (another agent's check). Its lowest-energy fold forms none of the four central-helix pairs at any of the three possible starts (37 °C). At the best of the three starts the central helix forms with probability about 0.09, at both 37 °C and 20 °C; at the other two it is 0. 1–2 % of random 31-base stretches reach that level, about 4 % once the choice among three starts is allowed for. The second gene's site does fold as a BHB (above). The endonuclease is assumed at both sites, not shown.
  - **Whether the joined protein is made: the authors' proteome cannot say.** I searched all 8 fractions of the authors' B35 proteome (PRIDE PXD059733; 424,558 MS2 spectra) with Sage, by a rule fixed before the search (`work/archaeal-introns/protein/PREREG_PROTEIN.md`; four amendments. The last changed only the file gate. It came after a test run on 57 % of one fraction, during which I saw 10 PSM rows from this locus; none was a junction or unspliced peptide. Protein from this locus is there: 12 peptides, all present in both forms. Seven are from UYP47557 1-166 (6 if FGIEGIR is set aside, see below) and 5 from UYP47558. Neither form's own peptide was identified: no PSM for the one tryptic peptide over the junction, and none for UYP47557's own C-terminus. No MS2 spectrum was taken of the junction peptide's precursor beyond chance. One spectrum fits it, only at isotope error 2, at 0.88 min in the void volume. The gradient's identified peptides elute from 14.6 min, and those of 30 aa or more at 75-99 min. By the rule this reads "protein detected, junction not seen", which the rule calls weak. These scans stop at 2000 m/z, and the search counted only singly charged fragments. So only 10 of the 33 junction-spanning ions (b9-b18) could have been matched. With no precursor sampled, though, a search allowing doubly charged fragments would have had nothing to score. The sample is the enrichment culture, and the partner genomes are not public, so the database is B35 only; about 9 % of spectra are identified.
  - **Which strand the RNA is.** Neither RNA library is stranded. The spliced reads keep the gene's reading frame, but their strand is not shown.
  - **The other cuts.** The ledger lists many more ~30-base stretches in Asgard and some other archaeal genomes. In a wetland soil dataset the RNA was too thin to answer by the rule I fixed beforehand ("too few"). Outside that rule, 6 message-strand molecules at 4 sites lacked exactly the predicted stretch.
- **For Ryan:**
  - **Checked by someone other than me:** (1) a blind second reader of the culture reads (`agents/vega/kills/culture_recount/REPORT.md`); (2) an independent agent that tried to kill the conservation claim (`agents/vega/kills/intron_conservation_check/REPORT.md`). It corrected the 88 to about 66, by a best-hit split rather than a tree (10 of the 66 and 12 of the 22 are ambiguous), and found B-36. (3) Altair's circle test (`agents/altair/kills/culture_circles/`); (4) a prior-art search (`work/archaeal-introns/culture/PRIOR_ART.md`). (5) This text was fact-checked sentence by sentence by a fresh agent (`agents/vega/kills/outbox_intron_factcheck/REPORT.md`, 14 Nov 01:15); its corrections are applied here, and it checked the prior-art citations. (6) Altair's strain and deletion check (`agents/altair/kills/vega_strain/RESULT.md`), which gave risks (1) and (2) their numbers (two of its eight libraries are counted from my alignments). A second fresh agent fact-checked risks (1) and (2) (`agents/vega/kills/outbox_strain_factcheck/REPORT.md`, 14 Nov ~04:05); its corrections are applied. (8) The protein readout was rechecked with separate code by a fresh agent (`agents/vega/kills/protein_readout_recheck/REPORT.md`). Every count agreed. It found two statements the rule requires missing from the read-out: the "weak" label and the observability limit. Both are now added. It also found the one-spectrum sampling figure, which I recounted with my own script (`vega_recount.py` in the same folder). A second fresh agent fact-checked the protein sentences (`agents/vega/kills/outbox_protein_factcheck/REPORT.md`) and recounted that figure independently; its corrections are applied, but the corrected wording itself was not rechecked. One of the 7 peptides, FGIEGIR, also occurs in *Methanoculleus* RuBisCOs in UniProt (one species was formerly placed in *Methanogenium*, a partner's genus), so a partner could have made it. 6 remain with no UniProt hit in a culture partner's genus (still ≥ 5), but partner proteomes are not public. Not rechecked: the group II definition.
  - **Also checked by someone other than me (7):** Altair's motif and fold check (`agents/altair/kills/bhb_fold/`, with a fresh-subagent audit in `audit/REPORT.md` that covered the motif, fold and hairpin-scan numbers but not the survey of other removals). It reproduced the motif call with separate code and computed the folds in "What cuts it" and "A second gene". It also tested the other removals of 15–200 bases in the culture RNA (170 distinct sites; longer removals were not tested). Outside repeats (sites with at most 5 equivalent start positions), the only ones passing the strict motif are a tRNA-SeC intron (3170068 + 30) and one site seen in a single read (1428226 + 100). Of the 72 removals of 51 bases or less, the only one passing the motif (relaxed or strict) is the tRNA-SeC intron, against 3.8 % of random positions. So a relaxed pass is not typical of this RNA's short removals. But 3.8 % of random positions pass anyway, and L02 passes only through its inner helix, which is also what a skipped hairpin's own stem would give. On reverse transcriptase skipping a hairpin: no locus with a strong hairpin (ΔG ≤ −7) loses more than about 1 % of its reads, against 94 % (385 of 411 molecules) at L02. That scan needs 30 reads through a locus, so it cannot include L02 itself. It tests hairpins in general, not L02's own near-palindrome flanked by a TTAA repeat. Short removals are enriched at hairpins in this RNA (42 % against 16 % of random sites), so a low-rate skip artefact does occur here.
  - **Not checked by anyone else:** Altair's survey of the other removals (Altair only, after the audit); removals longer than 200 bases; any folding program other than ViennaRNA, and whether the RNA folds that way in the cell; whether L02's own sequence could be a reverse-transcriptase skip or template switch (the 19 Nov tests argue strongly against a PCR-stage skip and a general hairpin skip, not against an event specific to this sequence at the RT step; an independent agent's attempt, verdict "not killed", is in `work/archaeal-introns/confirm/INDEPENDENT_agent.md`; the new 19 Nov text was fact-checked by a fresh agent, `agents/vega/kills/outbox_confirm19_factcheck/REPORT.md`, and its corrections applied, but the 'independent attempt' sub-bullet came after that check and was not rechecked). The fold and check-(7) text (edited 14 Nov about 06:59) was fact-checked by a fresh agent (`agents/vega/kills/outbox_bhb_factcheck/REPORT.md`) and its corrections are applied; the corrected wording itself was not rechecked.
  - **Risks:**
    - (1) The culture is an enrichment, so the RNA could come from another member. Against that: identity to B35's closed genome is 99.6–99.9 %, and no DNA fragment lacks the stretch. But the DNA comes from a different sample than either RNA run, and SRR32148620's sample record names no strain, although the run sits in the B-35 genome BioProject. Altair then counted fixed differences from B35's genome (depth ≥ 10, one allele ≥ 80 %): 0 in 5,992,344 DNA positions; 1 in the 1,648,851 positions SRR32148620 covers (one further insertion call is a counting artefact); and in SRR27313660, which covers 757,706 positions, 99, all in five clusters that may be another organism's reads mapping there (not checked). So over the 1.65 Mb and 0.76 Mb of the 6.04 Mb genome that the RNA covers, it shows no strain-level difference from B35 outside those five unexplained clusters.
    - (2) A deletion that arose in the culture after the DNA was taken would also differ by about 0 SNPs. B35's BioSample gives a collection date of 2019-04-21, probably the sediment sampling; when its DNA was extracted is not stated. SRR32148620's sample is dated December 2024. Altair found no later DNA of this culture among 18 Lokiarchaeum runs in SRA; one later DNA project (PRJNA1091558) is the B-36 culture (its BioSample names strain B-36; checked 19 Nov). Against a deletion: if both genes are expressed equally per cell, one lineage carrying both deletions would give the same unspliced share at both sites within a library. Altair counted all eight RNA libraries, mapping the six I did not use (first 10 million spots each, as for mine), and the shares differ. In the 469 pair, pooled, 9 of 582 reads (1.5 %) read through at L02 against 96 of 1,448 (6.6 %) at the 859933 site (Fisher p 3.7e-07; 2.3e-04 with every count halved as a rough bound for read pairs). The 478 pair points the same way (0 of 53 against 8.2 %, p 0.021). The A202 pair leans the other way (4.9 % against 3.0 %, p 0.071), all of it from one library. That fits precursor RNA and argues against one lineage carrying both deletions. Not excluded: two separate deletion lineages, or a deletion that changes its own gene's expression.
    - (3) Novelty rests on searches of abstracts and open texts. Not opened: the B35 supplements, Berkemer et al.'s gene list, and Becker 2017 (blocked). A 2025 review (Liang et al., MMBR) reportedly still names cbf5 as the case, but nobody here could open its full text, so it is left out of the post.
    - (4) One intron in one strain is a curiosity until the machinery or the protein is shown.
  - **Recommend:** the proteome test is done and does not decide the question (above). A short post with the read figure, stated as an observation, with the protein result given as the pre-registered reading, "protein detected, junction not seen" (weak). In the authors' DDA proteome (PXD059733), the junction peptide's precursor was not fragmented beyond chance. The direct test would be a proteome with the junction peptide targeted, or made with a different protease, with multiply charged fragments allowed.
  - **Experiments that would decide it (for the collaboration, 19 Nov):** (a) RT-PCR across L02 with two different reverse transcriptases, read out by amplicon sequencing, with an in vitro-transcribed unspliced L02 RNA as template control (an RT skip would appear on it too) and a no-RT control; (b) RT-free evidence: nanopore direct RNA after poly(A) tailing (archaeal mRNA has no tail), RNase protection across the junction, or a northern with a probe inside the 31 bases (it sees unspliced RNA and the excised piece; a probe across L02 cannot resolve 31 bases); (c) detecting the excised piece itself, linear or as an RNase R-resistant circle; (d) targeted proteomics (PRM) with heavy standards for the junction peptide RPLLNNMIKPCTGHTAAVAGDLAYQVALGGVDVVK (PLLNNMIK…, if trypsin also cuts before proline) and the unspliced RPLLNNMIN, or a second protease; not finding it would not refute splicing; (e) DNA from the same culture batch as the RNA, against a deletion lineage.
  - **Would care:** the B35 culture's corresponding authors (Schleper lab, University of Vienna; Pilhofer lab, ETH Zürich); groups working on archaeal tRNA splicing endonucleases.
  - **The image** (made 15 Nov, `work/archaeal-introns/culture/src/fig_l02.py`) draws the same molecules the text counts: the 411 RNA molecules from `out/rna_mol_templates.tsv` and the 46 DNA fragments re-read from the saved DNA alignment with the same rules. The script stops unless the counts are 385 / 26 / 0 / 46. A fresh agent checked the image and recounted every number from the alignments (`work/archaeal-introns/culture/FACTCHECK_fig_l02_agent.md`); its fixes are applied. Two details from it: (a) of the 385, 383 have a read with exactly the predicted removal and 382 have it in every read; in the other 2 the aligner placed the gap 4 bases away, on a mismatch tie or one low-quality base. (b) The molecule key includes read orientation, which means nothing in these unstranded libraries; merging orientations gives 364 / 24 molecules, the same 94 %. The DNA lines are fragments (46; 34 distinct outer ends). The corrected figure and this bullet were not rechecked after those fixes. My own wrong "correction" of 383 to 382 earlier on 15 Nov was caught by that check and undone; 383 stands.
  - Nothing has been sent.

---

## 2. B35 artefact confirmation tests (`confirm/RESULT.md`)

# Asgard B35 L02 31-nt RNA cut: kill-or-confirm results (Vega, 15:24 19 Nov 2026)

Rules frozen in `RULES.md` (md5 397767c644613bcd109524021033d775 before the post-freeze notes were appended) before
any of these data were read. Every number below is printed by the named script into the named .out file.

| Test | Question | Result | Frozen word |
|---|---|---|---|
| C1 `c1.py` → c1.out | Skip during PCR? Duplicate families (same run, same outer ends) should be mixed if each copy skips independently | 391 families, 100 with ≥2 reads; removal share p = 497/529 = 0.94; expected mixed 13.47, observed 0 (Poisson P 1.4e-6) | **PCR-stage skip rejected** |
| C2 `c2.py`, `c2_read.py` → c2.log, c2.out | Polymerase skip over a DNA hairpin (2nd strand/PCR)? Altair's hairpin scan rerun with DNA folding | L02's 31-nt window folds weakly as DNA (−3.6; 21 % of all loci are at least as stable). 3,735 loci at least as stable (≥30 reads through, L02/S859/tRNA introns excluded): 0 have a removal share ≥ 0.10; 2 have any removal read (3 reads in 1,195,616 through). At ≤ −7: 154 loci, 0 removal reads in 45,416. Detector check passed: L02 563 removal / 24 through, S859 window 1,285 / 95 | **not general** |
| C3 `c3/run.sh`, `c3/count.py` → c3/c3.out | Excised piece as a circle? All 65.7 M spots of 469A+469B, Altair's probes unchanged | L02 1,917 joined, 34 through, **0 circle reads**, 0 in all shifted controls. S859 5,044 joined, 0 circle. Control: 3 tRNA-intron circle reads (all tRNA-Tyr, 54-nt intron, rolling-circle tandem repeats), only 3 mature-joined tRNA reads → E = 3/3 × 1,917 | **"Missing circle counts against"** — but see below |
| C4 `c4/c4_read.py` → c4/c4_x767.out, c4/c4_x696.out | Unspliced protein made? Two later B35 DIA proteomes, vendor tables | PXD074767 (15 runs, peptide-level Spectronaut table): UYP47557 has 5 peptides in 9–15 runs each; the unspliced C-terminal peptide (…RPLLNNMIN) is in none. PXD047696: protein-level tables only | 074767: **Not seen, weak**; 047696: **No information** |
| C5 → c5/c5.out | Other copies of the hairpin removed too? | No second copy of the 31-mer or the ±20 nt hairpin in B35 or B-36 (only 15–16 nt chance matches) | nothing to test |
| C6 → `INDEPENDENT_agent.md`, `indep/` | Independent kill attempt (fresh agent, raw inputs only, own controls) | Aligner-free k-mer count 0/51 DNA, 285/31 and 342/7 RNA; PCR skip excluded in SRR27313660 (1 mixed vs 8.29, P 0.0023); share equal by read-1 orientation (0.907/0.931, 0.992/0.993); no B-36 reads, no SNV fixed in both RNA cultures (385 kb); 0/5,401 matched hairpins reach share 0.01 (L02 0.946); third library 469A raw share 0.990; circles not testable | **not killed**; strongest open route: a B35 subline that lost the 31 nt after 2022; runner-up an RT skip unique to L02 |

## C3: why its frozen word should not be read as evidence (post hoc, stated as such)
The rule's expectation E divides tRNA circle reads by mature (joined) tRNA reads. Two things make that ratio
meaningless in these libraries, both found only after counting:
1. The denominator is a probe-edge artefact: 3 joined reads against 873 pre-tRNA (through) reads over 14 ARAGORN
   intron-tRNAs in 65.7 M spots, because the joined probes use ARAGORN's approximate intron edges (it gives tRNA-SeC a
   28-nt intron at 3170071; the RNA shows 30 nt at 3170068). With tRNA-SeC's 52 exact-edge joined reads added,
   E = 3/55 x 1,917 = 104.6 (fact-checker's recount), still "counts against".
2. The same computation on the culture's own tRNA-SeC intron (3170068 + 30 nt, the exact removal seen in RNA, strict
   BHB motif) gives E = 52 and 0 circle reads, so the rule would "count against" a tRNA intron too. The only circles
   seen are of one 54-nt intron, as rolling-circle repeats; none of a ~30-nt intron was seen at any site probed (only the probed sites were searched). Why the 46-nt S859 site
   (5,044 joined) shows none is not known; size alone does not explain it when a 54-nt circle is seen.
So 0 circles for a 31-nt piece is what these libraries also give for the 30-nt tRNA-SeC intron (one comparator; assumes it circularizes in the cell). C3 neither supports nor
counts against the cut. A test that could: RT-PCR or RNase R-resistant small-RNA sequencing for the 31-nt circle.

## Overall (by the reading frozen in RULES.md)
- Library-artefact alternative "killed as far as public data can" requires C1 rejected, C2 not general, C5 no copy
  removed: C1 rejected, C2 not general; C5 has no copy to test (cannot fail, adds nothing).
- C3's frozen word points against the cut, but the identical calculation fails on the culture's own tRNA intron
  (above). C4 adds a weak point for the spliced form.
- What public data cannot exclude: a first-strand (RT) event specific to L02's own sequence that no other locus
  shares (C2 tests hairpins in general, not this one), and a deletion lineage in the culture (Altair's 14 Nov share
  test argues against one lineage for both sites). These need new experiments:
  1. RT-free RNA evidence: nanopore direct-RNA of B35, or RNase protection/northern across L02.
  2. RT-PCR across L02 with two different RTs (e.g. a group II intron RT and an MMLV RT) at two temperatures, plus
     an in vitro-transcribed unspliced L02 RNA as the template control: an RT skip must appear on the synthetic
     template too.
  3. Targeted proteomics for the junction peptide RPLLNNMIKPCTGHTAAVAGDLAYQVALGGVDVVK (PRM), or a second protease.
  4. Fresh B35 DNA from the same culture batch as the RNA (rules out a deletion lineage).

## Corrections from the fact-check (15:38; `agents/vega/kills/outbox_confirm19_factcheck/REPORT.md`)
- RULES.md C3 calls 469A/469B "the two libraries with the most L02 reads"; wrong: in the first 10 M spots A202_B has
  the most (338; 469B 309, 469A 273, A202_A 274). Does not affect the test.
- C1 and C2 used the first 10 M spots of A202_A (SRR27313660) and 469B (SRR32148620) only. Duplicate groups are
  mostly, not all, copies of one molecule (9 of 385 end pairs recur across two runs). C1 holds with per-run p
  (E 12.67), orientation kept (E 11.89), and a simulated PCR genealogy (E 12.8-15.5).
- C2's 1,195,616 is read-throughs summed over loci. Loci need >= 30 reads through (L02 has 24); the fact-checker's
  direct scan of 15-51 nt removals at sites with >= 3 reads found share >= 0.10 only at L02, S859, tRNA-SeC and
  1291865 (foreign reads stitched by the aligner, MAPQ 23).
- PXD074767 was measured May 2024 (released 2026), Schleper group.

## Verdict (15:50 19 Nov)
**Not killed; not confirmed.** Every library-preparation or mapping route that public data can test fails to
reproduce a 92-99 % removal at L02 (C1, C2, C6 T1-T3, T5), in two submissions and three libraries checked at full
depth. Two routes stay open and need the lab: (1) a B35 subline that precisely lost the 31 nt after the 2022 DNA
sample (against it: no SNV fixed in both RNA cultures, and S859 removed at a different share within the same
library); (2) an RT skip driven by something unique to L02's sequence. Experiments (a) RT-PCR with two RTs plus an
in vitro-transcribed unspliced template and a no-RT control, and (e) DNA from the same cultures as the RNA, would
settle both. The protein side is weakly in favour (C4) and needs targeted PRM for the junction peptide.

---

## 3. 3,666-genome archaeal DNA census (`archaeal-intron-census/RESULT.md`)

# RESULT: BHB-shaped insertions in protein genes are common in Asgard archaea, but the rate is set by class, not by being Asgard (Mintaka, 19 Nov 2026)

The frozen DNA verdict is **met**: 373 of 641 placed Asgard cuts carry the BHB motif, against 12.5 expected from
decoys, and all 11 testable Asgard classes pass. Several TACK classes carry the same kind of insertion at the same or
higher rates, and the Euryarchaeota side shows no excess over decoys. Asgard as a whole has about twice TACK's rate
once the shared cbf5 site is set aside, but Asgard and TACK classes interleave from about 0.1 to 2 per 1,000 genes
read. The RNA half of the claim is **not met and could not be tested**: the two Asgard cultures with public RNA carry
no BHB cut. The one Asgard protein-gene insertion known to be removed in RNA (Vega's L02) does not fold as a BHB in
Altair's check, and its mechanism is unknown. So removal by the tRNA-splicing machinery is not shown for any Asgard
protein gene. In TACK, only the known cbf5 introns have RNA coverage, and they are spliced.

Rules: `RULES.md`, frozen at 15:26 before any cut or control score was read (md5 477813327c76c4636097d68fd64b051f,
posted in commons). Amendments 1–7 are in the file, each with its reason; the final md5 is
afc423b558f6f14fa0ff9041b4d62e0f, and every amendment md5 is in `notes/md5.txt`. Code and write-up were audited three
times by fresh subagents (`notes/AUDIT_1.md`, `notes/AUDIT_2.md`, `notes/AUDIT_3.md`). Audits 2 and 3 recomputed every class
row, CI, stratum, R7 site and K1 row from the raw cuts and matched. Audit 3 also checked this file's numbers and wording,
and its fixes are applied below. Commands are in `README.md`; full output is in `out/summary/summary.txt`.

## What was counted
The units are GTDB archaea in the atlas cache plus the GTDB Asgard species reps: 3,666 genomes after dropping 16
GCA/GCF twins (amendment 5). Each genome was read with the reading ledger v0.7 against the anchor protein set. A
**cut** is a stretch of 15–300 nt that the ledger must skip to join two in-frame pieces of one protein (R1). A cut is
**placed** when an in-frame removal beats both the frameshift and the stop-read-through models (R2). Each placed cut
gets a BHB label (R3, Vega's checker, unchanged) at every position its microhomology allows (R4). It also gets the
same label on same-shape decoys, the same length moved along the same gene (R5). Expected counts use the decoy rate
(amendment 4). The main tier is used for the verdict; the wide tier is reported beside it.

These are **DNA candidates for introns, not proven introns**. Recall on 25–40 nt plants is 13/76 (C4), and counts
are not scaled up.

## Controls (all pass)
| control | bar | result |
|---|---|---|
| C1, scorer on tRNAscan-verified tRNA introns | BHB share ≥ 0.60, decoys ≤ 0.08 | 0.614 (n 11,009) vs 0.005 |
| C2, known cbf5 introns found blind | placed BHB TRUB cut in ≥ 3 of 4 species | 3/4 on the main tier (*A. pernix* L38, *S. solfataricus* rep L22, *M. sedula* L19); *S. acidocaldarius* L22 wide tier only |
| C3, Vega's L02 (*L. ossiferum* B35 RBLL, removed in RNA) | cut within 10 nt of 4487605, L 31 ± 3 | placed at 4487601, L 31, BHB, main tier |

C3 tests placement only. R3 labels L02 BHB, but Altair's fold check (cited by Vega, commons 15:50) finds that L02
does not fold as a BHB. So R3's sequence label and a fold call can disagree on a real removed insertion.

## Verdict (frozen bar: Asgard obs ≥ 20, obs/exp ≥ 3, ≥ 2 Asgard classes passing)
**Met.** On the main tier, 373 of 641 placed Asgard cuts are BHB-labelled against an expected 12.5 (ratio 29.9).
All 11 Asgard classes with ≥ 5 placed cuts pass the Bonferroni class test (m 40). The weakest are Thorarchaeia
(p_bonf 0.005), Heimdallarchaeia (0.024) and Njordarchaeia (0.032). Lokiarchaeia: 98 of 161 genomes carry at least one BHB cut, across 50 protein families.

- **MAG quality does not explain it (amendment 6).** Asgard obs/exp is 36.5 in compl ≥ 90 / contam ≤ 5 genomes (117/184 placed), 27.4 in the rest, and 35.0 in genomes of ≤ 10 contigs.
- **MAGs versus RefSeq genomes (post hoc, `src/src_strat.py`).** In TACK + Asgard, non-cbf5 BHB cuts are rarer in RefSeq (GCF) genomes than in MAGs of the same class (GCF obs 4 vs 19.4 expected from GCA rates per locus). That gap is lineage mix: cultured genomes come from low-rate orders and families. Matched by family, it is 4 vs 4.8 (P 0.47), and matched by genus, 1 vs 4.1 (P 0.09). Only 108 GCF genomes take part, so a twofold MAG excess cannot be excluded. Asgard has one GCF genome (no cut), so Asgard rests on MAGs and cultures deposited as GCA.
- **Not one family or one genome.** No genome carries more than 9 BHB cuts (*Lokiarchaeia* GCA_004524535.1). TRUB (cbf5) accounts for BHB cuts in 19 of the 173 Asgard genomes with any.

## Asgard compared with other archaea (post hoc, descriptive)
The class tests are frozen; this grouping is not. It is built by `src/group_compare.py` (output `out/group_compare.txt`).
Asgard = Asgardarchaeota; TACK = Thermoproteota + Korarchaeota; DPANN = Nanobdellota, Micrarchaeota, Aenigmatarchaeota,
Altiarchaeota, Undinarchaeota, Huberarchaeota, Nanohalarchaeota and Iainarchaeota; Euryarchaeota side = all other
phyla. Main tier, k ≤ 3, kept, placed cuts.

| group | genomes | BHB cuts (TRUB) | obs/exp vs decoys | per 1,000 loci | same, TRUB removed | genomes with ≥ 1 | families with a BHB cut |
|---|---|---|---|---|---|---|---|
| Asgard | 400 | 373 (19) | 29.9 | 0.74 | 0.70 | 43 % | 84 |
| TACK | 896 | 517 (250) | 26.7 | 0.64 | 0.33 | 41 % | 37 |
| DPANN | 1,068 | 74 (21) | 10.7 | 0.16 | 0.11 | 6 % | 13 |
| Euryarchaeota side | 1,302 | 32 (0) | 0.93 | 0.02 | 0.02 | 2 % | 23 |

- **cbf5 makes TACK look like Asgard.** Half of TACK's BHB cuts (250 of 517) are the one TRUB site, against 19 of 373
  in Asgard. The Asgard/TACK per-locus ratio is 1.16 with everything counted (genome bootstrap 95 % 0.97–1.38), 2.14
  with TRUB removed (1.72–2.58), 1.50 with Bathyarchaeia removed (1.22–1.82) and 2.93 with both removed (2.28–3.77).
- **But the rate belongs to classes, not to Asgard.** With TRUB removed, BHB cuts per 1,000 loci (classes with ≥ 8 genomes):
  Jordiarchaeia 2.12, CALKCG01 (TACK) 2.09, Hermodarchaeia 1.76, Methanomethylicia (TACK) 1.14, Lokiarchaeia 0.86,
  EX4484-205 (TACK) 0.85, Korarchaeia 0.62, Bathyarchaeia 0.58 … Thorarchaeia 0.11, Heimdallarchaeia 0.09,
  Thermoprotei_A 0.06. Two large Asgard classes sit below the DPANN Aenigmatarchaeia (0.24).
- **Classes that pass besides Asgard (11):** Bathyarchaeia (148 of 225 genomes, 20 families, p_bonf 4e-190),
  Thermoprotei_A, Nanobdellia, Methanomethylicia, Nitrososphaeria, Thermoprotei, CALKCG01, Nitrososphaeria_A,
  Aenigmatarchaeia, Korarchaeia and EX4484-205.
- **Classes with p_bonf = 1 (18):** Halobacteria, Methanosarcinia, Micrarchaeia, Thermoplasmata, Methanobacteria,
  Methanococci, Archaeoglobi, Methanomicrobia, Iainarchaeia, E2, EX4484-52, Bog-38, Syntropharchaeia, Thermococci,
  Nanosalinia, Penumbrarchaeia, Poseidoniia and SW-10-69-26. No tested class falls between 0.05 and 1.

So "Asgard archaea carry BHB-shaped insertions in protein genes" holds. Asgard as a group has the highest non-cbf5
rate and spreads it over the most families. But the rate varies 20-fold between Asgard classes, and several TACK
classes match the richest Asgard ones. So "more than other archaea" holds only at group level and not for every
class. The Euryarchaeota side shows no excess over decoys at all.

## Sites shared across phyla
Sites were clustered on reference-protein columns (R7). The output has 61 sites, 30 of them in ≥ 2 classes. R7 splits
one site into two rows when the microhomology slide changes the phase, and nested windows give separate rows. Merging
both (Audit 3) leaves **50 distinct sites, 21 in ≥ 2 classes**.
- **TRUB (cbf5), columns 80–81:** 287 genomes (236 at phase 1, 15 classes, plus 51 at phase 2, 8 classes). It spans
  Thermoproteota (196 at phase 1), Nanobdellota (19), Asgard (18, none in Lokiarchaeia) and Korarchaeota (3). This is
  the known crenarchaeal cbf5 intron position (Watanabe et al. 2002).
- **PIMT column 70–71:** 86 genomes in 7 classes at phase 1 (Asgard, Korarchaeota, Thermoproteota).
- **IF2P (aIF2-γ) columns 44/45:** one site split by phase (nt 134/135); rows of 20 and 17 genomes, 6 classes each.
- **Mostly one class:** DNAG 43 is 32 Lokiarchaeia genomes plus single Jordiarchaeia and Nitrososphaeria_A genomes,
  with near-identical intron sequence among the Lokiarchaeia (alike12 0.94). Y1416 151 is 22 genomes, alike12 0.91.
- **What this does and does not show.** A shared column across phyla fits an old intron kept in place. It equally fits
  insertion at a preferred proto-site, or gene transfer. Intron sequence does not settle it: no Asgard TRUB intron
  shares a 12-mer with another phylum's, and Asgard TRUB introns are longer (median 31 nt vs 19–23).

## The RNA test: BHB-shaped is not the same as spliced
- **K1 survey.** In the only census genomes with RNA (the 81-genome survey), 6 BHB cuts are covered. Two are removed in RNA, and both are the known cbf5 introns in Thermoproteota (*A. pernix* 28 templates; *S. acidocaldarius* 25). The RNA reads straight through the other four, all in Euryarchaeota-side genomes: *Halorubrum lacusprofundi* Y1395, *Methanobrevibacter ruminantium* RLSX, *Methanothermobacter marburgensis* ACDA, and *Methanothermococcus* Y1642 (k 24, outside the k ≤ 3 statistic). None of the 23 covered non-BHB cuts is removed.
- **What that means for the census.** The read-through BHB cuts sit in classes whose BHB rate matches the decoys (p_bonf = 1), where a BHB label is expected by chance. So they are the background the decoys model, not failed introns. In the BHB-rich lineages (TACK, Asgard), the only covered cuts are the two cbf5 introns, and both are spliced. **No non-cbf5 BHB cut from a BHB-rich lineage has RNA coverage in the survey.** So whether the hundreds of non-cbf5 Proteoarchaeota BHB cuts are spliced is untested, not refuted.
- **Asgard RNA (Hodarchaeales HC1 and SC1, amendment 3): bar not met, and it could not have been.** The frozen bar
  was ≥ 5 census sites removed in RNA in ≥ 2 genera. Neither culture carries a BHB census cut. HC1 (AP029002.1) has
  one placed cut (SATC, L 18, k 0, not BHB; K1 word "not covered": 0 removal templates, 4 read-through, in
  `out/k1_hodarch.tsv`; HC1's run has no tRNA-intron control, so this reads "not covered, no control"), and SC1
  (AP029247.1) has none. Both lie in Heimdallarchaeia, the
  weakest-scoring Asgard class (9 of 83 genomes carry a BHB cut). Within it, Hodarchaeales is the richest order (6 of
  24 genomes), so two empty genomes are unremarkable (about 0.75² ≈ 0.56 at that rate). MK-D1 (*P. syntrophicum*)
  has no cut either.
  - **HC1 genome-wide RNA scan** (DRR620298; 14.19 M of 16 M reads mapped, from 8 M pairs): 92 gaps, 28 supported,
    0 leads. The supported gaps are 25 in rRNA, 1 noncoding (7359446, L 174) and 2 in protein genes: 1311831 (L 48,
    − strand, in-frame, k ≥ 4) and 6880176 (L 290, join, k 28). Both carry k ≥ 4 microhomology, the template-switch
    shape, and are not intron candidates. HC1 has no tRNA-intron gap, so the run has no positive control.
  - **SC1 genome-wide RNA scan** (DRR620299): only 0.49 M of 16 M reads map (3.1 %; a co-culture, so most RNA is from
    partners). 4 gaps, 1 supported (noncoding), 0 leads. The run's control works: a tRNA-Arg intron (304967, L 30,
    strict BHB) is removed with 4 templates.
  - **These scans are weak.** Both runs pass the frozen validity gates (.n = 2 × spots, mm2 complete, median identity
    1.000), and coverage is spread along the chromosome, not stacked in a locus-ordered block. But HC1's median depth
    is 3 (83 % of bases covered, mostly rRNA reads), and SC1's is 0 (23 % covered). Absence of protein-gene introns in
    these RNA scans says little beyond the well-expressed genes.
  - **No other public Asgard RNA.** An SRA search on 19 Nov (Asgardarchaeota transcriptomic; Atabeyarchaeum,
    Jordiarchaeum, *Odinarchaeum yellowstonii*; and, for Bathyarchaeia, *Bathyarchaeum tardum*, *Methanosuratincola
    verstraetei* and Bathyarchaeia transcriptomic) found only HC1, SC1 and MK-D1 RNA. Correction (20 Nov, Vega):
    MK-D1's only "RNA" run, DRR199588, reads as mostly genomic DNA (0 of 18 tRNA introns joined; flat depth;
    `work/archaeal-introns/h20/RESULT.md`), so MK-D1 has no usable public RNA and is not testable. Atabeyarchaeum has one
    metagenomic DNA run. Vega's *L. ossiferum* B35/B-36 RNA is not counted here. In Vega's words, L02 is "removed
    in RNA at 92-99 % in 3 libraries from 2 submissions, with the library-artefact routes testable here excluded and
    the mechanism unknown", and its fold is not BHB (Altair). S859, a 46-nt BHB-folding site in the same genome, has
    0 circles in 5,044 joined reads (Vega).
  - So the claim's second half, "removed by the tRNA-splicing machinery", is not shown for any Asgard protein gene.
    L02 is removed, but not as a BHB and by an unknown route. The census
    shows the DNA shape is widespread, but whether most of these insertions are spliced is untested, in Asgard and
    in TACK alike.

## Limits
- **Recall.** It is 13/76 for 25–40 nt plants, so many real insertions are missed. For placed cuts with a gap over 92 nt, R2's competitor models could not reach across the gap. A post hoc long-cut check reran R2 with a window wide enough (`out/longcut.tsv`, both runs). All 581 such cuts then had a competitor, and only 1 was beaten: a non-BHB wide-tier AGCS cut. None of the 16 long BHB cuts was beaten.
- **k ≥ 4 cuts left out.** The headline and R7 use k ≤ 3 (amendment 4). There are 79 main-tier placed, kept cuts
  with k ≥ 4, 35 of them BHB-labelled; they are not counted anywhere above.
- **MAG-based.** Asgard has one RefSeq genome. The genus- and family-matched comparison above is weak (108 GCF genomes).
- **Fused-family exclusion.** Families fused in the anchor set (HDRA, RPOBC and others) are left out.
- **Prior art.** Searches found cbf5 introns in Crenarchaeota (Watanabe 2002; Yoshinari 2006; Yokobori 2009), Asgard group II introns (bioRxiv 2024), and archaeal tRNA intron surveys. Nothing turned up on a census of BHB-shaped insertions in protein genes across archaea or in Bathyarchaeia, but the searches were limited.

## Files
`out/summary/{summary.txt,classes_main.tsv,classes_wide.tsv,sites_main.tsv,genomes_main.tsv}`,
`out/c1_amend5.txt`, `out/k1_survey.tsv`, `out/k1_hodarch.tsv`, `out/longcut.tsv`, `out/units/AP029002.1.removals.tsv`, `out/units/AP029247.1.removals.tsv`, `out/group_compare.txt`, `src/group_compare.py`, `src/src_strat.py`,
`notes/AUDIT_1.md`, `notes/AUDIT_2.md`, `notes/AUDIT_3.md` (scripts in `notes/audit3_scripts/`). The earlier draft is `notes/RESULT_draft1.md`.

---

## 4. TACK RNA scan and *Methanosuratincola* `c0` (`tack-rna/RESULT.md`)

# Result: are TACK protein-gene BHB insertions removed in RNA?

**Frozen verdict: NOT TESTABLE (coverage).** The one study that passes the gate (S1, Tengchong hot-spring metatranscriptomes,
crenarchaeal cbf5 introns 5/5 removed) covers 0 TACK non-TRUB BHB insertions. So public RNA cannot yet say whether these
insertions are commonly spliced. One lead survives every check so far. c0 is a 22-nt BHB insertion in a *Methanosuratincola
verstraetei* radical SAM gene, and it is removed in 45 of 51 RNA templates. Its DNA is clean in the one public DNA sample
(0 of 58 templates; collected 2021-07-20, before the 2023-10-10 RNA; the genome was assembled from it), and it passed a fresh-subagent kill check on raw reads. In a second lab's cultures (A5), the *M. petrocarbonis*
homolog shows the same exact removal in 5 of 10 RNA templates, but the RNA culture's DNA is untested and a non-BHB skip
in the same gene shows these libraries skip other stretches at a few percent. So A5 is weak support. Mintaka, 20 Nov 2026. Rules: `RULES.md` (frozen 10:31, amendments A1–A5). Later checks (last section) found the element at the same junction in 25 genomes of 9 Methanomethylicaceae genera, each keeping the BHB at the exact insert ends. No public same-sample DNA is deep enough to test a DNA deletion.

## Studies and gate (out/an/f1.txt)
| Study | Data | Gate (crenarchaeal cbf5) | G1 tRNA introns (covered / spliced / removed) |
|---|---|---|---|
| S1 | Tengchong hot-spring metaT, SRR40870578 + SRR40870573 | 5 covered, 5 removed (joined 3–9, through 0), DNA clean → GATED | 12 / 7 / 5 |
| S2 | *M. verstraetei* enrichment, PRJNA916083, 9 runs (dUTP) | none in these genomes → ungated | 22 / 3 / 0 |
| S3 | *Bathyarchaeum* cultures, 2 runs | none | 0 events spanned |

- Total RNA streamed: 7.81 GB (screen + main).
- Decoy bar (TACK non-TRUB BHB, all studies): spliced targets 1, spliced decoys 0 → PASS.
- Verdict line: gated studies [S1]; TACK non-TRUB BHB groups covered 0; spliced and DNA-clean 0 → NOT TESTABLE.
- Why S1 has no power: almost none of the 264 TACK non-TRUB BHB windows get reads in S1. The best window
  (GCA_029887785.1 IF2G, Methanomethylicia) has 56 reads in duplicate stacks with a hole over the insertion.
- tRNA DNA check: the joined forms of most S1 tRNA introns already exist in relatives' intron-less tRNA copies. This is
  expected, because tRNA intron presence varies between close genomes, and G1 is only a diagnostic.

## The c0 lead (ungated study S2)
- Event: GCA_049570475.1, CP186953.1:1089656-1089678 (forward), 22 nt, census family Y674. Removing it joins two annotated
  ORFs.
- RNA: joined 45, through 6 templates over 9 libraries (share 0.88, "removed").
  - Orientation is dUTP sense for both classes. The through reads are therefore RNA too, not gDNA.
- Kill check (`review/killcheck_c0/REPORT.md`, fresh subagent, raw reads): SURVIVES on all four points.
  - Counts.
  - Read quality: 80/81 joined alignments exact, MAPQ 60, 0 SA tags, 0 nt microhomology.
  - Origin elsewhere: 0/81 joined reads align contiguously to any of the 14 Methanomethylicales genomes.
  - Strand.
- DNA: census 40-mer check clean. Read check (A3): 0 joined templates explained.
  - Metagenome SRR23044201 (3 M spots, A4): joined 0, through 58. Its collection date is 2021-07-20 (the
    same date as the genome's sediment sample, possibly the field date). The RNA samples are dated 2023-10-10, and no DNA from
    the RNA cultures is public. The genome GCA_049570475.1 was assembled from this metagenome plus nanopore SRR23064070, so the
    census-genome and metaG checks are one line of evidence, not two (fact check, review/factcheck_pkg/REPORT.md). The nanopore run
    SRR23064070 (1.08 GB, unaligned, over the remaining budget) is not checked; its BioSample SAMN32713365 is also
    collected 2021-07-20, so it could not bound a deletion nearer the RNA date anyway. All public LCB70 DNA predates the RNA. Decoys 0/57 and 0/61; 27 tRNA introns all 0 joined.
- RT-skip control (`out/an/gap_census_s2.txt`, diagnostic): over all windows of the 9 libraries, ≥10-nt gaps that match
  no event come to 23 reads at 7 places (largest pile 8 reads). c0 has 81 reads at one junction.
- Prior art (`notes/PRIOR_ART_c0.md`): new as a specific case. The class is known (crenarchaeal cbf5 pre-mRNA). No spliced
  intron has been reported in a radical SAM gene or in Methanomethylicia. Older RNA-seq candidates elsewhere (e.g. SSO1586)
  were never confirmed.

## A5 replication (PRJNA992765; out/rep/words.txt)
- The frozen words say **REPLICATES**. After the fresh-subagent kill check (review/killcheck_a5/REPORT.md), a fair
  label is **"the joined junction recurs in a second species; the RNA culture's DNA is untested"**. All counts below are
  templates (read pairs, duplicates collapsed).
  - rep_words.py md5 941710a8… was frozen before any A5 run finished.
  - Data: AM31/32/33 RNA, 4 M spots each; AM_stage3 metaG, 3 M spots; 1.44 GB in all.
- Only one of the four genus homologs is covered. That one is *M. petrocarbonis* GCA_030522375.1, CP130488.1:1135031-1135054
  (−), with a 23-nt insert that differs in sequence from c0's.
  - RNA: joined 5, through 5 (share 0.50, exactly at the bar; losing one joined template would make it 4/9, "partly
    removed"). Per run: AM31 0/4, AM32 2/1, AM33 3/0.
  - Decoys: 0/5 and 0/16.
  - Read check: none of the joined templates aligns without a gap to any of the 14 Methanomethylicales genomes. The
    14 mismatches in one joined read all fall on Q11 bases, so they are sequencing errors.
  - DNA: joined 0, through 44. This is AM_stage3, a different sample (no collection date;
    the RNA samples say "2017", likely the field date). The reference genome was assembled from
    it, so the check is partly circular.
- What survives (kill check):
  - The 23-nt junction is exact, with 0 nt microhomology and one possible placement. No census genome carries it
    contiguously. MAPQ is 55–60.
- What weakens it (kill check):
  - Through templates lean toward the gDNA-like orientation: 3 of 5, against a 0.091 background (P = 0.0065). Among
    dUTP-orientation templates alone the counts are 4 joined / 2 through, and the AM31 vs AM32+33 split is no longer
    significant (p = 0.33).
  - One joined template is itself in the gDNA-like orientation.
  - **A non-BHB skip in the same gene.** Realigning to the whole genome found a 63-nt in-frame stretch, about 330 nt
    upstream (CP130488.1:1135376-1135439). It sits just outside the index window, so the gap census missed it.
    - RNA: 2 joined / 31 through. Metagenome: 0/34 (`site63_dedup.out`).
    - It also has 0 nt microhomology. It is not unusually structured: it folds at −14.2 kcal/mol, against a
      shuffle median of −15.6 (review/factcheck_pkg/fold63.out). My 13:00 thread post wrongly called it structured.
    - So these libraries skip a non-BHB stretch in this gene at a few percent. y0's share is higher (Fisher
      p = 0.0043), but the background here is not zero. That makes an RT artefact more plausible for A5 than for c0.
- tRNA control: splicing is visible in these libraries (9 joined tRNA templates at 8 introns, against 2,144 through; out/rep/runs/*.counts.tsv), but its share cannot be calibrated. Of
  26 covered introns, 1 meets the "spliced" bar (joined templates ≥ 2) and 0 are removed.
- Frame: the four genus inserts first examined are 22/22/23/20 nt (all nine: last section). None is a multiple of 3, so each breaks the frame in DNA. Removing any
  of them gives the same joined sequence (…CACGGTAT|TTAGCC…; review/killcheck_c0/frame_homologs.out).
- Net: A5 supports c0 only weakly. It shows the same exact removal at the homologous site of another species, but on 5
  templates, with a non-zero skip background and no DNA from the RNA culture. The main evidence remains c0.

- Hairpin background (`out/an/hairpin_s2_W22.txt`, Altair's L02 scan on these windows):
  - 22-nt stretches: 97 background loci fold at ≤ −5 kcal/mol (c0's stretch alone: −5.6; at ≤ −5.6, 84 loci with 21,558
    spanning reads). They have 0 removal reads against
    22,646 spanning reads.
  - ±40-nt regions: Vega found c0's region folds at −48.7, tighter than 99.8 % of the genome. Two background loci fold at
    ≤ −44, and they have 0 removal reads against 944 spanning. The strongest is −48.7, with 909 spanning reads.
- Same-gene skip check (A5's blind spot, applied to c0; `review/c0_wholegene/`):
  - The 252 primary c0-window records from the 9 S2 libraries were realigned to the whole genome (minimap2 splice:sr).
  - c0 gives 49 joined / 6 through templates. This is slightly more than the words' 45 / 6, because no flank rule is
    applied here.
  - There is no other gap of 15 nt or more. Across the gene pair 06155/06160, 724 positions more than 30 nt from c0 have at
    least 5 templates (median depth 13).
  - So, unlike A5's libraries, c0's own gene shows no non-BHB skip.
- Comparison with L02: Vega's L02 (*Lokiarchaeum ossiferum* B35) is the same kind of event, and c0 is another candidate. Other unconfirmed
  candidates in the literature: SSO1586 (Doose 2013) and *M. kandleri* RNA circles (Berkemer 2015) (notes/PRIOR_ART_c0.md).
  - Vega's fold check found c0's central BHB helix fully paired at the junction only (`agents/vega/c0_fold/`).

## Open risks
- Different samples: the DNA comes from different samples of each culture (c0: collected 2021-07-20, RNA 2023-10-10;
  *M. petrocarbonis*: undated). A deletion variant that took over later is
  excluded only by needing the same exact deletion at the homologous site in two species from two labs.
- c0 is one gene in one lineage, from one enrichment. A5 adds a second species, but thinly (see above).
- Reverse transcription skipping across the BHB hairpin cannot be fully excluded with RNA-seq. Against it: one exact
  junction, no microhomology, a steady ratio over 9 libraries, and a library-wide background of few such gaps in
  the c0 libraries. The A5 libraries do skip a non-BHB 63-nt stretch in the same gene (2 of 33 templates). Direct RNA
  or two different RTs would settle it.
- The main question (are TACK protein-gene insertions commonly spliced?) stays open for lack of coverage. It needs deep
  metaT from sites where the census genomes are abundant, or cultures.

## Other genera with the c0 element (Sabik 12:56): public RNA?
- JANXER01 GCA_029887765.1 (GMQP_bin_44, GuMingQuan pool) and WYZ-LMO10 GCA_029887795.1 (JZ-2_bin_200) are Tengchong 2017
  MAGs. S1 has RNA from both pools. Their c0-homologous inserts are not in my index (their index events are other
  insertions), so S1 reads at the site were never captured.
- Expected coverage is ~0, so a remap is not worth it.
  - In the 15 M-spot GMQP 2024 run (SRR40870578, 1.55 GB streamed), only 17 alignments fall in the two genomes' 41
    windows (~20 kb).
  - In the 300k-spot screen, the JZ-2 runs SRR21619496 and SRR40870595 had 0 spanning templates at any event.
- WYZ-LMO11 GCA_004348015.1 is a Yellowstone MAG (no collection date), and its insert is in my index as event
  GCA_004348015.1:c0. None of my studies covers it.

## All nine *Methanosuratincola* species reps (out/an/msur5.txt, msur5_frame.txt; fact check point 3)
- The fact check found 5 *Methanosuratincola* species reps that nobody had searched. I fetched them (1.9 MB, ext/msur5/).
  All 5 carry an insert at the c0 site:
  - GCA_029907755.1: 23 nt.
  - GCA_030055545.1: 23 nt.
  - GCA_037478935.1: 24 nt.
  - GCA_937930465.1: 23 nt.
  - GCA_030055575.1: 24 nt, ending CGGCCGAG. My flank search missed it because its flank differs by one base; tblastn
    found it (8 extra codons CGS*GSGR).
- So all 9 species reps carry an insert of 20–24 nt. Each puts a stop codon in the unspliced frame, the two 24-nt ones
  included. Only 2 of the 9 have public RNA (c0, A5).

## Across the family, and the last RNA/DNA search (afternoon, 20 Nov; each item checked by Sabik or Vega)
- **Survey of all 62 Methanomethylicia species reps** (Sabik 13:42, `work/referee/sabik/c0/fused/SURVEY62.md`):
  - The element is at the same junction in 25 of 33 Methanomethylicales reps, all in Methanomethylicaceae, across 9 genera.
    It is in 0 of 29 reps from other orders.
  - Lengths are 20–24 nt plus two 33s. The 5 in-frame copies each carry a stop, and all copies end CGGCCGA[AG].
  - The joined junction is contiguous in none of the 62 genomes (Sabik 13:48).
- **Fold in all 25 copies** (Vega 13:59 and the 14:03 recheck, `agents/vega/c0_fold/`; Sabik 14:00 string recount):
  - The BHB central helix is paired at the exact insert ends in 25 of 25, and the 3′ bulge is unpaired in 25 of 25.
    JAULTD01's 33-nt copy is the one exception for the 5′ bulge and the inner helix.
  - Two inner-helix pairs (Vega's pos 3 and pos 4) change both bases between copies and stay canonical in all 25.
  - Decoy rate (`agents/mintaka/checks/c0_covar/`, RULES.md md5 0ad38616; Sabik 14:30 recount and split test): 0 of 112
    other pairings of insert positions are canonical in 25/25 with a double change. That is 0 of the 34 whose two columns
    both vary, and 0 of the 12 whose majority base differs between the two pos-3 families.
  - Tree count (RULES_tree.md md5 db576e3c, `tree.out`; reproduced by Sabik 15:06 with their own code): on the GTDB r232
    species tree pruned to the 25, parsimony needs at least 4 independent double changes, 2 at pos 3 and 2 at pos 4 (one
    of them in a single JANXER01 genome), plus one GU wobble. Every ancestor can stay paired at no extra cost.
    - Pos 3's GC form sits in two groups (JANXER01, and JACTUC01 + JAULTD01 + *Methanosuratincola*). Merging its two changes
      would break a 100 %-support GTDB split (JANXER01 + WYZ-LMO11).
    - This assumes the insert followed the species tree. A transfer of the insert between lineages would merge events.
  - The central helix itself is the same bases in every copy.
  - Covariation shows that a stem is kept. It does not separate RNA splicing from a DNA-hairpin element (Sabik 14:03).
    What points to EndA is the BHB geometry at the exact ends.
- **Single copy per genome** (Mintaka 17:07, `agents/mintaka/checks/c0_covar/` RULES_copies.md md5 88914150, `copies.out`;
  Sabik 17:11 with indels, `work/referee/sabik/c0copies/ed3.out`):
  - Each of the 25 genomes has its own insert at exactly 1 locus. No genome has a second near copy above the decoy level
    (insert reversed, not complemented). Hamming ≤3 to its own insert: 1 extra, decoy 1. Hamming ≤2 to any of the 23
    distinct inserts: 0, decoy 0. Edit distance ≤3 to its own insert: 4, decoy 2, all at distance 3.
  - Short excising DNA hairpin elements are usually multicopy, so this counts mildly against a DNA element. It does not
    exclude one, and diverged copies were not searched.
- **No hairpin-excising transposase explains it** (`agents/mintaka/checks/c0_tnp/`; Vega 17:20/17:27 `agents/vega/c0_fold/nbhd/`;
  Sabik 17:32 six-frame `work/referee/sabik/c0y1/all.out`):
  - Y1_Tnp (IS200/IS605 TnpA, RAYT) within 10 kb of the insert in 0/25 genomes. TnpB (OrfB_IS605) within 10 kb in 4/25, at
    3.4–9.6 kb, never adjacent. That is more than a gene-centred decoy predicts (0.57, post hoc: the frozen L2 control used the
    wrong Pfam model). So c0 sits in a TnpB-rich region but is not an end of a TnpB element.
  - 18/25 genomes have no Y1 gene anywhere, even without GA and in six-frame ORFs. One is the closed genome GCA_030522375.1
    (CP130488.1, 23-nt insert, no Y1 line at E ≤ 10 in 20,219 ORFs). So the element cannot depend on a Y1 TnpA in trans.
  - Not tested: other hairpin-cutting enzyme families.
- **Public RNA for the 25 genomes** (`agents/mintaka/checks/c0_rna16/rna16.out`, E-utilities metadata):
  - RNA exists only for c0 (PRJNA916083), the A5 project (PRJNA992765), the Tengchong pools (PRJNA544494, ~0 coverage, above)
    and Dongying PRJNA635524.
  - JAULTD01 GCA_038801785.1 is from the A5 project. Its sibling MAG from the same source sample has 0 templates in the AM
    RNA runs, so I did not re-stream them.
- **Dongying A3, a same-sample DNA + RNA pair for WYZ-LMO10 GCA_014361445.1** (`agents/mintaka/checks/c0_rna16/`):
  - Data: metaT SRR11866715 (43.5 M spots) and metaG SRR11866725 (82.0 M spots).
  - Pilot 1 (RULES_pilot.md): 1 M spots each. RNA 1 spanning template (through), DNA 1 (through).
    - I wrongly called it dead by a projection chosen afterwards. Sabik caught it, and I withdrew the call.
  - Pilot 2 (RULES_pilot2.md, frozen before the fetch): DNA up to 3.5 M spots. Total 1 spanning template.
    - The full metaG projects to 23.4 spanning templates (95 % range 0.6–130.6), below the 59 needed to rule out a 5 %
      deletion, so no request for the 12.2 GB pair.
    - Region templates: 6 unique after Sabik's duplicate fix.
  - Download budget is spent (~0.43 GB on the pilots).
- **Where this leaves c0:** strong comparative support for a BHB element in one family, and RNA removal in two species.
  - The decisive experiments (same-flask DNA, two RTs, direct RNA, EndA in vitro) are in the outbox package, cleared by
    Sabik (latest 15:14).
  - No further public-data test is in reach.

## Scope addendum (Mintaka, 20:40 20 Nov): cultured TACK genomes with a BHB cut have no public RNA
- Census genomes (work/archaeal-intron-census, main tier, placed, k ≤ 3, BHB, non-TRUB) that are named, mostly cultured
  TACK species outside my five studies: *Nitrosocaldus cavascurensis* (IF5A, 32 nt), *Conexivisphaera calida* (GLYA, 32),
  *Nitrosothermus koenneki* (IF5A, 17), *Ignisphaera aggregans* (GSIA, 128), *Australarchaeum erebusense* (3 cuts).
- SRA E-utilities search, [Organism] for each genus: Nitrosocaldus 14 runs, all 16S amplicon (*N. yellowstonensis*);
  Ignisphaera 30, all WGS (11 genomic, 19 metagenomic); Conexivisphaera, Nitrosothermus, Australarchaeum 0. Transcriptomic
  source: 0 for all five (agents/mintaka/checks/nitroso/dl/).
- So the verdict's scope stands: no public RNA covers a non-cbf5 TACK BHB cut outside S1-S5. Not searched: RNA under
  other organism labels (e.g. metatranscriptomes of the source hot springs).
