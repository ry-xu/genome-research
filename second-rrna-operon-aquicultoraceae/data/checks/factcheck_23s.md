<!-- Copied from the lab record. Paths that point to files in this repository were rewritten to repo-relative ones; any other path refers to the original working folder and that file is not included (see data/README.md). -->

# Fact-check: divergent 23S draft (the 23S post)

Fresh subagent (task a0bdc98d9bf1d8317), finished 14 Nov ~19:29. Its scripts are copied to a scripts folder (not included). Saved verbatim by Deneb; all fixes applied to post.md and src/frozen_core/fig_div23.py.

---

I found 18 claims in the draft that are wrong or overstated, and one figure defect: the footnote is cut off at the right edge, ending "…viral-fraction metatranscripto". The fix is to wrap the footnote text or move it left of x=0.47 in `src/frozen_core/fig_div23.py`. All ten citations check out on authors, year, journal, volume, pages and quotes, though two need rewording. The For-the-curator section has three misstatements of who checked what.

Every number below comes from scripts I ran or text I quote. Scratch scripts and outputs are in the checker's session scratchpad (not kept). I sent no email or personal details in any request and modified no files.

## 1. Wrong or overstated

| Line | Claim | What the source says | Suggested fix |
|---|---|---|---|
| L7 | 15,739 genomes, "which are the lab's cached genomes" | STATE line 255: the census is 15,739 genomes, about 8% of GTDB bacterial reps, a subset of the ~50,600 cached. This also conflicts with L50 "all ~50,600 cached". | "a 15,739-genome census (one 23S per genome)" |
| L7 | "a dozen PTC positions" | STATE line 432 lists 14 changed positions (2062, 2251–2253, 2447, 2449, 2451–2453, 2503, 2504, 2583–2585). | "14 PTC positions" |
| L7 | "soil MAG" | GTDB `isolation_source` is "none". STATE line 444 says the reads are PRJNA620721, a UBC soil methanol microcosm; I did not verify that. | "MAG from a soil methanol microcosm (per its read project)", or drop "soil" |
| L1, L41, L97 | "about 72%" identical; "28% divergent, roughly three times the largest case" | 72% is BLAST local identity (dc-megablast HSPs 72.4%/919 nt, 69.4%/1183). Identity over aligned columns is 0.792 (`pairing_free.out`), which is 21% divergent. Yap's 10% is a whole-gene figure, so like for like it is about 2×, not 3×. | Use 79% (21% divergent, about twice Yap's 10%), or say the 72% is local BLAST identity and drop "three times" |
| L9 heading, L13 "absent outside it", L42 | Copy is absent outside the family | Whole-genome search covered only the 65 Aquicultoria genomes (`loci.py` uses aq65). Everywhere else only one census 23S per genome was searched, and the second-copy scan is about a third done (L50 itself says so). | "not found in the one-23S-per-genome census or in the 65 Aquicultoria genomes" |
| L13 | No other census 23S "comes closer than 69–74% over at most 1.3 kb" | `q_vs_lsuall.tsv`: kmv_div23 reaches 80.45% over 578 nt (GCA_031458695.1) and cad_div23 reaches 79.4% over 597–658 nt. Only dx_div23 stays within 68.9–74.5%. The 1.31 kb maximum length is right. | "none closer than ~80% over ≤0.66 kb or ~74% over ≤1.3 kb" |
| L14 | The loop is "gone in every divergent copy" | I checked 6 sequences from 5 genera (dx, kmv, cad, cad2, dqq, aqA): none contains GATAACAGG. JASEIR01 was not checked. | "in the 6 copies checked (5 genera)" |
| L16 | "13 functional 23S genes … score 87–93%" | 87–93% (0.869–0.930) comes from 16 controls (13 s_ctl plus jausch, yes70, uba9087). The 13 s_ctl alone score 0.894–0.923. The kill report's "controls at most 0.91" is also wrong: GCA_050945135.1 is 93/97 = 0.96. "About twice the controls' rate of pair-breaking changes" is loose: divergent ratios are 1.35–1.85 against a control median of ~0.63. | "16 controls score 87–93%"; give the ratios instead of "about twice" |
| L10 | "every read agrees" | True only for reads at ≥97% identity. `pile.out` has 4 lower-identity reads reading A at 2452. | "every read at ≥97% identity agrees" |
| L19 | The RNA carries cellular rRNA "only as a small by-catch" | Read halves aligning to the 23S references are 9.8% (SRR11558425), 9.4% (542) and 6.8% (543) of all read halves, against 0.24% in DNA, about 30–40× more. Only this species' share is tiny (86–1,237× rarer). | "23S is 7–10% of read halves, but this species' share is 86–1,237× lower than in DNA" |
| L21 | "151 of 163 … match exactly" | Reproduces only as "at least one passing alignment with 0 mismatches" (151/163 RNA; 533/577 DNA). If all passing alignments must be exact it is 128/163 vs 422/577; with the NM tag DNA is 531. Too weak to rule out a relative on its own. | State the definition; lean on the audit's spacer reads as the evidence against a relative |
| L24 | "7,466 read pairs fall on protein genes" | 7,466 is a per-gene sum, so pairs overlapping two genes count twice. Distinct read pairs = 7,063 (`genecount.py`, same filters and exclusions as the audit's `magcount.py`). | 7,063 read pairs |
| L26 | "5 read pairs cover the 3' end of the divergent 16S" | 4 pairs cross the 16S 3' end (23,831) into the spacer. One (SRR11558425.14562365, spanning 23,876–24,012) stops ~45 nt short of the end. | "5 pairs at the 16S 3' end, 4 running across it into the spacer" |
| L31 vs L35 | "every case I found is described as functional or expressed" | The list then includes Pei's *B. afzelii* copy that "appears to be a pseudogene". | "every case except one *Borrelia* 16S copy…", or move the Pei bullet out of that list |
| L37 | G2252 and G2251 pair with C74 and C75 (Samaha 1995; Feinberg &amp; Joseph 2006) | Samaha 1995 shows only "A Watson-Crick G-C pair between G2252 … and C74". Feinberg &amp; Joseph state both pairs as background, but their experiments test G2252–C74. They find it "not essential for peptide bond formation" but important for peptide release. | Cite Samaha for G2252–C74 only; attribute G2251–C75 to structural work. Feinberg &amp; Joseph's finding also supports the post's "not a dead ribosome" caveat. |
| L50 | UBA6777: "nobody has checked it" | STATE line 441 has Deneb's own analysis: 5 MAGs from 4 projects, identical divergent 23S, pairing 0.927, no normal copy binned. Deneb posted it to the thread at 17:3x (STATE line 440). | "checked by me only (MAG count, pairing score); no reads or independent check" |
| L82 | "Restricted to the region all 16 loci cover, the bound is 0.040" | AUDIT's row is 558–2122, but the PEXG locus starts at 575. | "the region nearly all loci cover", or recompute from 575 (minor) |
| L94 | "PRJNA622041–4" | The three runs are PRJNA622041 (SRR11558542), 622042 (543) and 622044 (425). 622043 is not one of them. | "PRJNA622041, 622042, 622044" |

## 2. Verified, and how

**Expression counts and statistics**
- 545:577 (DNA), 0:129, 0:11 and 0:23 (RNA): recounted by distinct read name, identity ≥0.98, alignment length ≥60 (`recount.py`). These match the figure's right panel.
- Expected 79.18, P = 8.4e-48, per-copy bound 0.0196, and 86–1,237× (85.6, 480, 1237): all recomputed.
- Sample dates: RNA 25, 26 and 27 May 2015, DNA 28 May; same lat_lon; "1–3 days before" holds (run XMLs).
- RNA runs are titled "Metatranscriptome of freshwater viral communities…", RNA-Seq with RT-PCR selection. The DNA run is SAMN07631077, PRJNA405986.

**Strandedness and spacer reads**
- The libraries are dUTP-stranded: on genes, R1 is opposite in 6,445 and same-strand in 293; on the normal 23S, R1 is on '−' 20/20.
- The spacer reads are therefore sense: all R1 are on '+' while the operon is on '−'. Their identity is 97.4–100%, and the mapq-60 pairs span 23,711–24,012.
- No alignment touches the divergent 23S. Normal-operon reads stop by about 11,830, none run into the spacer, and genes on both sides of the divergent operon carry reads.

**Figure and census**
- Every left-panel conservation value matches `columns.tsv`, and the 99.7–99.98 range is right.
- Crystal Geyser PTC bases (C, UUU, U, A, A, A) agree across all 3 routes the post names. The route through the genome's own normal copy (aqA_norm, which has a ~985-nt insertion) disagrees, but the post does not claim it.
- `ref.fa` has 51 loci + 17,329 decoys. The 5 shared PTC changes are in `ptc_pairwise.out`.

**Kill-report and audit results**
- Prevalence is 23/18/15 by genus across six genera.
- Synteny: flanks hit divergent loci in 5 other genera and 0/9 at normal loci.
- Anti-Shine-Dalgarno: all 6 divergent 16S (5 genera) lack it; all 4 normal have it.
- Structure: 0.832/0.843/0.864 against a null maximum of 0.743; "2.2× faster" holds.
- GTDB: 22 of 23 genomes are "derived from metagenome"; GCA_902810715.1 is "none".
- The For-the-curator items attributed to the kill agent and the audit match `REPORT.md` and `AUDIT.md`. That includes: the remap was run on the 27 May sample only; overdispersion bounds 0.038/0.093/0.28; DASXVE01 reads were checked by Deneb (STATE line 432) and the kill agent (`REPORT.md` line 22).
- The thread times 18:51 and 18:57 match the lab discussion thread (not included) lines 4735 and 4781. STATE's "18:52" is the thing that is off, not the post.

**Citations** (via Europe PMC, PMC, Crossref and LPSN)
- **Yap, Zhang &amp; Wang 1999**, J Bacteriol 181:5201–5209: "The coding sequences for the 16S and 23S rRNA genes differ by approximately 6 and 10%". Functionality is "concluded", and they "propose that T. chromogena acquired rrnB operon … via horizontal gene transfer". Correct.
- **Johansen et al. 2017**, PLoS One 12:e0186393: "(7.3-9.0%)"; 23S "highly divergent … (8.7%), and likewise had nearly identical secondary structure". Correct. Optional: *S. hyalinum* is now *Kalymmatonema* (Bohunická et al. 2025, J Phycol 61:1225).
- **López-López et al. 2007**, J Mol Evol 65:687–696: operon B shows "an expression level four times higher than the other two together at 50 degrees C and three times lower at 15 degrees C". Correct.
- **Pei et al. 2010**, Appl Environ Microbiol 76:3886–3897: "the diversified gene in B. afzelii (20.38%) appears to be a pseudogene". Correct.
- **Vignolini, Carril … Schwartz, Baumgarten**, bioRxiv: Crossref lists the DOI prefix 10.64898 under openRxiv, posted 2026-06-18. doi.org redirects to biorxiv.org. Abstract: "appears to elongate inefficiently, likely due to divergences in the peptidyl transferase center obstructing the exit tunnel"; "a dominant-negative repressor of translation". The quote is exact, truncated before "obstructing".
- **Samaha, Green &amp; Noller 1995**, Nature 377:309–314. Correct for G2252–C74 (see L37 above).
- **Feinberg &amp; Joseph 2006**, J Mol Biol 364:1010–1020 (see L37 above).
- **Green, Samaha &amp; Noller 1997**, J Mol Biol 266:40–50: "Mutations at either position confer a dominant lethal phenotype when the mutant 23 S rRNA is coexpressed with the endogenous wild-type". Correct.
- **Polacek et al. 2001**, Nature 411:498–501: "large ribosomal subunits with mutated A2451 showed significant peptidyl transferase activity". Correct.
- **Jiao et al. 2021**, ISME J 15:3005–3018: LPSN gives "*Candidatus* Aquicultor" Jiao et al. 2021, not validly published, type species "Ca. Aquicultor primus". So L38's attribution is right.

## 3. Could not verify
- That DASXVE01's reads come from a soil methanol microcosm (PRJNA620721). STATE says so; I did not check the BioProject.
- The kill report's NCBI nt BLASTs, which never returned.
- The remaps of SRR11558542 and SRR11558543, and the SRR647* runs. I downloaded no reads.
- The loop motif in JASEIR01.
- UBA6777 beyond what STATE line 441 records.
- The library protocol behind "RT-PCR" selection (L97), and which paper the metatranscriptomes belong to.
- PRIORART Q4 ("no such copy reported"), which rests on web snippets only.

## 4. For the curator: misstatements of who checked what
- **L93 "Not checked by anyone: UBA6777"** (and L50 "nobody has checked it"): Deneb did a first-pass analysis (STATE lines 440–441) and posted it to the thread. It should read "checked by me only".
- **"Checked by me only: the RNA identity pattern (151/163 exact)"**: it reproduces only under the loosest definition (any exact alignment). Under all-alignments-exact it is 128/163. Name the definition.
- **"Checked by me only"** for the PTC bases, the 151/163 pattern and the spacer-read coordinates: I have now independently reproduced all three, so they could say "checked by me and the fact-check agent", subject to the 151/163 caveat.
- **L94**: the BioProject range "PRJNA622041–4" is wrong, as in the table.
- The attributions to the kill agent and the audit, and the 18:51/18:57 times, are accurate.

Sources:
- [LPSN genus "Candidatus Aquicultor"](https://lpsn.dsmz.de/genus/aquicultor)
- [LPSN species "Candidatus Aquicultor primus"](https://lpsn.dsmz.de/species/aquicultor-primus)
- [Europe PMC REST search](https://www.ebi.ac.uk/europepmc/webservices/rest/search)
- [Jiao et al. 2021, PMC8443620](https://pmc.ncbi.nlm.nih.gov/articles/PMC8443620/)
- [Candidatus List No. 6, PMC12281948](https://pmc.ncbi.nlm.nih.gov/articles/PMC12281948/)
- [Crossref record for 10.64898/2026.06.17.732804](https://api.crossref.org/works/10.64898/2026.06.17.732804)
