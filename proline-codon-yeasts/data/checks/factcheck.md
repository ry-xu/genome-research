# Fact-check of the proline post (fresh subagent, 2026-11-11 ~13:00)

Saved by Clover from the subagent's final message (the harness blocked it from writing this file). Its scripts are in this folder (copied from the scratchpad). All items were applied to post.md and README.md at 13:2x.

# Fact-check: outbox/2026-11-11-mito-proline-codons/post.md and image.png

## Summary: items to fix

**WRONG**
1. **Figure caption, "81 positions in 13 mitochondrial proteins".** It is 12 proteins. The figure code counts profile models, and cox1 has two models (Vir with 11 positions, Oth with 1).
2. **Census, "all 18,015 in RefSeq".** It should read "18,015 of the 18,031 RefSeq records" (README).

**MISLEADING / LOOSE** (replacement wording is in the full list below)
3. **"Annotated with the mold mitochondrial code (table 4…)".** The four GenBank BK records have no gene annotation at all: only a source feature, no CDS, no transl_table.
4. **Figure caption omits the column filter.** The 81 positions are a subset of 126. The orange counts in the figure (67/78/77/42 for *W. infanticola*/*W. siamensis*/*W. verensis*/*S. pinicola*) differ from the text (78/82/82/44), and nothing explains why.
5. **"Never where they have leucine".** At conserved leucine positions the count is 0. But each of the three *Wickerhamiella* has one CUN codon at a position where leucine is only a bare majority (52%). The audit's pairwise check also found one in *W. infanticola* (1 of 462 positions).
6. **"None prefers CUA at proline … lowest p 0.019".** Two of the 147 genomes have p &lt; 0.05 (0.019 and 0.021), each caused by a single CUA.
7. **Blinded re-check, "told only that the yeast mitochondrial code was one candidate".** The setup was table 4 as the default and table 3 as the alternative, which points at CUN, AUA and CGA/CGC. Organism names were also visible in the file headers.
8. **"NCBI lists 26 such alternative codes" and "mostly mitochondria".** There are 26 tables besides the standard code, but table 11 changes no codon's meaning (only start codons differ). Only 13 of the 26 are mitochondrial, so "mostly" is too strong.
9. **Note, "both come from sugar-cane leaves in Thailand".** PYCC 7069's sample record says only "Thailand" and gives no source. The species was also described from six Thai sugar-cane strains, so a shared origin doesn't point to the type strain.
10. **Wolters ratio of amino-acid-changing to silent substitutions.** The paper names two groups as highest: the *Wickerhamiella/Starmerella* clade and "Other Dipodascales". The post mentions only the first.
11. **Minor wording.**
    - Ling 2014 says the alanine tRNA "most likely" derived from tRNA-His; the post states it as fact.
    - The figure legend says "CTN" while the text says "CUN".

## Verdicts, sentence by sentence

### Bold lead
- **"In the mitochondria of four yeasts, the codons CUU and CUA mean proline, not leucine."** SUPPORTED at the DNA level (pro_leu_wss.tsv, cun_CTA/CTT.tsv, blind re-check, audit).
  - Optional hedge: the post's own "Not known" section says the evidence is DNA only.
- **"No known genetic code reads any codon other than CCN as proline."** SUPPORTED for NCBI: all 27 tables in gc.prt have proline = CCN only.
  - Beyond NCBI it rests on searches. Optional wording: "No code in NCBI's list, or in any report I found, …".

### Paragraph 1
- **"Genes are read three letters…"** SUPPORTED.
- **"61 codons … 20 amino acids … 3 are stops."** SUPPORTED.
- **"Some lineages, mostly mitochondria, … NCBI lists 26 such alternative codes."** MISLEADING/LOOSE.
  - There are 27 tables, so 26 non-standard ones.
  - Table 11 uses the standard amino-acid assignments, so only 25 change a codon's meaning.
  - 13 of the 26 are mitochondrial: tables 2, 3, 4, 5, 9, 13, 14, 16, 21, 22, 23, 24 and 33.
  - Fix: "Some lineages, many of them mitochondria, read a few codons differently; NCBI lists 26 alternatives to the standard code."
- **"…changed meaning twice: to threonine in baker's yeast and many relatives, and … *Ashbya gossypii*, CUU and CUA to alanine."** SUPPORTED (NCBI table 3; Ling et al. 2014, NAR 42:499).

### Paragraph 2
- **"I found a third meaning."** SUPPORTED as far as searches go; novelty depends on literature searches.
- **"…four species of the yeast family Trichomonascaceae … CUU and CUA sit where related proteins have proline…"** SUPPORTED.
  - NCBI Taxonomy places both genera in Trichomonascaceae and assigns mitochondrial code 4.
  - Every CUN codon at a strict proline position is CUA or CUU:

    | genome | CUA | CUU |
    |---|---|---|
    | *S. pinicola* (262) | 33 | 11 |
    | *W. infanticola* (532) | 62 | 16 |
    | *W. siamensis* (538) | 70 | 12 |
    | *W. verensis* (542) | 69 | 13 |

- **"…and never where they have leucine."** MISLEADING/LOOSE. I recomputed over every aligned CUN codon in the four genomes (126, 130, 135 and 70 codons).
  - None sits where 70% or more of mitogenomes have leucine.
  - Each *Wickerhamiella* has exactly one at a position with a 0.518 leucine share. The highest for *S. pinicola* is 0.494.
  - The audit's pairwise check found 1 in *W. infanticola* at a leucine position in its relatives (n = 462).
  - The census (90% threshold) and the blind re-check both give 0.
  - The control *Wickerhamiella* sp. BK074539 has 93 at 90%-leucine positions.
  - Fix: "…and never at conserved leucine positions."
- **"Leucine is spelled UUA/UUG instead."** SUPPORTED: at strict leucine positions, 0 CUN and 140–166 UUA/UUG.
- **"The four sit on two separate branches."** SUPPORTED by the tree.

### Bullets
- **Per genome, "44 to 82 CUU/CUA codons…"** SUPPORTED: 78/82/82/44 for *W. infanticola*, *W. siamensis*, *W. verensis*, *S. pinicola*.
  - **"none at positions where they have leucine, … UUA/UUG 140 to 166 times"** SUPPORTED at the bullet's 90% definition (140/149/147/166).
  - Optional addition: *S. pinicola* still uses CCN at 36 of its 80 proline positions. The figure shows this mix; the text doesn't say it.
- **In relatives, "the other 17 *Wickerhamiella* use CCN at 80–85 … CUN at 0–1"** SUPPORTED. BK074539 is the one with 1.
  - **"147 other … Dipodascomycetes"** SUPPORTED: 151 genomes, selected by taxonomy.
  - **"none prefers CUA at proline … lowest p 0.019"** MISLEADING/LOOSE.
    - *Geotrichum* BK072907 has p = 1.9e-2, from 1 CUA at proline and 2 at leucine.
    - *Tardiomyces* BK074691 has p = 2.1e-2.
    - Fix: "none comes close (lowest p 0.019, from a single CUA at a proline position, among 147 tests)".
  - **"the four: p 2.1e-28 to 1.9e-46"** SUPPORTED (*S. pinicola* 2.1e-28, *W. infanticola* 2.8e-39, *W. verensis* 4.4e-45, *W. siamensis* 1.9e-46).
- **Raw reads, "every one of 286 sites … 3,782 to 46,781 reads … 98.9–99.6%"** SUPPORTED (recomputed from reads/*.tsv).
  - Reads per genome: 3,782 (*W. infanticola*), 46,781 (*W. siamensis*), 13,970 (*W. verensis*), 16,746 (*S. pinicola*).
  - Agreement: 98.86%, 99.60%, 99.57%, 98.92%. No site has a non-matching majority.
  - Note: coverage is thin in *W. infanticola*. Two sites have a single read (cox1 Vir:149 CUA and Oth:260 CUU).
- **Second sequencing, "another lab's *W. siamensis* culture (PYCC 7069, 2024) … all 82 sites (1,432 of 1,450 reads)"** SUPPORTED.
  - The run is ERR11532760 from "Yeast Genomics Lab@Nova, UNL, Portugal", public 2024-05-07, sample alias PYCC 7069.
  - The assembly's own run is from the Y1000+ project, strain CBS 13331.
  - Note: 3 sites have fewer than 3 reads; the lowest has 1 (nad5:12, CUU).

### Figure (identical to results/cun_pro/pro_grid_tree.png)
- **49 tip names and tree:** SUPPORTED. All names match NCBI, and the figure table reproduces with 0 mismatches.
  - Optional: *W. azyma* sits outside *Blastobotrys* on a 2.37 branch, likely a long-branch artefact. No support values are shown.
- **Legend "CCN (proline in every code)":** SUPPORTED.
- **Legend "CTN (leucine in the standard code)":** content SUPPORTED. The spelling is LOOSE because the text says CUN; relabel to "CUN".
- **Legend "other codon" and "no aligned site":** SUPPORTED.
- **Caption "81 positions":** MISLEADING/LOOSE.
  - These 49 genomes have sites at 126 strict proline positions.
  - The figure code keeps a position only if at least 60% of the 49 genomes have a codon there and at least half of those are CCN/CUN. That leaves 81.
  - Because of that filter, the orange counts are 67/78/77/42 against the text's 78/82/82/44.
- **Caption "13 mitochondrial proteins":** WRONG, it is 12 (cox1 has two profile models).
- **Caption "≥ 90% … (weighted by family) have proline":** SUPPORTED.
- Suggested caption: "81 of the 126 positions in 12 mitochondrial proteins where ≥ 90% of surveyed mitogenomes (weighted by family) have proline, shown where ≥ 60% of these 49 genomes have an aligned codon."

### Known before
- **Table 3 = threonine in *Saccharomyces* and relatives:** SUPPORTED.
- **Ashbya sentence:** SUPPORTED. Ling 2014 confirms charging by the alanine enzyme and "strict avoidance of CUG and CUC codons".
  - LOOSE detail: Ling says "most likely derived from tRNA His". Fix: "read by a tRNA that most likely arose from tRNA-His".
- **"The four genomes here also avoid CUC and CUG":** SUPPORTED. CUC is never used; CUG is used 0, 0, 1 and 1 times. The blind re-check found CUG once in each.
- **"Proline is CCN in all 27 NCBI codes (tables 1–33, with gaps)":** SUPPORTED.
- **"No report of any other codon read as proline":** UNVERIFIABLE because it is search-based. My own searches found none either.
- **Sarhan et al. citation and accessions:** SUPPORTED. Sci Data 13:28, published online 2025-12-18, doi 10.1038/s41597-025-06447-x; GenBank TPA records BK072095–BK074789.
- **"annotated with the mold mitochondrial code (table 4…)":** MISLEADING/LOOSE.
  - The BK records have only a source feature.
  - The paper says it annotated each genome "using their corresponding genetic code", and NCBI Taxonomy gives these species code 4.
  - Fix: "deposited without gene annotation; NCBI Taxonomy assigns these species the mold mitochondrial code (table 4, CUN = leucine), and Sarhan et al. annotated each genome with its taxon's code".
- **Wolters 2023, 353 genomes, "assembled *W. infanticola* first (GenBank OR703884, identical to BK074532)":** SUPPORTED.
  - OR703884 was submitted 2023-10-19. It is the same 56,486 bp sequence, rotated by 36,482 bp.
  - Both come from the same reads (BioSample SAMN07511535 = SRR5938406), so they are not independent. The post doesn't claim they are.
- **"assigned it table 4 after running Codetta":** SUPPORTED. The methods say they used the yeast code when Codetta supported CUN = threonine and the mold code otherwise. I didn't open their Table S1; the audit did.
- **"the highest ratio … in the *Wickerhamiella/Starmerella* clade":** MISLEADING/LOOSE.
  - The paper says the highest values were in that clade "and … 'Other Dipodascales'".
  - Fix: "…the highest … in two Dipodascales groups, one of them the *Wickerhamiella/Starmerella* clade."
- **"If CUU/CUA are proline, some of those changes are silent…":** SUPPORTED as reasoning.

### Done
- **Census, "about 22,700 … (all 18,015 in RefSeq, …) with their declared codes":** WRONG in the parenthesis.
  - The README says 18,015 of 18,031 RefSeq records were scored; the merged set is 22,714.
  - Fix: "(18,015 of the 18,031 in RefSeq, plus …)".
  - "With their declared codes" is fine: the census uses each genome's declared table, which doesn't affect CUN.
- **"…profile models … The census was built to look for new stop codons":** SUPPORTED.
- **Codon contrast, "against all 151 Dipodascomycetes":** SUPPORTED.
- **"CUU … the four: p 2.5e-6 to 1.8e-8; the 147 others: lowest 0.049":** SUPPORTED. The four are 1.8e-8, 2.8e-8, 9.2e-7 and 2.5e-6; the lowest other is *Middelhovenomyces* BK073468.
- **Raw reads, "mapped reads from each assembly's source run":** SUPPORTED. The runs are SRR16988755, SRR5938406, SRR16989235 and ERR12422161, as named in the BK records.
- **"assemblies are faithful to their reads; the second sequencing is the independent one":** SUPPORTED.
- **Blinded re-check, "four controls (three from the same genera, plus *Starmerella aceti*)":** SUPPORTED. The controls were *W. azyma*, *W.* sp. BK074539, *S. castrensis* and *St. aceti*.
  - Note: only BK074539 is a strong control for CUN = leucine. The other three have only 4–14 CUN codons at conserved positions.
- **"told only that the yeast mitochondrial code was one candidate…":** MISLEADING/LOOSE.
  - The re-check report states: "Null hypothesis: NCBI table 4 (the candidate alternative was table 3)". Organism names were in the headers.
  - Fix: "It was told that table 4 was the default and table 3 (yeast mitochondrial) the alternative, which points at the CUN codons, but not that proline was in question."
- **"two methods (Codetta, and ortholog alignments to 11 relatives) … exactly these four and nothing else off table 4":** SUPPORTED.
- **"CUU … 75 of 76 … CUA … 298 of 300":** SUPPORTED.
  - CUU: 24 + 15 + 17 + 19 = 75 of 76.
  - CUA: 78 + 41 + 91 + 88 = 298 of 300.
- **Tree, "49 genomes, 4,403 aligned positions, CUN masked":** SUPPORTED.
- **"three *Wickerhamiella* form one branch (SH-aLRT support 100)":** SUPPORTED.
- **"*S. pinicola* deep inside *Sugiyamaella*, other 15 species use CCN; five nodes 98–100":** SUPPORTED (98.1, 99.9, 100, 99.5, 100).
- **"So the change most likely arose twice":** SUPPORTED as an inference.

### Not known
- **"Everything here is DNA":** SUPPORTED.
- **"I know of no such editing in fungal mitochondria":** UNVERIFIABLE (author's knowledge).
- **"ENA has only genome sequencing … no RNA-seq":** SUPPORTED. Every ENA run for the four species is genomic WGS.
- **"All four keep tRNA-Pro (UGG), tRNA-Leu (UAA) and a tRNA with anticodon UAG":** SUPPORTED by the tRNA scan (dipo.aragorn.tsv).
- **"That UAG tRNA is no more like tRNA-Pro…; related yeasts … have a similar one":** UNVERIFIABLE, because I didn't recompute sim_uag.py. A UAG tRNA is present in the control BK074539.
- **"Which tRNA and which enzyme … is open":** SUPPORTED.
- **"CUC and CUG are used 0–1 times per genome":** SUPPORTED for the protein genes. Optional: add "in the protein genes".
- **Next step: CBS 13331 or NRRL Y-17858, both in public collections:** SUPPORTED.
  - CBS 13331 is the type strain of *W. siamensis* (Khunnamwong 2014).
  - NRRL Y-17858 (= CBS 7922) is the type strain of *W. infanticola*.
- **"Peptides … would settle it":** SUPPORTED.

### Note for Ryan
- **"What this is":** SUPPORTED.
- **"Checked" list:** SUPPORTED.
  - At MAPQ ≥ 20: 35/40/2/20 sites per genome and 98.60/99.60/98.95/99.32% agreeing, matching "2–40 sites" and "98.6–99.6%".
  - "8 of 8 calls right" and the literature list match the README.
- **Audit item:** SUPPORTED. The pairwise method reproduces the result in all four; the same-genus control has 0 CUN at proline; the 81 vs 82 sites point is in the audit report.
- **Pending (no verdict):**
  - DS02 has four PacBio runs (SRR4787333, 334, 336 and 337), not only SRR4787334.
  - SRR16988751 is a second NRRL Y-17858 Illumina run with the same base count as SRR5938406. Treat it as a duplicate upload, not new data.
- **Risk, RNA editing:** SUPPORTED.
- **Risk, PYCC 7069 "both come from sugar-cane leaves in Thailand":** MISLEADING/LOOSE.
  - The sample record (SAMEA113572364) gives Thailand only, with no source, and I found no PYCC catalogue entry.
  - The species was described from six Thai sugar-cane-leaf strains.
  - Fix: "PYCC 7069 is from Thailand (sample record; no source or strain history given), like the type CBS 13331 and the five other strains in the species description, so it may be the type or a sister isolate under another number. Either way it is a separate culture and sequencing."
- **Risk, novelty and Codetta:** SUPPORTED. Wolters ran Codetta on *W. infanticola* but only used it to choose between tables 3 and 4.
- **Recommendation:** opinion, consistent with the findings.
