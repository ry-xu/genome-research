# Tweet check: mito-proline thread (21 Nov 2026)

Thread: `~/.local/share/lab-dashboard-v3/posts/mito-proline/thread.md`. Figure: `~/.local/share/lab-dashboard-v3/figs/mito-proline/public/proline_rows.png`,
script `work/curation/src/fig_proline_rows.py`. Line numbers below are for the files as of 21 Nov 12:15. `post.md` = this folder's post.md.
"## 1" is 222 words after the recheck fixes (217 before) without the image line (limit 280).

## Tweet 1, line by line

The tables quote draft v1. Four lines were reworded after the fresh-agent recheck; the final wording and the sources for it are in the last section.

| Thread text | Source | Verdict |
|---|---|---|
| "four yeasts read the codons CUU and CUA as proline in their mitochondria" | post.md:1 ("In the mitochondria of four yeasts, the codons CUU and CUA mean proline"); post.md:5 | OK. Inferred from DNA; caveat in tweet 2 (post.md:24). |
| "the first genetic code known to read anything but CCN as proline" | post.md:1 ("No known genetic code reads any codon other than CCN ... as proline"); post.md:15 ("Proline is CCN in all 27 NCBI codes ... no genetic code in any organism or organelle is known to read another codon as proline"); pre-post check fix, post.md:45 (Streptomyces is mistranslation, not a code); Saiph's 20 Nov prior-art pass found nothing new (work/proline-codon-yeasts/notes/prior_art_github.txt, empty hits; agents/saiph/STATE.md:104) | OK, supported by the checked prior-art search. Not weakened. Rests on searches (post.md Note, "Risks": someone may have seen it and not published). |
| "Genes are read three letters (a codon) at a time, and each codon names one amino acid." | post.md:3 | OK (simplified: 3 of 64 are stops). |
| "Proline is spelled CCN (CCU, CCC, CCA or CCG) in every known genetic code." | post.md:1, post.md:15 | OK. |
| "CUU and CUA normally mean leucine." | post.md:1 ("not leucine"), post.md:3 (CUN leucine in the standard code) | OK. |
| species list W. siamensis, W. infanticola, W. verensis, Sugiyamaella pinicola | post.md:5 | OK. |
| "CUU/CUA sit where related proteins have proline: 286 of 290 times at strongly conserved positions" | post.md:5 ("of the four genomes' CUU/CUA codons at strongly conserved positions, 286 of 290 are at proline positions"); post.md:47 (recount, 3 Phe, 1 Ala); 286 = 44+78+82+82 lines in work/new-stops-euk/results/cun_pro/audit/BK0742{62,532,538,542}.1.census_ctn_p.txt | OK. |
| "These yeasts spell leucine UUA/UUG instead." | post.md:5 ("Leucine is spelled UUA/UUG, as in their relatives"); post.md:6 (UUA/UUG 140 to 166 times) | OK. |
| "the other 17 Wickerhamiella use CCN at 80–85 such proline positions, and CUU/CUC/CUA/CUG at 0–1" | post.md:7 ("CCN at 80–85 such proline positions and CUN at 0–1") | OK (CUN spelled out). |
| "It most likely happened twice: the three Wickerhamiella form one branch, and S. pinicola sits inside Sugiyamaella, whose other 15 species use CCN." | post.md:22 ("support 100%"; "deep inside Sugiyamaella, whose other 15 species here use CCN"; "most likely arose twice") | OK. "15 species" = those in the 49-genome tree. S. pinicola is the thinner origin (post.md:26), not stated in the tweet. |
| "raw reads carry the codon at all 286 sites (98.9–99.6% agree)" | post.md:8 | OK. |
| "another lab's W. siamensis culture has it at all 82 sites" | post.md:9 (PYCC 7069, 1,432 of 1,450 reads) | OK. |
| Prior work: "these codons had already changed meaning twice: to threonine in baker's yeast, and to alanine in Ashbya gossypii" | post.md:5 ("already changed meaning twice"; post.md:3 says the same for the four CUN codons), post.md:14 (all four CUN = Thr in Saccharomyces; CUU/CUA = Ala in Ashbya, Ling et al. 2014) | OK. Baker's yeast: all four CUN; Ashbya: CUU and CUA. |
| "Some Streptomyces bacteria occasionally put proline at other codons, but as mistranslation, not a code." | post.md:15 ("partial mistranslation ... sometimes put proline at alanine, asparagine or threonine codons", Vargas-Rodriguez 2021, Schuntermann 2023) | OK. |

## Tweet 2

| Thread text | Source | Verdict |
|---|---|---|
| "everything here is DNA" | post.md:24 | OK. |
| "There is no RNA or protein from these four yeasts yet" | work/proline-codon-yeasts/RESULT.md:35 (R0: none in SRA, ENA incl. strain searches, PRIDE, ProteomeXchange); outbox/2026-11-20-proline-codon-rna-request/request.md:8 (NCBI STAT and Logan not searched) | OK as "no public data found"; "yet" keeps it from claiming none exists anywhere. |
| "RNA editing (U turned into C at each of these codons) is not fully excluded" | post.md:24; request.md:7 | OK. |
| "In the one relative with enough mitochondrial RNA in public data (W. versatilis)" | RESULT.md:42-48 (R2: four other runs not testable; SRR30194995, W. versatilis PYCC 7102, testable); request.md:9, :18 | OK. |
| "we found no U-to-C editing at the 4,456 coding positions it covers well" | RESULT.md:48 (4,456 coding positions at depth >= 20; 5 U>C, all fixed in every read, RESULT.md:49-52); notes/audit_R/rdd_with_ont_filter.txt:2,6 (same-strain ONT DNA SRR38099989 removes all 5 U>C: "U>C count 0"); notes/audit_R/ont_check.txt rows 12142, 14947, 16794, 16914, 21582 ("DNA variant (strain)"); Sabik's check (work/referee/sabik/proline_req/CHECK.md, "The ONT DNA filter removes all 5"; commons 2026-11-20-thirteen.md:1368-1369) | OK, scoped. Caveats not in the tweet: the audit's relaxed-mapping step (d) never finished (agents/saiph/STATE.md:117, 127, 132); only 25% of coding nt reach depth 20 (notes/audit_R/coverage_SRR30194995.txt); W. versatilis itself spells proline CCN (84 CCN, 0 CUN at strict Pro columns, work/new-stops-euk/results/cun_pro/pro_leu_wss.tsv:53), so this tests editing in the family, not in the four. "It covers well" carries the coverage scope. |

## Figure

Cox2, census profile `mito_Fun_cox2` columns 100-108, the census alignment (best hmmsearch domain of each gbmitoorf ORF,
work/new-stops-euk/data/gbmitoorf.genes.tsv.gz + data/mito_profiles.hmm; posterior '*' at all 27 codons). Each 27-nt
string checked as a unique forward-strand match in work/proline-codon-yeasts/data/genomes/all.fa:

| Row | Accession | First base | DNA (27 nt) | Translation (table 4; CUA = P in the four) |
|---|---|---|---|---|
| W. brachini (relative) | BK074528.1 | 6574 | TTAATAGCATTCCCATCATTTATATTA | L I A F P S F I L |
| W. siamensis | BK074538.1 | 31423 | TTAATCGCATTCCTATCATTTATCTTA | L I A F P S F I L |
| W. infanticola | BK074532.1 | 10766 | TTAATTGCATTCCTATCATTTATATTA | L I A F P S F I L |
| W. verensis | BK074542.1 | 22393 | TTAATTGCATTCCTATCATTTATATTA | L I A F P S F I L |
| S. americana (relative) | BK074251.1 | 238 | TTAATAGCTTTCCCATCATTTGTATTA | L I A F P S F V L |
| S. pinicola | BK074262.1 | 1086 | TTAATAGCTTTCCTATCATTTGTATTA | L I A F P S F V L |

- Column 104 (the orange CUA) is in each target's census set: `CTA mito_Fun_cox2:104` at audit/BK074262.1.census_ctn_p.txt:8,
  BK074532.1:11, BK074542.1:20, BK074538.1:21. Reads: work/new-stops-euk/results/cun_pro/reads/BK074538.1.tsv:42 (697/698 CTA),
  BK074532.1.tsv:58 (47/47), BK074542.1.tsv:34 (183/184), BK074262.1.tsv:82 (225/226); consensus P, share 1.00, strict.
- Other relatives at column 104 (same extraction): W. sergipiensis CCA, W. spandovensis CCA, S. grinbergsii CCT, S. japonica CCA,
  Starmerella bombicola CCT, so the two relatives shown are typical.
- Relatives were chosen by the tree: W. brachini is one of Saiph's three sisters of the target clade (RESULT.md "sisters
  BK074528/37/40"); S. americana is next to S. pinicola in the old figure's tree (image.png).
- Caption number "286 of 290" = post.md:5. Caption "Orange: CUA, leucine in the standard code, where relatives have CCN (proline)" matches the figure.
- No title, white background, orange only on the four CUA codons and their P; one pale column band; names plus a grey "relative" sublabel.

## Recheck by a fresh agent (21 Nov 12:19)

A fresh general-purpose subagent got the thread, the figure, its script and the raw project files (post.md,
work/proline-codon-yeasts, Saiph STATE, request.md, Sabik CHECK.md, commons, work/new-stops-euk/results/cun_pro), but not
this file. It edited nothing. Its findings: no number in the thread was wrong. Figure correct: all six 27-nt stretches match
all.fa at the stated first base, forward strand, and are unique on either strand. Translating ~110 codons either side gives
clean Cox2 with no stops (YQWYWSYEYSDF motif). The orange codon is CTA at mito_Fun_cox2:104 in all four (census_ctn_p.txt;
reads codon indices 119/100/117/119). Table-4 letters are right, names match the headers, and 286 of 290 is exact. No layout
problem. Its check script: scratchpad figchk.py (session scratch, not kept).

Applied to the thread:
- ## 2 rewritten. "There is no RNA or protein" became "No RNA or protein ... is in public databases yet" (STAT and Logan were
  not searched, RESULT.md:35-36). "is not fully excluded" became "can't be ruled out" (post.md:24). The relative test is
  scoped: "(which spells proline CCN)", which is 84 CCN and 0 CUN at strict Pro columns (pro_leu_wss.tsv:53; RESULT.md:16).
  "it covers well" became "read 20+ times" (4,456 of 17,526 coding nt; coverage_SRR30194995.txt). It adds "but that doesn't
  test the four themselves" (request.md:10).
- "These yeasts spell leucine UUA/UUG instead." became "Leucine is spelled UUA/UUG, as in their relatives." "Instead"
  implied the four switched their leucine codons, but relatives already use UUA/UUG (post.md:5, :47).
- "each codon names one amino acid" became "... (or a stop)".
- "whose other 15 species use CCN" became "where 15 other species use CCN". One of the 15 is Sugiyamaella sp. yHKS617.

Not applied (optional, for Ryan):
- "appear to read" in the opening line: kept "read", with the caveat in ## 2, matching post.md:1.
- "first genetic codes" (plural, two origins): kept the singular, which states the claim as a class.
- The RNA result is preliminary. Saiph's Arm R audit stopped before its relaxed-mapping step (d), RESULT.md's verdict is
  "TBD", and Saiph's 11:38 post says "nothing is claimed here yet" (commons 2026-11-20-thirteen.md:791). Sabik checked the
  5 U>C and the ONT filter (CHECK.md), but only as part of the request's wording.
