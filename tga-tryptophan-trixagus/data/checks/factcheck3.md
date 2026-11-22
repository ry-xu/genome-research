# FACTCHECK3: revised sentences of post.md

Checked 20:02 2026-11-11. Scope: only the sentences revised after FACTCHECK2. Post line numbers refer to
`scratchpad/post/post.md`. Project root abbreviated `P = lab/work/qc-rejected-codes`; literature
extracts in `scratchpad/factcheck2/`. One new analysis, on reads FACTCHECK2 had already retrieved:
`factcheck3/ttur_sbost_reads.fa`, `ttur_sbost_all.tsv` (built with `pull.py`, `gc.py`).

9 non-OK items: 2 WRONG, 3 OVERSTATED, 1 UNSUPPORTED, 3 UNCLEAR. Most serious first.

## Non-OK

### 1. WRONG: "none matches a Shikimatogenerans RF2" (l.14); "none found" (Note l.43)
- **Evidence.** The RF2 database of the read scan contains five Shikimatogenerans RF2s, including the bostrichid
  symbiont's (`S_bostrichidophilus|0_75`), plus E. coli (`P/results/rf2_readscan/rf2_scan_summary.py`). Best-hit
  recount: **Ttur 125 reads and Tder 6 reads best-hit the S. bostrichidophilus RF2**. That makes the
  sentence literally false. README l.643-647 counts only reads "best-hitting a throscid RF2", which is correct.
  README l.649-651 and FACTCHECK2 l.84-85 do not mention the 125 reads.
- **Are they the symbiont's?** No. The 125 Ttur reads are 78 unique sequences with median GC 0.245 (0.185-0.305).
  The Ttur genome's GC is 0.127 (`P/results/symbiont_profile.tsv`). Re-searched with `--max-target-seqs 0`, they
  hit only the bostrichid RF2: they hit no throscid RF2 and not E. coli. A Ttur-symbiont RF2 would sit closest to
  the throscid copies. Only 5 of the 125 are in FACTCHECK2's 275 low-GC set. The 6 Tder reads have GC 0.39-0.51.
  So the deletion stands: there are no throscid-RF2 reads, and the gene's neighbours lie 113/124 bp apart. But
  125 is exactly Ttur's 10th-percentile expectation (README l.645: 177 / 125). Anyone who reruns the scan will
  see "125 Shikimatogenerans RF2 reads" where the post says none.
- **Fix.** Post: "none matches the RF2 of any throscid symbiont, where 27-206 would be expected ... In *T.
  turgidus*, 125 reads match the bostrichid symbiont's RF2 best, but at twice the genome's GC (25% vs 13%)
  they come from another bacterium." Note l.43: "(none from a throscid symbiont)". Add these reads to README
  l.649-651.

### 2. WRONG: "gTranslate 2026 also found standard-code genomes nested among TGA readers, in *Stammera* and in a Minisyncoccia genus ... the *Trixagus* symbionts are another case, in beetle symbionts." (l.26)
- **Evidence.** gtrans.txt l.92: the five code-11 *S. capleta* "form a monophyletic clade basal to all other
  Ca. S. capleta ... reassignment may have occurred after the basal clade diverged". That is a basal
  standard-code lineage, consistent with a single change, not nesting. gtrans.txt l.7 (abstract) also says
  basal. Only the Minisyncoccia genus GCA-2747955 is nested (gtrans.txt l.132-134: code-11 genomes in a lineage
  otherwise of code-4). *Stammera* is a symbiont of tortoise leaf beetles (gtrans.txt l.10), so "in beetle
  symbionts" is not a new setting: a code split inside a beetle symbiont is already published. README l.417-418
  has the same "nested ... in *Stammera*" error. README l.422 has it right ("basal code-11 lineage").
- **Fix.** "gTranslate 2026 found standard-code genomes nested among TGA readers in a Minisyncoccia genus, and
  a basal standard-code lineage in *Stammera*, a tortoise-beetle symbiont. So a patchy code within one group is
  not new; the *Trixagus* symbionts are another case." Drop "in beetle symbionts". Correct README l.417-418.

### 3. OVERSTATED: "also under a second model that lets amino-acid preferences vary across sites (LG+C20)" (l.16); "both models prefer" (l.35)
- **Evidence.** The AU numbers are right. `P/results/throscidae_tree/c20au/test.iqtree` l.133-137 gives
  constrained-readers p-AU 2e-53 (Δ3503.5). But the fitted mixture puts weight 0.929 on the single LG+F
  (empirical-frequency) component, and each of the 20 C20 profiles gets at most 0.014 (`c20au/test.iqtree`
  l.32-55; `c20/concat.faa.log` l.49). The fitted model barely lets preferences vary across sites. It is close
  to the first model, LG+F+G4 (`concat.faa.iqtree` l.39). So it is a weak independent check on the
  compositional bias that Note l.49 flags (6 of 8 sequences fail the composition test). The sentence implies
  more robustness than the fit delivers.
- **Fix.** "... also under LG+C20, though that fit put 93% of sites on one frequency profile, so it differs
  little from the first model." Alternatively, drop the C20 clause from the post and let the identity check
  carry it, as Note l.49 already does.

### 4. UNSUPPORTED (LG half): "neither rejects it as sister to the other *Trixagus*" (l.35)
- **Evidence.** Under LG, that tree was never AU-tested. `au/trees.tre` holds only 2 trees (ML and
  readers-constrained; `au/test.iqtree` l.108-111). The only LG evidence on the *T. turgidus* placement is
  UFBoot 93 on the node that excludes it (`concat.faa.treefile`). That is moderately strong support against the
  alternative, not a non-rejection. Only C20 tested the alternative: p-AU 0.109, ΔlogL 46.0
  (`c20au/test.iqtree` l.135).
- **Fix.** "C20's AU test does not reject it as sister to the other *Trixagus* (p = 0.11), and under LG the
  split that excludes it has 93% bootstrap support." Or run the LG AU test with the third tree before posting.

### 5. UNCLEAR: "C20's partial search kept the LG placement" (Note l.52)
- **Evidence.** `c20/concat.faa.log` l.48 shows a starting logL of -310811.498. That matches the alternative
  tree (tree 2 in `c20au/test.iqtree` l.135, -310811.40), so the search seems to have started from the
  alternative, not the LG placement. Iteration 1 is -310765.610 (l.71) and iteration 4 is -310765.122 (l.74).
  Those match the LG ML tree's -310765.39. The search therefore moved to a tree with the LG tree's likelihood;
  it did not "keep" a placement. The tree was not saved, so the topology is inferred from logL only.
- **Fix.** "A partial C20 search started from the alternative and reached the LG tree's likelihood (tree not
  saved), but the alternative is not rejected (AU p 0.11)."

### 6. OVERSTATED (minor): "in the third, the assembly of that stretch is not supported by reads" (l.14)
- **Evidence.** README l.652-656: the Tder junction is 514 bp (39,333-39,846). Of that, 298 bp
  (39,457-39,754) has no reads. The 10-bp assembly gap is spanned by 9 reads. Note l.43 states this correctly.
- **Fix.** "in the third, 298 bp of that stretch has no read support, so its length is uncertain."

### 7. OVERSTATED (minor): "a second fresh agent checked this draft and found 20 more. All are fixed here." (Note l.45)
- **Evidence.** FACTCHECK2 checked the previous draft, not this one. Not yet fixed:
  - The post still embeds `image.png` (l.20), but the post directory holds only post.md and
    commons_correction.md (FACTCHECK2 #25). The figure is `P/results/throscidae_tree/figure.png`, and its
    captions match the post.
  - The README still has unfilled `@@T@@` placeholders at l.21, 416, 492, 588, 590, 619, 625 and 633.
  - The README calls the *Nardonella* assembly "7-contig" at l.293 and l.398, where the post says 26 contigs
    in 7 scaffolds.
  - The *Stammera* "nested" wording (item 2) is in both the post and the README.
- **Fix.** "checked the previous draft and found 20 more". Then either fix the items above or drop "All are
  fixed here".

### 8. UNCLEAR (minor): "the largest earlier one (Shulgina & Eddy 2021) took about 250,000 GenBank assemblies" (l.3)
- **Evidence.** The number is right: se.txt l.1, 3 and 183 say "over 250,000". But gTranslate 2026 ran on GTDB
  R09 (596,861 genomes, 356,020 after deduplication) and GlobDB, with 250,583 test genomes (gtrans.txt l.17,
  61, 92, 132, 139). "Largest" is true only if "earlier" means "before the GTDB-based screens", and a reader
  may take it as "largest screen to date".
- **Fix.** "an earlier screen not tied to GTDB (Shulgina & Eddy 2021) took about 250,000 GenBank assemblies".

### 9. UNCLEAR (minor): "these do not, like *Zinderia*, *Stammera* and *Providencia*" (l.32)
- **Evidence.** The content is right: choi.txt l.82 gives the trnW anticodon of *Zinderia*, *Stammera* and
  *Providencia* as CCA, and Parks l.593 gives a UCA tRNA in Eggerthellaceae. But the comma construction reads
  as "unlike them".
- **Fix.** "these do not; nor do *Zinderia*, *Stammera* or *Providencia* (Choi et al. 2026)."

## OK (brief)

| # | Claim (post line) | Evidence |
|---|---|---|
| 1 | Intro, l.3: usually lost RF2; recent screens use GTDB; GTDB tools switch to code 4 when coding density clearly rises; "predates the genomes below" | gtrans.txt l.12 (CheckM rule in GTDB-Tk/CheckM2/SingleM); parks2026.txt l.6; reads SRR2803xxxx (2024) vs 2021 |
| 2 | 103 list, l.7 | `P/results/labels.tsv`: Hodgkinia 49, Nasuia 24, Mollicutes orders 17, code-25 Gracilibacteria/Absconditabacterales 12, Pinguicoccus 1 = 103; CHECK 38 |
| 3 | Codetta re-check, l.13 | README l.337-357 (corrected 18:11): 15 accessions, 3 readers + Nardonella + 2 controls called, 9 not |
| 4 | "first 4 million read pairs", "27-206 expected", "66" | rf2_readscan scripts `fastq-dump -X 4000000 --split-spot`; README l.643-647 |
| 4 | "113 and 124 bp"; every relative's RF2 lost SPF; two broken at the front; bostrichid keeps SPF | README l.463-497 (SPF to F--, F--, S--, SSI; Tduv 5' end gone; Tmey in-frame TGA at codon 129; S_bost SPF) |
| 5 | Tser sister to Tduv; last sentence of the at-least-twice bullet (parsimony, under 5%) | `concat.faa.treefile`; README l.548-632; usage bullet l.10 (1.4-4.8%) |
| 6 | GTDB, l.18: code is not why GTDB rejected them; Parks note for other TGA readers | README l.658-663; parks2026.txt l.627 |
| 7 | *Nardonella* paragraph, l.22 | README l.292-327 (26 contigs / 7 scaffolds, 25 of 141 TGA sites covered, 72% / 70%) |
| 8 | TGA sense-codon list; Choi six; only *Cerorastrobacter* in Bacteroidota; ciliate/dinoflagellate symbionts | choi.txt Table 1 (l.55-82) and title (l.1); gtrans.txt l.10 |
| 8 | Parks: same argument for Eggerthellaceae, with a lineage that still stops at TGA between the readers | parks2026.txt l.6, l.471, l.593, l.597 |
| 9 | Done, l.29-30 | README l.232 (criteria 16:19), l.255-263; ortholog/read files 16:20-16:25; AU rule after the tree (see Note) |
| 10 | How TGA is read (content); order of events; "either placement makes a third loss" | choi.txt l.82; parks2026.txt l.593; tree topology |
| 11 | Note times 16:19 / 17:03 / 17:27 / 18:15; C20 tested the LG-constrained tree | README l.232; treefile mtime 17:03; au/CRITERION.txt 17:27; c20au/CRITERION.txt 18:15 (note appended 19:07, after the 18:15:37-18:34:58 run); c20au/test.iqtree |
| 11 | Prior work; throscidiphilus risk; at-least-twice risk; Recommend matches the brief's bar | README l.412-425; brief-codon-free.md l.9 |
