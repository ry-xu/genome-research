# Fact-check: fig_l02 image and the 15 Nov 383 -> 382 edits

Independent fact-checker (fresh agent), 2026-11-15 03:14. I did not make any of the checked items. Everything below
was recomputed with my own code (scripts and outputs in
`/private/tmp/claude-501/-Users-ryan-lab/40ff2367-a964-4475-b57b-539035c4acd3/scratchpad/figcheck/`):

- `recount.py` re-reads `out/sam/SRR22273278.sam.zst` and the two RNA SAMs with its own CIGAR walker. It applies the rules in the
  junctions.py docstring (identity ≥ 0.97, one gap ≥ 10 inside [P−15, P+D+15), 15 aligned bases and ≤ 1 mismatch on each
  side; through = no gap and crossing P or P+D by 15). It builds the templates (outer ends over all mates, clips added back),
  and writes per-fragment detail (`dna_frags.json`, `rna_frags.json`).
- `seqcheck.py` checks the genome sequence and GFF (frame, stop codon, translation).
- `placement.py` compares mismatches for each "shifted" junction read under the called and the predicted removal.
- `dnacmp.py`, `rnadraw.py` and `orient.py` compare the figure's inputs with the recount and measure unsequenced spans and
  duplicates.

P = 4487605, D = 31 throughout (0-based).

## 1. Numbers and labels in the image

| item | verdict | evidence |
|---|---|---|
| "RNA: 385 molecules lack the 31 bases" | OK | `rna_mol_sites.tsv` L02/0/joined: 497 fragments, 385 templates. The recount gives 497 fragments and 385 templates, and its 385 joined rows are identical to the table's (ends, orientation, removal). |
| "26 read through" | OK | Table: 32 fragments, 26 templates. Recount: 32 / 26, identical rows. |
| "DNA: 46 fragments read through, 0 lack them" | OK as fragments | `dna_sites.tsv` L02 joined 0, through 46. The recount gives 0 / 46, and its fragment set and ends equal `out/fig_l02_dna_templates.tsv` (46 rows, 0 differences). |
| Footnote "Each line is one molecule (read-pair outer ends)" | **WRONG for the DNA panel** | The 46 DNA lines are fragments. They have only **34 distinct outer ends** (10 groups of 2-3 identical templates, several with near-consecutive spot numbers, e.g. SRR22273278.2411825/.2411875/.2412063 and .6244303/.6244312). By the molecule rule used for the RNA (src/molecules.py), the DNA has 34 molecules. So the figure compares 385 RNA *molecules* with 46 DNA *fragments*. The 0-lack claim is unaffected. |
| Footnote "382 of the 385 lack exactly the predicted 31 bases" | OK as aligned, incomplete | 382 templates in `rna_mol_templates.tsv` have removal 4487605+31; the recount agrees. But see section 2: the 383rd molecule has one read at exactly 4487605+31. |
| Footnote "3 lack 31 bases shifted by 3-4" | **WRONG / overstated** | As aligned: 2 at 4487601 (−4), 1 at 4487608 (+3); 4487605-4487607 are equivalent placements, so 4487608 is 1 base past the equivalent range. But (a) the "+3" molecule (SRR27313660.5352211) has mate 1 at exactly 4487605+31 (`21=31N124=1X4=`); only mate 2 is at 4487608. (b) Every one of the 3 shifts rests on a tie or on a single low-quality base. Mismatches per read, called vs predicted placement: SRR27313660.4934966.2 is 3 vs 3 (a tie); SRR27313660.5352211.2 is 1 vs 2, the deciding base Q11, while its mate favours the predicted placement on a Q37 base; SRR32148620.9948378.1 is 0 vs 1, the deciding base Q12. The reads do not show that these 3 molecules have a different removal. |
| "RNA: SRR27313660 and SRR32148620, first 10 M read pairs each" | OK | `map_culture.sh`: `fastq-dump --split-spot -X 10000000`; `map.log` / `.n`: 20,000,000 reads each. |
| "neither library is stranded" | OK (matches post/RESULT) | RESULT.md: mate-1 sense share 0.62 / 0.72 (rRNA), 0.51 / 0.56 (CDS), against the 90 % rule. |
| "DNA: SRR22273278, first 8 M read pairs" | OK | `-X 8000000`; `map.log`: 16,000,000 reads. |
| Ticks 4,487,405 / 505 / 605 / 736 / 836 | OK | = P−200, P−100, P, P+D+100, P+D+200 on a linear axis. Their pixel positions match (e.g. 4487736 at 1411.7 px of 2414 predicted, drawn at about 1411). The uneven gap between 605 and 736 is real (131 bases). |
| 31-base band at 4487605-4487636 (0-based, half-open) | OK | `Rectangle((P, ·), D)`. Pixel edges are 987 and 1088 predicted (of 2414), drawn at about 988 and 1088. The tick at 4,487,605 sits at the band's left edge, consistent with "0-based". |
| UYP47557 extent | OK | GFF CDS 4487084..4487611 (+, 1-based) = [4487083, 4487611) 0-based, as in `GENES`. Drawn from the window edge (4487345). |
| UYP47558 extent | OK | GFF 4487703..4488440 (+) = [4487702, 4488440). Drawn from 4487702 (about 1301 px predicted and drawn) to the window edge, with no head because it runs past the window. |
| UYP47557 drawn end | Minor, misleading | The arrow tip (4487601-4487611) is under the orange box (zorder 3, from 4487605), so the gene *looks* as if it ends where the band starts. It really ends 6 bases inside the band, with its TAA stop at 4487608-4487611. The subtitle states this, but the drawing hides it. |
| "Ca. Lokiarchaeum ossiferum B35, CP104013.1" | OK by consistency | The GFF gives CP104013.1, taxon 2951803, isolate Loki-B35 and old-name "Candidatus Lokiarchaeum sp. B-35". The current species name is not in any local file (no network used); it matches the post. |
| "Positions 0-based on CP104013.1" | OK | |

## 2. Is 382 right, and is the new bullet's explanation right?

- **382 is right under one definition.** 382 molecules have every junction read at exactly 4487605+31; 383 have at least one
  read there. The old 383 was not simply an error. It was the "at least one read" count.
- **The explanation in the new bullet ("383 came from the genome-wide gap table, which has no flank filter") is WRONG
  as a cause.** The gap table (`rna_mol_gaps.tsv`) does have no flank filter, but here that filter removes nothing.
  All 495 fragments in the gap table's 4487605+31 row pass every site rule and are site-joined. Its 383 templates are
  the 382 plus SRR27313660.5352211.
  - In that fragment, mate 1 shows 4487605+31 and mate 2 shows 4487608+31, and both pass the flank rules.
  - The gap table counts the fragment under both removals. Its 4487608+31 row has 1 template, the same molecule.
  - The site table keeps one removal per fragment, the last joined read seen, because
    `at_site[(site, sh, "joined")][fr] = ...` overwrites. Here that is mate 2, so the molecule appears as "+3".
  - `junctions.py` does the same at fragment level (494 at 4487605, "4487608+31:1").
- New bullet, "the other 3 are shifted by 4 or 3 bases": **WRONG for 1 of 3, and weak for all 3** (section 1).
- New bullet, "draws the same molecules the text counts: the 411 RNA molecules … and the 46 DNA fragments re-read …
  with the same rules": OK (verified row-for-row). "The script stops unless the counts are 385 / 26 / 0 / 46": OK
  (`assert`; the script is run without `-O`). "Corrected above and in RESULT.md": OK. Both places were edited
  (03:04:45), and the image (03:05:32) was written after the script's last edit (03:05:08).
- RESULT.md line 138, "382 of 385; 383 in the unfiltered gap table": the number is OK. **"unfiltered" gives the wrong
  reason** (see above).
- RESULT.md line 197, "382 at 4487605+31 (was 383, corrected 15 Nov)": the number is OK as an every-read count. "corrected"
  implies 383 was wrong; it was the at-least-one-read count.
- Related, not edited today: RESULT.md line 142, "185 of 188 share the predicted removal" (message strand). The mixed molecule is
  on the sense side (gap table `templates_plus` 186), so it is 186 of 188 by the at-least-one-read count.
- Post line 9, "382 of them lack exactly the predicted 31 bases": defensible, but the honest summary is stronger. All 385 fit
  the predicted removal: 382 as aligned, and the other 3 within a tie or one low-quality base (Q11, Q12), one of them
  with its other mate exactly at the predicted place.

## 3. Are the DNA fragments classified like src/molecules.py / src/junctions.py?

**OK for classification.**
- Pass 1 of `fig_l02.py` repeats the joined/through conditions of molecules.py character for character (same `parse`,
  MIN_ID, FLANK, FLANK_MM).
- "A fragment with a joined mate counts as joined" is reproduced: `cls[fr] = "joined"` overrides `setdefault(…, "through")`.
- Templates use all mates with clips added back. The reference has one contig, so the contig filter changes nothing.
- My independent recount gives 0 joined / 46 through and the same 46 fragments with the same ends.
- **Difference in unit:** molecules.py merges fragments with identical templates into molecules; the figure does not
  for DNA (46 fragments, 34 molecules). RESULT.md line 141 ("Counts are molecules (templates)") sits under a table
  whose DNA column (0 / 46) is fragments. This was not edited today, but it has the same unit problem.

## 4. Title and subtitle against the post

- "31 bases cut out of a protein gene's RNA, not its DNA": **OK against the post's lead**, which says "It is cut out
  at the place and length I predicted". Both the lead and the title are stronger than the post's own "Not known" and
  "Recommend" sections: endonuclease assumed, not shown; an RT skip at L02 itself not tested; "stated as an observation".
  "missing from" would match the evidence. Optional.
- "the shaded 31 bases hold the first gene's stop codon": **OK.** UYP47557's last codon is TAA at 4487608-4487611,
  inside [4487605, 4487636).
- "Removing them joins the two genes in frame": **OK.** (P − 4487083) = 522 = 174 codons, so the cut falls between
  codons. (4487702 − 4487636) = 66 = 22 codons. The spliced 1326-nt ORF has its first stop at codon 442, giving 441 aa,
  and reads …RPLLNNMI|KPCTGH… across the junction, as in the post.

## 5. The drawing

- **Solid lines over unsequenced middles:** minor and disclosed ("read-pair outer ends").
  - RNA joined: 98 of 497 fragments have a gap between mates, 1.9 % of the drawn length, none over the band.
  - RNA through: 5 of 32 fragments, 1.0 % of the drawn length; 1 fragment's unsequenced middle overlaps the band.
  - DNA: 3 of 46 fragments have their unsequenced middle overlapping the band. They count as "through" because one
    mate reads ≥ 15 bases into the stretch from one edge, but the line is drawn solid across all 31.
- **Clipping:** disclosed. 17 DNA templates start before the window and 9 end after it.
- **Colours:** consistent. Orange means the stretch is missing; grey means read through, in both RNA and DNA. The two
  "+3"/"−4" stubs at the band edges show aligner placements that the reads do not support (section 1).
- **Gene track:** UYP47557's end and stop codon are hidden under the orange box (section 1).
- **RNA molecule count, a caveat outside today's edits.** The molecule key includes an inferred transcript strand, which
  is arbitrary for these unstranded libraries.
  - Ignoring orientation gives 364 joined and 24 through, not 385 / 26. In 20 of the 23 merged pairs both orientations
    come from the same run.
  - Conversely, 5 same-orientation duplicate groups span both runs, so they cannot be one molecule.
  - So "molecule" counts are approximate in both directions. The share is unchanged (364/388 = 93.8 % against
    385/411 = 93.7 %).

## Required fixes (most important first)

1. **New "The image" bullet in the post, and RESULT.md line 138: fix the stated cause.** The 383 was not caused by the
   gap table lacking a flank filter; that filter removes nothing here. One molecule (SRR27313660.5352211) has one mate
   at exactly 4487605+31 and the other at 4487608+31. The gap table counts it under both, and molecules.py keeps the
   last mate's removal. Say that 383 = molecules with at least one read at the predicted removal and 382 = all junction
   reads there. Drop "unfiltered" from RESULT line 138 and "error"/"corrected" framing where it implies 383 was false.
2. **Image footnote "3 lack 31 bases shifted by 3-4"** (and the bullet's "the other 3 are shifted by 4 or 3 bases"):
   rewrite. Suggested: "382 lack exactly the predicted 31 bases in every read; in 3 the aligner placed the gap 3-4
   bases away, on a tie or one low-quality base (one of them has its other read at the predicted place)."
3. **Image footnote "Each line is one molecule": wrong for DNA.**
   - Either deduplicate the DNA panel by outer ends, giving 34 molecules, and label it "DNA: 34 molecules (46 fragments)
     read through, 0 lack them".
   - Or reword: "RNA lines are molecules; DNA lines are fragments (46 fragments, 34 distinct)".
   - If deduplicating, change the script's assert, and make RESULT.md line 141's "Counts are molecules" exclude the DNA
     column or give 34.
4. (Minor) Draw UYP47557 above the orange box or mark its stop codon, so the gene visibly ends 6 bases inside the band.
5. (Optional) The title's "cut out of" could read "missing from", to match the post's "Not known" / "stated as an
   observation". The post's own lead uses "cut out", so this is not a mismatch between image and text.
