# Kill-check A5: 23-nt removal in *M. petrocarbonis* (GCA_030522375.1:y0, CP130488.1:1135031-1135054)

Scripts (`*.py`) and outputs (`*.out`, `*.tsv`, `*.paf`, `sra_meta/`) are in this folder. Every number below comes from one of them, named in brackets.
How it was done:
- Raw SAM lines for the y0 window were pulled from the four `RUN.sam.gz` files (`y0win.RUN.sam`) and classified with my own code (`sampass.py`, `templates.py`, `tally.py`). `reads.tsv` was not used.
- All reads in the window were then realigned to the whole genome (`win_reads.py`, minimap2 -ax splice:sr), giving `y0win_wholegenome_splice.sam`.

**1. Recount: SURVIVES (5 joined / 5 through), but the "removed" call sits exactly on the threshold.**
- Per run, after collapsing duplicates on both mates' coordinates [tally.out]:
  - AM31: 0 joined / 4 through.
  - AM32: 2 / 1. The raw count was 3 joined: SRR28822790.2442967 and .2443108 have identical coordinates on both mates and nearby spot numbers, so they are one molecule.
  - AM33: 3 / 0.
  - Metagenome: 0 / 44 (51 reads before collapsing).
- Allowing a 10-nt flank instead of 20 adds no RNA templates.
- The whole-genome realignment gives the same 5/5 [site63_dedup.out].
- The share is 0.50, exactly the "removed" bar. Losing any one joined template gives 4/9 = 0.44, which is "partly removed". "Spliced" (joined ≥ 2) is robust.
- AM31 vs AM32+AM33: Fisher two-sided p = 0.048; Freeman-Halton exact test over the 3 runs p = 0.048 [tally.out]. Most of this split comes from orientation (check 3):
  - 3 of AM31's 4 through templates are in the minority orientation, the one gDNA would give.
  - With only majority (dUTP-like) templates, AM31 0/1 vs the others 4/1 gives p = 0.33 [extra_stats.out].

**2. Read-level: SURVIVES.**
- Placement: there are 6 joined alignments, all with exactly a 23-nt gap at window position 383. Only one placement is possible (0 nt microhomology) [getwin.out, gaps_all.tsv].
- Alignment quality: MAPQ is 55-60. No soft clips at the junction (only 4S at a read end), no SA tags, and the shorter flank is 21-54 nt [y0_templates.tsv, joined_tags.out].
- SRR28822789.2494549:
  - All 4 mismatches in mate 1 and all 10 in mate 2 fall on Q11 bases.
  - By chance that would have probability 4.6e-4 and 9.1e-6.
  - All 18 mismatches in the joined templates are at Q11, and 15 of them are C>A or G>T [extra_stats.out, qualcheck.out].
  - These are sequencing errors, not another strain. Ignoring them, the read is the M. petrocarbonis joined form.
- Whole-genome realignment against the whole genome plus all 14 Methanomethylicales census genomes:
  - minimap2 -x sr (no splicing) and blastn (pid ≥ 97, no gap): 0 of 6 joined reads align contiguously across the junction ± 20 nt [wholegenome.out, map_sr.paf, blast_rna_templates.tsv].
  - The best hit is always CP130488.1, with a 23-nt deletion.
- Reading frame: removing the 23 nt (23 mod 3 = 2) joins the two ORFs into one of 364 codons [frame.out].

**3. Library behaviour: SURVIVES. The libraries are dUTP-stranded with about 9% gDNA-like reads, and there is a measurable skip background, but it is lower than at y0.**
- Strandedness: on CDS, the dUTP orientation (mate 1 antisense) is 0.914 / 0.924 / 0.929 in AM31 / AM32 / AM33, against 0.502 in the metagenome. On tRNA it is about 0.70 [strand_summary.out].
- In the y0 window, 17 of 187 RNA reads are in the minority orientation (0.091).
- Through templates are enriched in the minority orientation: 3 of 5, P = 0.0065. Joined templates: 1 of 5, P = 0.38 [tally.out].
  - So the RNA through reads look partly like gDNA.
  - Majority-orientation templates only: 4 joined / 2 through.
  - Minority-orientation templates only: 1 joined / 3 through. One joined template (SRR28822789.2556028) is in the gDNA-like orientation.
- Gap background across the index:
  - RNA: one place with a ≥ 15-nt gap that matches no event (53 nt, MAPQ 19, 2 alignments of one duplicate pair), out of 49,852 primary alignments.
  - Metagenome: 0 [bg.out].
- A blind spot in that census: the y0 window is cut off at 1135437.
  - Whole-genome realignment found a 63-nt gap (1135376-1135439) in the same gene, about 330 nt upstream. It is in frame and runs GT…AG on the gene strand.
  - It has 2 templates in 2 runs (AM31 and AM33), 0 of 34 in the metagenome, and is conserved in all 3 relatives [site63.out, site63_dedup.out, rel63_blast.tsv].
  - Its rate (2/33) is below y0's (5/10): Fisher p = 0.0043.
  - So recurring non-BHB junctions do occur in these libraries, at about a third the rate of y0's joined share.
- Partly supportive: 3 of the 8 internal soft-clip boundaries in the RNA reads in this window sit exactly at position 406, the insertion edge. The metagenome has 0 of 7 [clipends.out].
  - These are 2 templates, both in AM33, both in dUTP orientation. They read across back-junctions 692→406 (+TGCC) and 849→406 (supplementary alignment confirms) [backjunc.out, y0_SA_tags.out].
  - They could be circularised 5'-exon fragments, but 2 templates is anecdotal.

**4. Same culture as the RNA: CANNOT TELL. The metagenome is not the same sample, and the reference genome was assembled from it.**
- RNA: SAMN41101933/4/5, "log phase … AM culture incubated with methanol and yeast extract", biological replicates 1-3.
- DNA: SAMN36464503, AM_stage3, "the more simplified AM culture at the isolation stage3", yeast extract + casamino acids + methanol, 55 °C.
- The MAG GCA_030522375.1 (SAMN36373251, strain LWZ-6) is "derived-from-SRA SRR25244881, SRR25392297" [sra_meta/*.txt].
- So the 0/45 DNA result only says the stage-3 population, the one the genome came from, carries the insertion. It does not test the DNA of the RNA culture.
- No AM metagenome from the RNA time point exists in PRJNA992765. The only other AM DNA run is the 2.2 GB nanopore run, which was not downloaded [sra_meta/PRJNA992765.runinfo.csv].
- Inside the RNA libraries, the gDNA-like templates are 1 joined / 3 through. That is too few to exclude a deletion subpopulation.

**5. "Replicates" wording: CANNOT TELL / weak.**
- The positive signal is 5 templates from 2 of 3 runs. These are biological replicates of one culture, sequenced together.
- In favour:
  - The c0 study is independent (Montana State University vs Chengdu) [sra_meta].
  - The decoys stay clean on my own recount: RNA y0m 0/5, y0p 0/15; metagenome 0/55, 0/56 [decoys.out].
- Against: the A5 "removed" word depends on a single template, and the DNA condition rests on a different sample.

**Verdict: the joined junction itself SURVIVES.**
- It is exact, unique, not present contiguously in any census genome, and above the local non-BHB background.
- **"Replicates (removed, DNA clean)" is NOT ESTABLISHED.**
  - The share sits on the 0.50 bar.
  - The DNA check is a different sample, and the genome was assembled from it.
  - One joined template is in the gDNA-like orientation.
  - A recurring 63-nt non-BHB skip in the same gene shows that the background is not zero.
- Call it "spliced in 2/3 RNA replicates, DNA of the RNA culture untested".
