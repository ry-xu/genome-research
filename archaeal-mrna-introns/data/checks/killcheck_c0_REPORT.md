# Kill-check c0: 22-nt removal at CP186953.1:1089656-1089678

Scripts and outputs (`*.py`, `*.out`) are in this folder.

**1. Counts: SURVIVES.** Templates with exactly a 22N gap (one possible placement, 0 nt microhomology), runs 81-89: 1,3,5,5,8,5,8,7,8 = 50. Through: 0,2,1,0,1,1,0,1,0 = 6. Collapsing identical mate coordinates gives 45 joined and 6 through, matching the claim. No template mixes a joined and a through mate. Ten more templates are clipped at the insertion edges. All 10 clips match the joined form and none match the insertion, so 45 is conservative. One through template is clipped 46S/93S, probably a chimera.

**2. Read quality: SURVIVES.** 80/81 joined alignments match the 22-nt-removed genome exactly; one has 1 mismatch. All are MAPQ 60, with 0 SA tags and 0 adapter hits. The shorter flank is ≥21 nt (median 52), and ≥25 nt in 44/50 templates. Five reads carry long unrelated clips, likely ordinary chimeras.

**3. Origin elsewhere: SURVIVES.** With minimap2 -x sr against the whole genome plus all 14 Methanomethylicales genomes, 0/81 joined reads align contiguously; each needs a 20-23 nt deletion. Junction 16-40-mers: 0 hits (genomic control 16-mer: 2). Three other Methanosuratincola genomes carry the insertion in their DNA: 23/22/20 nt, ends ACGG…CCGGCCGAA. In all four genomes, removing it fuses the frames into one ORF of 368-510 codons (06155+06160 gives 364 aa, ending at the 06160 stop). An assembly error is therefore very unlikely.

**4. Strand: SURVIVES.** In this window, mate 1 is minus-strand 123/123 and mate 2 plus-strand 129/129. Joined: m1− 45, m2+ 49. Through: m1− 6, m2+ 6. This is dUTP sense for the forward gene, so the through reads are RNA too, not gDNA. All windows: 41,184/41,184 pairs have mates on opposite strands; median mate-1 strand share 0.990 over 34 windows.

**5. Other explanations: CANNOT TELL.** The SAMs cannot rule out a culture subpopulation whose DNA has a precise deletion. The conserved insertion ends could even mark an element that excises at the DNA level. They also cannot rule out reverse-transcriptase skipping, though the junction has 0 nt microhomology and the ratio is steady across 9 runs. The decisive test, the same culture's metagenome reads, has not been done.
