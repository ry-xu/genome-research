# Independent check: P. freudenreichii 23S SRL tetraloop GCAA

Written by a fresh subagent (20 Nov, Zosma's kill check). The subagent could not create this file, so Zosma saved its
returned report here unchanged in substance. Scripts and outputs are in this folder.

**Verdict: NARROWED.** GCAA is real and fixed in every P. freudenreichii genome, and the sister species read GAGA. But
the same natural sequence is already a row in Blanchet et al. 2014 (Table 1 #15, organism unnamed), and the same double
change was built and studied in E. coli by Kudrin et al. 2018. What is new is naming the carrier and showing the change
is fixed across the species.

## Method (no 15-mer regex)
- `getref.py`: E. coli rrnH 23S (NC_000913.3:225759..228662) and 16S. 2652-2667 reads CTAGTACGAGAGGACC.
- `srl_scan.py`: blastn -task blastn of E. coli 23S/16S against each assembly (e ≤ 1e-10, dust off); loci merged within
  4 kb; bases at E. coli 2659-2662 read from the alignment, then re-read from the contig on the right strand.
- `fragcheck.py`: relaxed (e ≤ 1e-5) search with a P. freudenreichii 23S for missed copies.
- `plant.py`: planted variant loops in control genomes. `readloop.py`: 27-nt anchor ending at C2658 in short reads.

## Results (tally_out.txt, strains_out.txt, meta_out.txt, sra_loop_counts.txt)
- Positive controls: E. coli K-12 GAGA 7/7; C. acnes GAGA 3/3 (16S loci match). Planted: C. acnes GAGA 2 + GCAA 1,
  E. coli GAGA 6 + TTCG 1, read back correctly.
- P. freudenreichii: 55 assemblies, 53 with a 23S (the 2 without are a 39,309-bp phage record GCA_900087275.1 and a
  24,909-bp fragment GCA_900095035.1). 120 copies, all GCAA, alignment read = genomic read for all 120; context
  TTAGTACGCAAGGACC in all. Every chromosome record has 2 copies (7 Helsinki JS assemblies carry two near-identical
  chromosome records, so 4). Relaxed search: 0 extra hits.
- Duplicates: 53 assemblies → 49 after removing identical ones (P.UF1 = MGYG-HGUT-02401, TL18 ×2, PB4 ×2, J117
  re-assembly): 112 copies, all GCAA; about 48 distinct strains. Some FAM strains look clonal.
- Type strain DSM 20271 (GCA_000940845.1, the only NCBI type-material flag): GCAA 2/2. CIRM-BIA1: GCAA 2/2.
- 12 submitters; PacBio, ONT, IonTorrent, Illumina hybrids, 454/SOLiD.
- Short reads: ERR12342576 60/61 anchor hits GCAA; SRR27779104 (TL19) 247/261 GCAA; GAGA 0 in both. Three runs labelled
  P. freudenreichii (ERR3261331, ERR3261339, ERR2221169) have 40-42% GC (species ~67%): mislabelled, excluded.
- Sisters (P. australiense, P. acidifaciens, P. cyclohexanicum, P. ruminifibrarum, including 3 type strains): 9
  readable copies, all GAGA. Thin: several are drafts; one P. ruminifibrarum assembly has no 23S.

## Prior art
- Europe PMC: no title/abstract reports natural GCAA or A2660C in Propionibacterium.
- Blanchet et al. 2014 NAR (PMC4176186), Table 1 #15 `UUAGUACGCAAGGACCG`, Freq 1, from a bacterial 23S alignment
  (E. Westhof); predicted viable, not tested; organism not named (identification is an inference).
- Kudrin et al. 2018 NAR (PMC5829649): E. coli A2660C/G2661A grows normally at 32-42 °C, is cold-sensitive, weakens
  RelA activation by vacant ribosomes. No natural carrier named.
- St-Onge et al. 2007 (PMC1865062): names P. freudenreichii, but for other sarcin-ricin motif sites, not this loop.

## Caveats
- The change is a double one: A2660C + G2661A, not "A2660C".
- 48 strain names are not 48 independent lineages.
