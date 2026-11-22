# Bilabrum cf. latius GCA_965178125.1: TAG -> Gln, evidence log

Written 2026-11-11 21:56 by Myrtle. Genome: data/dl_bil/GCA_965178125.1_puBilCfla1.1_genomic.fna (48.57 Mb; 43 chromosomes OZ237735-OZ237777, mito OZ237778, 67 unplaced CBCXTV contigs, 2.39 Mb).

## Screen numbers (results/stage2/GCA_965178125.1_TAG.txt, frozen rules of 21:17)
- TAG: n=190 inner columns in 121 proteins, rate_al 1.389, top meaning Q (pm 46), rel 0.56, powered. TAA n=0, TGA n=0.
- Stage 2 vs control CAA: dispersion 123/145 obs/exp (CAA 202/227.7); GC3 with/without 0.79/0.84 (CAA 0.79/0.85); usage distance 0.054 (CAA 0.080); pident with/without 47/50.7 (CAA 46/53.2); Q share at conserved columns 0.468 of 62 (CAA 0.548; own-codon identity 0.595).
- Gene ends (Clover's ends.tsv, send == slen): TGA 152, TAA 16, TAG 0 (rows); after removing duplicate loci (review, 22:21): 149 genes, TGA 136, TAA 13, TAG 0.
- Relatives (same pipeline): Cafeteria roenbergensis GCA_008330645.1 TAG n=1; Caecitellus GCA_965644315.1 n=2; Halocafeteria GCA_057816405.1 / GECG01 n=0; uncultured Cafeteriaceae GCA_964026335.1 n=0. Their gene ends, duplicate loci removed (review/terminal_codons_relatives.tsv): Cafeteria TAG 11 of 145 (7.6%), Caecitellus 35 of 172 (20.3%), Halocafeteria 36 of 149 (24.2%), uncultured Cafeteriaceae 39 of 160 (24.4%). Bilabrum 0 of 149 vs Cafeteria 11 of 145: one-sided Fisher p = 3.4e-4. The contrast is real but its strength comes from Cafeteria's low TAG use, so it is supporting, not decisive.

## Organism identity (blastn of NCBI 18S, results/bilabrum/ssu_queries.fa)
- Chromosome 12 (OZ237746.1) carries two rDNA copies matching Bilabrum latius 18S: PX596818.1 99.8% over 1535 nt, MN315515.1 99.6% over 1415 nt, OR016443.1 99.8%.
- The host dinoflagellate's 18S (Galeidinium rugatum = Dinothrix rugata AB195668.1, 99.4% over 1741 nt) is on one unplaced contig, CBCXTV010000052.1, not on a chromosome.
- All 123 TAG bodies in the stage 2 list lie on chromosomes (42 distinct OZ sequences); none on unplaced contigs.

## eRF1 (tblastn human P62495; results/bilabrum/erf1_aln.fa)
- One locus, OZ237769.1 (chromosome 35) ~455.5-457.4 kb, 63% identity to human over residues 35-422, introns.
- The eRF1 gene itself has an in-frame TAG at human L69, inside domain 1: RVNR[TAG]SVL (context AGGGTCAATCGA TAG TCTGTGCTAAGT, OZ237769.1:456067). The alignment has no gap there. Tetrahymena (TAA/TAG=Q) has Q in this column.
- Domain 1 motifs: E55 kept; TASNIKS -> TAVNIKS (S60V); YxCxxxF kept (YMCDSHF); GTS kept (upstream exon, RRDGTSMVSLILPPK).
- S60 prevalence: 3,189 UniProt reference-proteome eRF1s, 2,776 with the motif: S 2,714, T 22, K 13, A 12, G 7, F 3, N 2, L 2, Q 1 (Stylonychia, TAA/TAG=Q), V 0. S60V is unseen elsewhere.

## Literature (searched 21:58)
- Ledger (notes/known-nuclear-stop-codes.md section 1.2): no nuclear code variant known in stramenopiles. Swart2016 screened MMETSP (incl. Cafeteria and Bicosoecida sp. CB-2014) and found none outside ciliates. Labyrinthulea have mitochondrial UAG=Tyr (mito only).
- Web searches for Bilabrum / bicosoecid / stramenopile UAG glutamine found nothing.

## Residues at TAG columns (results/scan/GCA_965178125.1.stops.tsv.gz, all 191 TAG columns)
Q 55, K 20, E 18, R 11, A 10, L 10, T 10, D 9, G 9, N 7, P 7, S 6, H 5, M 3, F 2, I 2, V 3, C 1, W 1, Y 1, * 1.
Q leads, followed by Q's usual substitutes (K, E, R). Y (an amber-suppressor Tyr) is 1 of 191 and H is 5.

## tRNAs (tRNAscan-SE 2.0.13 -E, run 21:47-22:04: trna.out, trna.ss, trna.stats)
- A tRNA with anticodon CTA (reads UAG only), OZ237757.1:907679-907607 (chromosome 23), score 70.6, not flagged as a pseudogene.
  CGCGACATAGTGTAGTGGTTAGCACAGTGGACTCTAACTCCATCAGACGCAGGTTCGAATCCTGCTGTCGCGG
- A second CTA hit, OZ237737.1:1096901-1096970, score 25.5, flagged as a pseudogene.
- No nuclear TTA-anticodon tRNA, which fits TAA staying a stop. Sup-TCA is on the mitochondrion only (OZ237778).
- Isotype models (cmscan, TRNAinf-euk-iso): the CTA tRNA scores His 73.1, Gln 72.9, Glu 57.0, Tyr 51.6. The genome's
  tRNA-Gln(TTG) scores Gln 71.4, His 66.2. So the CTA tRNA is as Gln-like as the genome's own Gln tRNA. By sequence it
  is 68.5% identical to tRNA-Tyr(GTA) and 58.9% to tRNA-Gln and tRNA-Glu, so it is not a recent copy of the Gln
  tRNAs. The alignment data (Q 55 of 191, Y 1) say which amino acid it carries.
(added 22:06)

## Independent review (subagent, finished 22:22; review/, start with bil.analysis.txt and caf.analysis.txt)
Verdict: TAG -> Gln confirmed; no alternative explanation survives.
- Re-derived with a different aligner: tblastn of 76 conserved human proteins (standard code) against Bilabrum
  (76,835 aligned codons) and Cafeteria GCA_008330645.1 (77,777). Distinct genomic stop sites inside alignments:
  TAG 201 (166 aligned to a residue, 149 >= 10 from HSP ends, 91 clean), TGA 28 (4 aligned, 1 clean), TAA 6
  (1 aligned, 0 clean). Cafeteria: TAG 0, TGA 13, TAA 0.
- The 91 clean TAGs: Q 44, K 11, E 5, S 5, R 5, others <= 4. Best-fitting codon profile is CAG's (log-likelihood
  ratio 1.005 per site; K 0.120). 36 proteins on 30 chromosomes, incl. enolase, PCNA, BiP, Sec61A, RPB1/2/3, PRP8,
  SF3B1. 2.16 TAG per 1,000 aligned codons.
- Contamination: GC3 0.771 (HSPs with a clean TAG) vs 0.779 (without); 28 of 30 TAG chromosomes also carry TAG-free
  HSPs; none on unplaced contigs. A ciliate would also read TAA; there are no clean TAA sites.
- Point errors: at TAG's rate, 37.6 TAA and 39.6 TGA would be expected among aligned columns; 1 and 4 observed.
- Gene ends: duplicate loci removed, 0 of 149 (see screen numbers above). The reviewer's own tblastn gene-end test
  has too few loci to separate the genomes (Bilabrum 40 TGA, 6 TAA, 1 TAG; Cafeteria 43, 2, 1).
- Not excluded from DNA alone: occasional TAG termination, or RNA editing (needs RNA-seq).
- Reviewer's correction, applied above: gene-end statements now carry denominators and deduplicated counts.

## RNA: no Bilabrum cf. latius RNA-seq, but Bilabrum sp. RNA in 1KP algal transcriptomes (added 22:47)
- No reads under any Bilabrum taxid at ENA/SRA (taxids 2721177, 2604063, 3115832, 3503366, 3503367). The genome is an
  ASG cobiont of the dinotom dinoflagellate Dinothrix rugata (Galeidinium rugatum); the host has HiFi and Hi-C only.
- Vierock et al. 2022 (Sci Adv 8:eadd7729, WiChR; Supplementary Text S1; first author corrected 23:45, was misattributed to Govorunova) found Bilabrum sp. as contamination in four 1KP algal
  transcriptomes (LLXJ Chroodactylon ornatum, JIWJ Acrosiphonia, ALZF Halochlorococcum, JQFK Picochlorum); 18S 99-100%
  to B. latius HFCC35 (MN315515.1), which matches this genome's rDNA at 99.6%. They wrote that transcript
  1KP_LLXJ_2042068 (channelrhodopsin B1ChR2) "appeared to contain premature stop codons" and built the protein from a
  stop-free haplotype in ERR2041020 (Data S3, "Bilabrum sp..gbk"; scratchpad/dl_wichr_supp).
- That premature stop is a TAG at B1ChR2 codon 99 (ESTAFSLN*VWYSDY, nt TTCTCCCTGAAC TAG GTCTGGTAC). The stop-free
  haplotype has CAG (Q) there; the two differ at 32 of 1,171 nt, almost all third positions. The expressed protein and
  the PDB structures 8JXO/8JXP have Q there. The only other stop in the LLXJ transcript is the terminal TGA.
- This genome's ChR2 (OZ237746.1 ~86.8-88.8 kb, 84% aa identity to B1ChR2, so a different strain or species or a
  paralog) has CAG at codon 99 and TAG at codon 11, where B1ChR2 has Q (CVCQKEP). TAG and CAG swap in both directions.
- Control: HBQF01 (MMETSP re-assembly of another Chroodactylon culture) has no Bilabrum (no >= 97% hits, best 82%;
  tblastn 76-protein panel: TAG 0.05 per 1,000 aligned codons, none opposite Q in clean columns).
- Next: map 1KP LLXJ reads (ERR2041128) to the genome and read the codon at the 91 clean TAG sites.

## RNA result 1: 1KP LLXJ reads (ERR2041128, Chroodactylon culture with Bilabrum sp.; added 23:15)
Streamed from the SRA AWS mirror into minimap2 -x sr (rna_map/run_1kp.sh); per-site codons from rna_map/site_reads.py,
joined to the reviewer's sites by rna_map/site_table.py (ERR2041128.table.txt); strain comparison by
rna_map/codon_pairs.py (LLXJ.pairs.txt, LLXJ.codons.tsv.gz). The .sra file was deleted after mapping.
- 749,234 primary read alignments (mapq >= 10) on the chromosomes at 95.8% identity, so this is a related Bilabrum
  strain or species, not this genome's. 12,864 on unplaced contigs.
- The TAG sites are in mRNA. Of the 91 clean TAG sites, 58 have >= 2 reads: the majority read codon is TAG at 39,
  CAG at 18, CAA at 1 (reads: TAG 495, CAG 269). At clean TAG sites opposite human Q, 22 of 29 read TAG.
  All 166 aligned TAG sites: 92 covered, TAG 64, CAG 25, 3 others.
- TAG swaps with CAG between the two strains at the synonymous rate. Over the 76-protein panel (48,106 codons
  called from >= 2 reads, 5,014 differ from the genome):
  - genome CAG read as TAG at 21 of 1,677 sites (12.5 per 1,000); synonymous Leu CTG -> TTG 29 of 2,233 (13.0).
  - C->T changes that alter the protein are rare: His -> Tyr 1 of 957, Arg CGG -> Trp 0 of 547, Pro -> Ser 1 of
    1,941 (together 2 of 3,445, 0.6 per 1,000; one-sided Fisher vs CAG -> TAG p = 7e-9).
  - C->T changes that make TAA or TGA: CAA -> TAA 0 of 88, CGA -> TGA 1 of 230 (together 1 of 318; vs CAG -> TAG
    p = 0.11, too few sites to be decisive on its own).
  - genome TAG read as CAG at 25 of 93 called sites (269 per 1,000); synonymous Leu TTG -> CTG 69 of 340 (203).
  So in conserved proteins TAG <-> CAG behaves like a silent change in both directions, while protein-changing C->T
  changes are about 20 times rarer. A stop codon arising at CAG in the middle of conserved genes would be purged
  like the missense changes; it is not.
- Not excluded: the reads are a different strain, so they show TAG is in that strain's mRNA, not this genome's.

## RNA result 2: all three 1KP runs (ALZF ERR2041020, JIWJ ERR3487419 added; 23:34)
- Per run (rna_map/RUN.summary.txt, RUN.table.txt): ALZF 123,209 chromosome alignments at 95.6% identity, JIWJ
  256,742 at 95.8%. Clean TAG sites with >= 2 reads: ALZF 11 (TAG 5, CAG 5, AAG 1), JIWJ 33 (TAG 25, CAG 7, CAA 1).
- The three cultures carry one Bilabrum strain (rna_map/strain_compare.txt): read codons agree at 97.9-99.0% of
  shared panel codons, against ~89.5% agreement with this genome, and where LLXJ differs from the genome, ALZF and
  JIWJ share its codon at 88% and 89%. So the runs are pooled, not counted as replicates.
- Pooled, clean TAG sites (rna_map/all3.table.txt): 60 covered, majority read codon TAG 41, CAG 18, CAA 1
  (reads TAG 642, CAG 356). Opposite human Q: 29 covered, TAG 22. All aligned TAG sites: 99 covered, TAG 69, CAG 27.
- Pooled strain comparison (rna_map/ALL3.pairs.txt; 53,415 codons called, 5,628 differ from the genome):
  CAG -> TAG 23 of 1,875 (12.3 per 1,000); synonymous CTG -> TTG 36 of 2,506 (14.4; one-sided Fisher p = 0.32 for
  CAG -> TAG being lower). H -> Y, R -> W and P -> S together 3 of 3,945 (0.8 per 1,000; p = 4e-9 vs CAG -> TAG).
  CAA -> TAA and CGA -> TGA together 1 of 375 (2.7 per 1,000; p = 0.07). TAG -> CAG 27 of 100 (270 per 1,000);
  TTG -> CTG 85 of 395 (215).
- Per-run CAG -> TAG vs CTG -> TTG: LLXJ 12.5 vs 13.0, ALZF 5.7 vs 10.9, JIWJ 8.1 vs 15.8 (the smaller runs call
  fewer sites, mostly the same ones as LLXJ).
- Deleted after mapping (re-fetch from https://sra-pub-run-odp.s3.amazonaws.com/sra/RUN/RUN): ERR2041128.sra,
  ERR2041020.sra, ERR3487419.sra. The PAF alignments are kept.
