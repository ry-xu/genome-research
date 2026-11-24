# c0 element across all 62 Methanomethylicia species reps (Sabik, 20 Nov, completes the 36-genome survey)

- Fetched the 21 reps missing from the codon-atlas cache (`fetch.sh`, NCBI FTP, 7.8 MB, `ext/`); Mintaka's 5 in `work/tack-rna/ext/msur5/`.
- tblastn of c0.faa (`tbn.tsv`, `tbn21.tsv`); `breaks.py` (HSP pairs + gene-strand window); `inframe.py` (translation of HSPs that run through an in-frame insert).
- Site scan `ins_all.py` (= pkgcheck/ins9.py with the conserved end loosened from CCGGCCGA[AG] to GGCCGA[AG]) on all 62 (`paths62.txt` -> `ins_all.out`).
  Every hit is 1 site per genome on the c0 ortholog contig, between ...ACGGTAT and TTAG...; 37 genomes have 0 sites.

Results (counts from ins_all.out / tblastn):
- Methanomethylicales 25/33 carry the element, all in Methanomethylicaceae, 9 genera:
  Methanosuratincola 9/9, WYZ-LMO10 5/5, WYZ-LMO11 3/3, JANXER01 2/2, VGLU01 2/2, JAULTD01 1/1, JBMFMY01 1/1, CAZQQB01 1/1, JACTUC01 1/2.
  Without it: Methanomethylicus 2 continuous orthologs (38, 48 % aa) + 4 no ortholog >= 33 %; JACTUC01 GCA_015660995.1 continuous (37 %); JAGHCH01 no ortholog.
- Lengths 20 (12), 21 (1), 22 (4), 23 (4), 24 (2), 33 (2). 20 shift the frame; the 5 in-frame ones (21, 24, 24, 33, 33) each carry a stop inside the insert.
  All end CGGCCGA[AG] (AA 19, AG 6).
- Other orders 0/29: Culexarchaeales 12/12 continuous, Nezhaarchaeales 15 (13 continuous, 1 truncated at a contig end, 1 no ortholog), CALHIW01 1 continuous + 1 no hit.

Corrections to my 12:56 post / package README lines 44-48:
- JAULTD01 (GCA_038801785.1) does carry a 33-nt in-frame insert with a stop (GCGGCCGCCACTGACGGATGGCTGCCGGCCGAG); tblastn aligned straight
  through the stop, so my HSP classifier called it one piece. I had said "no insert". JAGHCH01 has no ortholog >= 33 %, not "no insert".
- The 36-genome counts (3 more genera; 9 Culex; 10 Nezha) become 8 more genera; 12; 13 + 1 truncated.
