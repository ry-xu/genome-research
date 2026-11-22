<!-- Copied from the lab record. Paths that point to files in this repository were rewritten to repo-relative ones; any other path refers to the original working folder and that file is not included (see data/README.md). -->

# GCA_025542975 divergent 23S: an assembly artefact (killed)

Candidate from the partial 23S scan: Aristaeellaceae UBA11524 sp000437595 MAG (gut, colon-model co-assembly of
4 metagenomes, PRJNA722408, BioSample SAMN18871268). One of 4 23S copies (JAHQGG010000009.1:296413-299287) is
94.2% to the other three, evenly along its length, with a 1-nt change in the P-loop G run (E. coli GACTGGGGCGG vs
GATTAGGCGG; the aligner calls it G2251A). No other scanned copy has the motif; the nearest census 23S is 0.942.

Against it, from the genomes:
- The same locus in two conspecific genomes (GCA_025993545, GCA_037069585; flanks syntenic, 93-97%) holds a
  normal 23S (99.4-99.6% to this MAG's normal copy).
- The divergent copy's own 16S, in the same operon, is 100% identical to the MAG's other three 16S copies.

Long reads (nanopore batch-1 runs of all four source samples: SRR14780132/31/30/29; *.paf here; reference
ref2.fa = the c9 locus as assembled, plus the same locus with the MAG's normal 23S swapped in, plus two normal
operons; script ont_swap.py):
- Reads anchored in the c9 locus's unique flanks and covering >=1 kb of the 23S: 7. All 7 carry the normal 23S.
  ont_swap.py called one of them (SRR14780131.162216) "div" from minimap2 AS, but its 23S portion blasts 89.3% to
  the normal 23S vs 86.4% to the divergent one; the AS gap came from alignment extent, not the 23S bases.
- Reads covering >=1 kb of a 23S at >=0.80 identity: 39 best on a normal copy, 1 "best" on the divergent copy
  (the same read 162216). If the divergent copy were 1 of 4 real copies, about 10 would be expected.
- So the divergent 23S is not in this organism's genome in these samples. Most likely a related taxon's 23S joined
  to this host's 16S/ITS by co-assembly. Not pursued further.

Lesson for the full triage: a lone divergent copy in a MAG, whose operon 16S is identical to the other copies,
should be checked against conspecific genomes at the syntenic locus before anything else.
