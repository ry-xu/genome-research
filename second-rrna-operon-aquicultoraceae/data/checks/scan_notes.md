<!-- Copied from the lab record. Paths that point to files in this repository were rewritten to repo-relative ones; any other path refers to the original working folder and that file is not included (see data/README.md). -->

# Multi-copy divergent 23S scan (src/frozen_core/div23_scan.py over shared/codon-atlas genomes, 50,607 list lines = 50,594 unique accessions, 4 lanes; all scanned except the 2 below)

## Two genomes nhmmer could not read (alphabet guess), rescued by hand
- GCF_001270025.1 (lane 3) and GCF_900199455.1 (lane 1) are clean ACGT (no IUPAC), but their first contigs are 128-137 nt
  of GC-rich sequence, so nhmmer's alphabet guess fails. The IUPAC patch does not help; `--dna` fails too, because
  RF02541.hmm is an RNA-alphabet model. `nhmmer --rna` works. Tables are in rescue/.
- GCF_900199455.1: 2 full copies (model 9-2922), equal scores 2053.6 -> nothing to flag.
- GCF_001270025.1: only fragments (1583-2922, 3-1454 on separate contigs, two 160-nt bits) -> no copy reaches
  span 2000, so div23_post.py could not flag it.
- So the scan covers the whole list with respect to flags. A rerun of div23_scan.py should pass `--rna`.

## Full-scan triage (22:30 14 Nov): no new PTC exception
- div23_post.py: 51 flagged copies (genomes with >=2 copies spanning >=2000: 4,641). All 29 partial-scan flags recur; 22 new.
- div23_check.py -> div23_class.py (data/scan/check_full.tsv, class_full.tsv): singleton 36 (18 new), ref-foreign 7 (2 new),
  foreign 4 (2 new), rrna-pile 4 (0 new), shared 0.
- div23_chim.py (data/scan/chim_full.tsv): chimera no 47, maybe 2, na 2; edge 18 of 51.
- PTC bases of the 18 new singletons (src/frozen_core/ptc_bases_pair.py vs own reference copy and E. coli rrlB; local alignments
  where the global route gave "-"; scratch bt2/):
  - canonical at all 12 checked positions, both routes: GCA_016720145.1, GCA_040766325.1, GCA_041389925.1,
    GCA_041390515.1, GCA_041392305.1, GCF_021596505.1 (Trichormus).
  - canonical where readable, foreign: GCA_963663115.2 (Paludibacter) and GCA_964536015.1 (Ignavibacteria) carry
    near-identical fragments (80% to each other), 75-84% to Alphaproteobacteria 23S.
  - DASQQM01 (GCA_035643475.1, Aquicultoria) is already a member of the Aquicultoraceae lead (cca/aq): a recall check.
  - P-loop lost, explained by sequence quality:
    - GCF_032654805.1 Kitasatospora brasiliensis (Unicycler, Illumina + MinION): operon 6 of 6 (16S too, 93%) is divergent.
      P-loop GGGG(2250-53) -> CAG with a 2-nt gap, 1-nt gaps near 2058 and 2500. 23S vs own normal copy: 33.7 mism/kb,
      26.1 gap-opens/kb; real census 23S pairs at 90-94% identity: 61.2 mism/kb, 3.3 gap-opens/kb (806 pairs). 8x the
      indel rate at half the substitutions, deletions 116 vs insertions 62: reads as an unpolished nanopore repeat copy.
      No public reads found for SAMN36026885, so not read-tested.
    - GCA_052174455.1 Arsenophonus (Illumina MAG): 80.6% to own normal, nearest census are Arsenophonus at ~81%; P-loop
      GGGG -> TGG, 2452 and 2584 changed. 19.0 gap-opens/kb vs 8.9 for census pairs at 78-84% (1,183 pairs). Reads as a
      decaying copy, not a functional PTC variant.
  - unreadable / foreign: GCA_016764715.1 Kordiimonas (78% complete; copy rearranged, nearest CALQRR01 Bacteroidia 83%),
    GCA_963717225.1 PALSA-350 (6.5% cont; 3,491 nt, best non-self 87% Streptomyces), GCA_023250515.1 Binatia (both
    copies archaeal: 85% Nitrososphaeria / Bathyarchaeia), edge-truncated GCA_028289345.1, GCA_041488375.1,
    GCA_963738755.1 (no PTC region on the divergent copy).
- Reading: the full multi-copy scan adds no PTC exception. Checked by me only.

## Two partial-scan leftovers closed (01:03 15 Nov)
Two partial-scan flags had loose ends in the 19:51 14 Nov triage (STATE). Both are indel-rich copies from nanopore-only
Flye assemblies (NCBI Datasets assembly_info: Flye 2.9, Oxford Nanopore GridION), so they read like Kitasatospora:
- GCA_046351315.1 (Polyangiaceae MAG, PRJNA1012837): odd copy 86.5% to own normal over 2,951 nt, 73.9 mism/kb and
  36.9 gap-opens/kb, against 122.2 and 9.7 for 121 census pairs (normal copy vs census, 83-90%, >2.5 kb). Local
  alignments (scratch bt2/loc.py) show 1-nt gaps near 2060-2066 and in the P-loop G run, and a gapped 2580-2588 region.
  Its "2583-2584 -> AA" depends on where the aligner puts those gaps. It is not a clean substitution, so it is not counted.
- GCA_041398565.1 (PRJNA974210): odd copy 91.4% to own normal, 76 gap-opens in 2,995 nt (25.4/kb, against 3.3/kb
  for census pairs at 90-94%). A 2-nt gap at 2061-2063 explains the earlier "-" at 2062.
- Sequencing tech for all 46 flagged accessions: data/scan/flag_tech.tsv. Nanopore-only: 16 (GridION 9, unspecified
  "Nanopore" 3, PromethION 2, MinION 2). No base rate over the 4,641 multi-copy genomes was computed, so no enrichment
  is claimed.

## Every out-of-family flag accounted for (01:17 15 Nov): data/scan/disposition_outfam.tsv
- 51 flagged copies = 46 genomes. 11 genomes are foreign / ref-foreign / rrna-pile (another organism's rRNA or a pile),
  2 are the family (DASXVE01, DASQQM01), 33 are out-of-family singletons. Nine partial-scan singletons had no PTC read on file;
  I read them now (scratch bt3/, ptc_bases_pair.py + local loc.py).
- 33 singletons: 16 canonical where read (13 at all 12 positions; 2 foreign fragments; 1 truncated before 2583),
  6 changed-looking copies all explained (3 indel-rich nanopore/hybrid assemblies; 2 decaying Arsenophonus copies;
  1 co-assembly artefact killed by long reads), 11 unread (6 truncated, 3 foreign, 2 with no dc-megablast hit to E. coli 23S).
- New: GCA_046530935.1 (Arsenophonus, Illumina MEGAHIT MAG, CheckM 100/0.54). Its third 23S (JBHEQW010000092, model
  834-2913) is 70.1% to its own normal copy, and its P-loop 2250-2251 reads GG -> AA (local alignment). The global
  ptc_bases_pair route had misplaced an A at 2251 in the NORMAL copy too; local alignment shows the normal P-loop intact.
  Decay call: copy GC 0.374 against 0.501 for the normal 23S and 0.326 for the genome (flanks 0.30-0.33, so the contig is
  host-like). The 5' third is missing. The best census hit is 70.3%. Its gap-open rate equals the family's real paralogs at
  the same identity (29.6/kb vs 29.9-34.8/kb), so the decay call rests on GC drift and the missing 5' third, not indels.
  It is 65.7% to the other Arsenophonus divergent copy (GCA_052174455) over 531 nt, so not the same copy. Two Arsenophonus
  genomes with decaying 23S copies is a pattern, not a PTC exception. Not read-tested.
- Reading: no out-of-family flag shows a PTC change that survives a check; 11 of 33 singletons could not be read.

## Fact-check caveat on the dispositions (01:36 15 Nov; FACTCHECK_div23_v2.md row 7, not redone by me)
- "canonical" means canonical at the 12 positions ptc_bases_pair.py reads. At least 5 "canonical" copies (GCA_041392305,
  GCA_016720145, GCA_046350195, GCA_041389925, GCA_047199095; also a swap in GCA_046351325), all long-read assemblies,
  89-94% to own normal with ~17-32 gap-opens/kb, show indels in or next to the P-loop G run or the 2580s, the same
  profile as the "changed-seq-error" GCA_041398565. The 16/6 split is uneven; the conclusion (no surviving PTC change)
  stands. GCA_013152945 has A2453U in both copies.
- Only 22 of the 44 out-of-family flagged genomes had any PTC reading. 4,641 is the automated run; the hand-rescued
  GCF_900199455.1 also has 2 full copies (4,642).
