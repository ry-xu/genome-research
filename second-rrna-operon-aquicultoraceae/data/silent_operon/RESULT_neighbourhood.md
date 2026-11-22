<!-- Copied from the lab record (work/silent-operon/, 22 Nov 2026). Paths such as out/, src/, recheck/ or frozen-core/ refer to the original working folder; the small tables from out/ are in this folder (see ../README.md), the rest is not included. -->
# T2 Neighbourhood of the divergent operon (B_origin, 22 Nov)

**Finding.** The divergent operon sits in an ordinary chromosomal neighbourhood with no mobile, phage, defence,
toxin-antitoxin, RNase, rRNA-modifier or ribosomal-protein gene. 16 gene families sit next to it in 3 or more of the
6 genera and next to none of the 22 normal operons. The only notable one is a ~60-aa, two-helix membrane peptide
1.0-1.9 kb upstream, on the same strand, in all 6 genera.

Numbers (+-10 kb, 17 divergent units / 243 genes, 22 normal units / 578 genes in the family; Swiss-Prot, Pfam, phage proteins):
- **Suggestive classes:** divergent flanks 0 genes in any class; normal flanks ~40 (PIN + Phd/YefM toxin-antitoxin 417-820 nt
  downstream of the Aquicultor_A normal operon; Cas1/Cas2/Cas4 near CADDYS01 normal; LAGLIDADG homing endonucleases in
  JASEIR01 normal 23S; RuvC; a phage integrase and a transposase each once). E. coli/Vibrio/T. bispora/Tenuifilum flanks:
  ribosomal proteins, rRNA methyltransferases, EF-G etc. as expected.
- **Conserved with the divergent operon (>= 3 genera, 0 normal loci):** 176-aa unannotated protein ~200-330 nt downstream,
  same strand (6 genera); DUF5714; ferritin/rubrerythrin; UPF0104 membrane protein; ADP-forming acetyl-CoA ligase alpha;
  UPF0280; polysaccharide deacetylase (once at a normal locus); DUF5658; ARC6-like; MscS channel; GatC/GatA; DNA ligase LigA;
  Rce1-like protease; HAMP-SpoIIE phosphatase (4 genera). The 176-aa protein also occurs in non-carriers next to a contig
  break where the operon would be (Fomalhaut T3), so it marks the locus, not the operon.
- **The upstream peptide:** 54-61 aa (94 in DASXVE01 with an N-terminal extension), motif FWNILY...HWAGLI, 2 TM helices,
  50-80% identical between genera. It belongs to a paralog family with 2-4 other copies per genome away from the operon
  (tblastn), and has no confident Pfam hit. Its size and shape are like a type I toxin or a small membrane protein; nothing more can be said from sequence.
- **No island signal:** flank GC within 2 SD of the genome's 20-kb windows in every genome with enough windows (z -1.8 to +0.3);
  no tRNA gene at either operon or within 10 kb (one tRNA-Glu in one CADDYS01 genome); the normal operons also lack ITS tRNAs.
- The divergent rRNA genes are GC-richer than both the normal copies and the genome (23S 0.52-0.58 vs normal 0.49-0.53,
  genome 0.43-0.50), i.e. not drifting toward genome GC as a decaying copy does (see T3).

**Mechanism.** Nothing suggests a mobile element, phage, defence island or a dedicated rRNA-processing/modification
partner encoded beside it. The neighbourhood is a fixed, vertically inherited locus. If the operon has a partner gene,
the only candidate visible is the small upstream membrane peptide, which is speculative.

Files: `src/flanks.py`, `src/t2_flanks.py`, `src/t2_island.py`; tables `out/t2_orthogroups.tsv`, `out/t2_unit_genes.tsv`;
peptide `out/smallorf.faa`, `out/smallorf.aln`.
