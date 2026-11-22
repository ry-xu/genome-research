# Fact-check of outbox/2026-11-14-choristocarpus-aua-met/post.md (fresh subagent, finished ~10:45 14 Nov)
Saved by Lyra (the subagent cannot write here); condensed. Its scripts: scratchpad/factcheck_chor/code/recount.py, rna.py; out/recount.out.

## Recounted with its own code (own tblastn, own NCBI downloads)
- WGS contig ATA 88/6, ATG 27/2, ATT 3/245, ATC 2/13; NC_085326 ATA 88/5, ATT 3/242; Discosporangium 3/51; Saccharina 1/116;
  Phaeostrophion 2/69. 20 of 21 genes. 6/264 = 0.0227, expected 2.11, P = 1.10e-137 (NC record 5/260, P = 4.6e-144).
  Saccharina ATG at 88/88.
- AB899175 (cox1 partial, KU-1152, Hanyuda/Kawai 2014, Sanger) is 100% identical over 608 nt to both genome sequences and
  holds 4 ATA@Met sites incl. the cox1 PRMNN motif (AGA ATA AAT AAT, NC_085326 pos 17541). Not an assembly artefact.
- RNA: ~1.1M read-1s of ERR12954301: reads starting inside genes 116 sense / 10 antisense -> stranded RNA library, DNA at
  most ~16%. 39 ATA@M sites: A 77, G 0, T 1. (Only 5 of 35 genes on minus strand; not split by strand.)
- Strict rule (>=5 refs, unanimous): Choristocarpus ATA 68 Met / 0 Ile; Discosporangium 0/28, Saccharina 0/42, Phaeostrophion 1/35.

## Required fixes
- WRONG: Xanthophyceae is not the sister class of Phaeophyceae; Schizocladiophyceae is (RESULT.md, Hwang 2026). priorart.md has the same error.
- INCONSISTENT: "second genome gives 88/6" -> the body already uses the Genoscope contig; NC_085326 gives 88/5.
- Strain: BioSample SAMEA115803497 says the Genoscope culture was first labelled Discosporangium KU-2346 (Hydra, Greece, 2004),
  "probably identical to KU-1152"; GenBank LC664917-47 (Akita 2022) label KU-2346 as C. tenellus. ON552980 has no strain
  qualifier. -> "probably the same culture".
- Overstated: "animals, yeasts" (echinoderms, flatworms, cnidarians read AUA Ile); Sphaerosorus a "known case" (Ehara 1997
  is class-level; Collins 2026 PV826132 gives no code) -> "expected from its class"; "20 of 21 protein-coding genes" (35
  genes; 21 = those with a conserved Met coded ATA/ATG); "five brown algae from other orders" (Phaeostrophion and
  Protohalopteris both Sphacelariales -> four orders); "does not edit" -> "no sign of editing"; "modified to read AUA" ->
  presumably; contaminant identity weak (Schizocladia also 73.2%) -> cite AB899175; P-value null ignores the easy AUG->AUA
  transition -> lead with Ile codons filling 67% of Met columns and ATA ~absent from Ile columns; "about 20,000 genomes" ->
  "records"; 13.7M reads ~4% of the 181.6M-pair run; Mediterranean -> "this strain is from".
- Headline: fair at the usual comparative bar (as Ehara 1997); add "by codon usage" / "has reassigned".
- "No brown alga was known to do this": safe as worded for "known"; "not previously reported, as far as we could find" safer.

## Could not check
Proteome; tRNA charging/modification (no tRNAscan/aragorn); the full RNA run; which culture Genoscope sequenced;
Denoeud 2024 Table S2; Ehara 1997 full text; non-English literature.

## Sources it opened
Ehara 1998 PDF (Kawai Kobe page), Hwang 2026 PMC13280774, Collins 2026 PMC13485189, Cornwall Red Data Book, Europe PMC records
(Ehara 1997, Žihala 2020, Barcytė 2025, Kawai 2007), ENA filereport ERR12954301, ENA XML SAMEA115803497, NCBI taxonomy
116064/116065/376270, nuccore NC_085326, ON552980, PV826132, AB899175, LC664940.
