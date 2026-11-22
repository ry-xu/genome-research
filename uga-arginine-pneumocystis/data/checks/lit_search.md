# Literature check: UGA = Arg (Pneumocystis, Craticula) and AGA = Ser (Sympodiomycopsis)

Fresh subagent, finished 2026-11-11 14:48 EST; saved verbatim by Clover. Scripts and downloads: scratchpad lit_new/ (session-local).

## 1. Direct prior reports
**Craticula vanensis UGA = Arg: already reported.** Yılmaz E, Gastineau R, ... Witkowski A. 2026. "Craticula vanensis
sp. nov., a new species of diatom with a complex mitogenome, inhabiting soda lakes in central Europe and Anatolia."
*Sci Rep* 16 (July 2026). DOI 10.1038/s41598-026-62093-2; PMID 42463532; PMC13672107.
- Abstract: "the use of an alternative genetic code with the putative reassignment of the UGA codon to arginine".
- Results: Codetta "suggested that it could instead be reassigned to arginine (Supplementary File S2) ... In the absence
  of confirmation by protein sequencing and taking into account that GenBank does not currently provide a corresponding
  translation table for this reassignment, code 4 was retained for annotation."
- Discussion: never documented for diatoms; suspected in Trentepohliales plastids and Oikopleuridae mitochondria.
- Other Stauroneidaceae (Fistulifera solaris, F. saprophila, Proschkinia staurospeciosa, Schizostauron trachyderma)
  "are decoded using the standard genetic code".
- Not covered: RNA editing; absence of a UCA tRNA; per-position conserved-Arg analysis.
- Count mismatch: its Table 2 gives 119 UGA, 148 UGG, 341 canonical Arg codons; we count 62 TGA at aligned sites
  (its gene set may include intron ORFs: 34 proteins, 29 tRNAs). GenBank PV031684 still has an "Unpublished" title.

**Pneumocystis UGA = Arg: nothing found** (neither reassignment nor editing).
- Ma L et al. 2013 *FASEB J* 27:1962-1972 (DOI 10.1096/fj.12-224444; PMC3633818) read UGA as Trp under code 4: "there
  is no tRNA gene containing the corresponding anticodon UCA. It is possible that the C in the wobble position of the
  tRNA-Trp with anticodon CCA present in all 3 genomes is modified to recognize the UGA codon". Also: "Six codons, AAC,
  GAC, CGC, CGG, AGG, and GGC, are not used at all in P. jirovecii... The UGA codon (Trp) is found in all 3 species."
- Cissé OH et al. 2018 *mBio* 9:e00381-18 (PMC5941068): no "genetic code", UGA, TGA, codon, arginine in full text.

**Sympodiomycopsis kandeliae AGA = Ser: nothing found.** PP995850/PP995849 unpublished (submitted 2024-07-06;
Brejová, Hodorová, Nosek). Brejová B, Hodorová V et al. 2026 *Nucleic Acids Res* (DOI 10.1093/nar/gkag566;
PMC13234504) sequenced these species but is a telomere paper; full text has no AGA/serine/genetic code/codon and no
PP99585x accessions. NCBI Taxonomy: mito code 4.

## 2. Closely related prior work
- Maturová E, Žárský V, Eliáš M. 2026. "An unprecedented codon reassignment and two new genetic code variants in
  ciliates." bioRxiv DOI 10.64898/2026.11.04.764322 (~4 Nov 2026). Licnophora: "UGA exhibits homonymy, encoding
  arginine when present in-frame and terminating translation at the coding sequence ends ... the mechanism of UGA
  decoding as arginine in Licnophora and certain other spirotricheans remains elusive." First claimed nuclear UGA = Arg.
- Fang J, Liu B, Liu G, Verbruggen H, Zhu H. 2021 *Front Plant Sci* 12:780054 (PMC8692980): Trentepohliales plastid
  in-frame UAG/UGA at conserved Arg positions, "may be reassigned as arginine".
- Pichon J, Luscombe NM, Plessy C. 2020 *F1000Research* 8:2072 (PMC7047913): tunicate mito, "B. stygius uses TGA at
  positions where all other species would encode Arg ... numbers are currently too small"; 3 codons.
- Barcytė D, Žihala D, Eliáš M. 2025 *PLoS Genet* (DOI 10.1371/journal.pgen.1011901): U-to-C editing UGA -> CGA in the
  plastid of a green dinoflagellate (TGD; also Matsuo et al. 2022 *Front Plant Sci* 13:918543), the only documented
  case outside land plants; also mito AGA/AGG Arg -> Ala in Marsupiomonadales.
- Gastineau R et al. 2022 *Front Mar Sci* 9:941506 (Amicula micronesica, raphid diatom): no tRNA-Trp(UCA); TGA only in
  mttB, nad2, rps10, rps11; four hypotheses, none Arg or editing. Worth testing.
- Gastineau R et al. 2026 *Sci Rep* (DOI 10.1038/s41598-026-58736-z): Plocamiomonas has tRNA-Trp(UCA), TGA 88 vs TGG
  36, all 11 cox1 Trp encoded by TGA (contrast case: code 4 really Trp).
- Hegedusova E, Brejová B, Tomáška L, Sipiczki M, Nosek J. 2014 *Curr Genet* 60:49-59 (Jaminaea angkorensis):
  "stop-to-tryptophan reassignment of the UGA codon was accompanied by structural alterations of tRNA-Trp(CCA)";
  full text not reached (AGA not known).
- Fungal mito editing: Teichert 2020 (plastid U-to-C only); Wang C et al. 2020 *Mitochondrion* 51:88-96 abstract
  claims "new models of RNA editing" in Beauveria/Cordyceps mitochondria (weak). No report found of mitochondrial
  U-to-C mRNA editing in fungi, diatoms or other stramenopiles.

## 3. NCBI tables (gc.prt 4.6, updated 2024-09-23; 27 tables)
No table has UGA = Arg. AGA = Ser only in tables 5, 9, 14, 21, 24, 33 (all animal mito). No fungal or other
non-animal table has AGA = Ser.

## 4. Not reached
Jaminaea 2014 full text; Sesterhenn et al. 2010 full text; Brejová 2026 supplement (corrupt zip); Semantic Scholar
(429); Cissé 2021, Pungan 2022, Ma 2024 (found, not opened).
