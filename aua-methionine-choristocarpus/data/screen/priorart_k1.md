# Prior-art check: AUA = Met in Choristocarpus mito and Chloroparvula plastid (k1_gb)

Written 14 Nov 2026 by a prior-art subagent for Lyra. This was a read-only literature search. Downloaded texts are in the
session scratchpad (`.../scratchpad/pa/`).

## Summary

| Item | Verdict |
|---|---|
| 1. *Choristocarpus tenellus* mito NC_085326.1, AUA = Met | **NOT FOUND** for *Choristocarpus*. The source paper is still unpublished. The only direct brown-algal test (Ehara et al. 1998) reports the **standard** AUA = Ile in 5 Phaeophyceae. No *Choristocarpus* or Discosporangiales was in it. |
| 1b. AUA = Met in any Phaeophyceae / ochrophyte mito other than Xanthophyceae | **NOT FOUND**. Xanthophyceae: Ehara et al. 1997 *J Mol Evol*. Brown algae and *Chattonella* (Raphidophyceae) were explicitly found to use AUA = Ile. Thraustochytrids (Labyrinthulea, not ochrophytes) use AUA = Met (Žihala et al. 2020). |
| 2. *Chloroparvula pacifica* plastid NC_042489.1, AUA = Met | **FOUND**. It is published in the very paper that describes the genome: Turmel et al. 2019 *GBE*. *C. pacifica* is analysed by name, together with trnI(cau)/tilS loss and a modified tRNA-Met(CAU). Barcytė et al. 2025 and Lukeš et al. 2025 restate it. |
| NCBI genetic code tables | No table covers either lineage. Only tables 2, 3, 5, 13 and 21 read AUA as Met, and all are animal or yeast mito. No plastid table reads AUA as Met. NC_042489 is still annotated `transl_table=11`. NC_085326 has no transl_table, and NCBI Taxonomy gives *C. tenellus* MGCId 1 (Standard). |

## Item 1: *Choristocarpus tenellus* mitochondrion (NC_085326.1 = ON552980.1)

**GenBank record (efetch, read in full).** ON552980.1 (updated 05-MAR-2025) and NC_085326.1 (01-FEB-2024) both cite
Choi S.-W., Choi J.W., Graf L., Kawai H., Jung S.W., Choi C.G. and Yoon H.S., "Organelle phylogenomic analysis reveals
isogamous ancestral polymorphisms in brown algal evolution", JOURNAL **Unpublished**. Neither record has a note on the
genetic code. There is no `transl_table` qualifier (the default is code 1). The tRNA annotations are 2 × trnM-CAU
(product "tRNA-Met") and 1 × trnI-GAU. There is no trnI-CAU and no lysidine or tilS note. The plastid NC_085297.1 (=
ON552957) cites the same unpublished title.

**Is the paper out? No match found.** I searched Europe PMC (title phrase, "ancestral polymorphisms", Choi/Yoon and
Graf/Kawai author queries), the web, and Discosporangiales/Choristocarpaceae + mitochondrial queries. Nothing with that
title, or anything like it, has been published as of 14 Nov 2026. The related papers by the same group are:
- **Choi S.-W., Graf L., Choi J.W., Jo J., Boo G.H., Kawai H., Choi C.G., Xiao S., Knoll A.H., Andersen R.A., Yoon H.S.
  2024. Ordovician origin and subsequent diversification of the brown algae. *Curr Biol* 34:740–754.
  doi:10.1016/j.cub.2023.12.069.** I read the full text of the published-version PDF from the VTechWorks repository. It
  is plastid-only: "Brown algal phylogeny and evolutionary timeline were established with plastid genomes". Its accessions
  are ON527003, ON552953–ON552976 and OR757060–OR757061. These include the *Choristocarpus* plastid ON552957 but not the
  mito ON552980. It has 0 hits for "mitochondri", "genetic code", AUA, isoleucine or methionine. So this is not the
  paper behind the mito record.
- **Organelle genome evolution in the early-diverging brown algal lineage, Ishigeales. *J Phycol* 2026,
  doi:10.1111/jpy.70180 (PMC13280774).** I read the full-text XML (searched, with key sections read). It uses
  NC_085326 as a comparison genome; Table 1 gives "Discosporangiales Choristocarpus tenellus Mitochondria 32,351 30.60% 0
  35 0 22 3 nad6b (loss) NC_085326". It also says: "given the strong structural conservation observed in mitochondrial
  genomes of early‐diverging brown algae (i.e., Choristocarpus tenellus)". It has 0 hits for AUA, genetic code,
  isoleucine, methionine, lysidine, tilS or trnI. It does not discuss the code.

**Other brown-algal organelle papers searched, full text, 0 hits for AUA / AUA-Met** (the code is either not
discussed or stated as standard):
- Starko et al. 2021, *GBE* 13:evab124 (PMC8290108), "Genomic rearrangements and sequence evolution across brown algal
  organelles".
- *Sargassum hemiphyllum* var. *chinense* organellar genomes, *IJMS* 2024 (PMC11354929). It says "The mtDNA of 14
  Sargassum species encoded 35 PCGs with the standard genetic code".
- *Coccophora langsdorfii*, *PLoS One* 2017 (PMC5695614).
- Phaeoexplorer: Denoeud et al. 2024 *Cell*, doi:10.1016/j.cell.2024.10.049. I read the bioRxiv v1 full-text HTML
  (doi:10.1101/2024.02.19.579948) and searched it; the HAL PDF downloads were truncated. It assembled 34 organelle
  genome sets (Table S2, not opened). There are 0 hits for AUA, codon reassignment or translation table. "Genetic code"
  appears once, in a generic ENC (effective number of codons) definition. *Choristocarpus* is mentioned once, as a
  nuclear-genome MDS outlier.

**The direct test of the brown-algal mito code. Read in full (4-page PDF hosted by Kobe Univ., Kawai lab):**
**Ehara M., Watanabe K.I., Kawai H., Inagaki Y., Hayashi-Ishimaru Y., Ohama T. 1998. Distribution of the mitochondrial
deviant genetic code AUA for methionine in heterokont algae. *J Phycol* 34:1005–1008.
doi:10.1046/j.1529-8817.1998.341005.x.**
- The study used COXI sequences (~1.1 kb). The taxa were *Chorda filum*, *Colpomenia bullosa*, *Ectocarpus* sp.,
  *Pseudochorda nagaii*, *Undaria pinnatifida* (new), plus *Pylaiella littoralis*, and *Chattonella antiqua*
  (Raphidophyceae). There is no Discosporangiales or *Choristocarpus*.
- Quote: "Comparative analyses of these sequences revealed that all the algae described above bear the universal genetic
  code, including the assignment for the AUA codon."
- Quote: "no sign of modified genetic codes was detected in the analyzed COXI genes of the five phaeophytes and one
  raphidophyte".
- Quote: "the most parsimonious scenario is that this genetic code change occurred in an ancestor of the Xanthophyceae
  after this group branched off from the Phaeophyceae."
- Note: H. Kawai is a co-author both here and on the *Choristocarpus* record.

**Xanthophyceae AUA = Met (known case).** Ehara M., Hayashi-Ishimaru Y., Inagaki Y., Ohama T. 1997. Use of a deviant
mitochondrial genetic code in yellow-green algae as a landmark for segregating members within the phylum. *J Mol Evol*
45:119–124. doi:10.1007/pl00006210. I read the abstract only. Quote: "We analyzed the mitochondrial cytochrome oxidase
subunit I (COXI) gene sequences of eight algae classified as xanthophyceans and found that six manifested the expected
deviant genetic code where AUA codes for methionine (AUA/Met)... all the examined eustigmatophytes utilized AUA for Ile."
The wording "expected" hints at an earlier mention, but I found none. Ehara et al. 1998 cite Ehara 1997 as the report.

**Other ochrophyte mito papers searched, full text, no AUA = Met:**
- Eustigmatophytes. Richtář et al. 2026 *GBE* doi:10.1093/gbe/evag169 (PMC13449250): "no indications for codon
  reassignments were detected in any of the eustig organellar genomes analyzed". Also Ševčíková et al. 2016 *GBE*
  (PMC4824035).
- *Plocamiomonas* (Pelagophyceae) and Dictyochophyceae, *Sci Rep* 2026 (PMC13582874). Reports code 4 (UGA) only.
- *Leucomyxa plasmidifera*, *Open Biol* 2024 (PMC11528492). Only a tRNA-annotation mention of the lysidinylated
  Ile-tRNA(cau).
- Ehara et al. 2000 *Curr Genet* 37:29–33 (diatoms; abstract only): UGA = Trp in Thalassiosirales and "No evidence of a
  deviant genetic code was detected in the other analyzed diatom species."
- The Labyrinthulea are outside Ochrophyta, but for completeness: Žihala et al. 2020 *Mol Phylogenet Evol* 152:106908,
  doi:10.1016/j.ympev.2020.106908, report AUA = Met (from lab notes; I did not re-read it).

**Verdict, item 1: NOT FOUND.** AUA = Met has not been reported for *Choristocarpus* or any Phaeophyceae. The one
published direct test found the opposite in 5 derived brown algae (Laminariales/Ectocarpales/Scytosiphonaceae), which
supports the 1998 inference that AUA = Met is a Xanthophyceae-only change. A *Choristocarpus* AUA = Met would contradict
that scenario for an early-diverging brown alga. It would need its own check: one genome, and NC_085297 (its own
plastid) reads AUA as Ile.
Caveat: the Choi et al. organelle paper is unpublished. If it appears, it could state this.

## Item 2: *Chloroparvula pacifica* plastid (NC_042489.1 = MK085995)

**Source. I read the full text (Europe PMC XML) and searched it, reading all AUA/trnI/tilS passages.** Turmel M.,
Lopes Dos Santos A., Otis C., Sergerie R., Lemieux C. 2019. Tracing the Evolution of the Plastome and Mitogenome in the
Chloropicophyceae Uncovered Convergent tRNA Gene Losses and **a Variant Plastid Genetic Code**. *Genome Biol Evol*
11(4):1275–1292. doi:10.1093/gbe/evz074. PMID 30937436, PMC6486808. The RefSeq record cites this paper, and its
comment says it is identical to MK085995 (strain RCC4656).

Quotes:
- Abstract: "Unexpectedly, despite the disappearance of trnI(cau) from the plastome in the Chloroparvula lineage, AUA
  codons (the codons recognized by this gene product) were detected in certain plastid genes... we discovered that the
  AUA codon was reassigned from isoleucine to methionine in Chloroparvula. This noncanonical genetic code has not
  previously been uncovered in plastids."
- trnI/tilS: "The fact that both trnI(cau) and tilS are either present or missing is not a coincidence, as tilS
  encodes the tRNA(Ile)-lysidine synthase..." and "Given that trnI(cau) occurs in all Chloropicon taxa but is absent in
  all Chloroparvula taxa, one would predict the AUA codon to be less frequent..."
- *C. pacifica* by name: "all five Chloroparvula species, with the exception of Chloroparvula pacifica, exhibit much
  fewer AUA codons compared to the Chloropicon species (fig. 4C and supplementary table S1)". Also: "the protein-coding
  genes of Chloroparvula japonica (RCC2339) and Chloroparvula pacifica (RCC4656) carrying AUA codons were aligned with
  their orthologs".
- Result: "the AUA codons occupying universally conserved positions at the protein level were found to invariably
  correspond to AUG codons in Chloropicon species as well as in the other green algae examined (fig. 6)... in
  Chloroparvula plastomes, the AUA codon is decoded as methionine and that two codons (AUU and AUC) instead of three, are
  read as isoleucine."
- Fig. 6 legend: "The 21 sites containing AUA codons are distributed among eight protein-coding genes; the number below
  the gene name indicates the position of the AUA codon in the Chloroparvula pacifica gene."
- Mechanism: specific mutations in the Chloroparvula elongator tRNA-Met(CAU) (Fig. 7, "variant nucleotides observed in
  Chloroparvula pacifica"). The paper frames this as codon capture ("also known as the codon-capture hypothesis"). It
  rules out RNA editing as very unlikely.
- Careful: a separate sentence says trnI(cau) "was retained solely in... Chloroparvula pacifica". That sentence is
  about the **mitogenome**, not the plastome.

**Later restatements:**
- **Barcytė D., Žihala D., Eliáš M. 2025. Expanded diversity of pedinophytes provides a window into the evolution of
  the genetic code in organelles. *PLoS Genet* 21:e1011901. doi:10.1371/journal.pgen.1011901, PMC12574857.** I
  searched the full-text XML and read the abstract and key passages. It mentions Chloropicophyceae once: "there is only
  a single other known instance of the AUA reassignment in plastids, namely in the Chloroparvula lineage of the green
  algal class Chloropicophyceae [45]" (ref 45 = Turmel et al. 2019). The codes it reports:
  - Pedinophyte mitochondria: UGA = Trp "evolved multiple times"; AGA/AGG Arg → Ala "apomorphic for the whole order
    Marsupiomonadales"; UUA and UUG → stop in Marsupiomonadales (mtRF1a mutations); "All pedinophyte mitochondria seem
    to decode AUA as methionine".
  - Pedinophyte plastids: "an analogous reassignment seems to be evolving also in plastids of two separate pedinophyte
    lineages" (trnI(cau) and tilS absent from *Pedinomonas minor* and three Resultomonadaceae). It calls this "A
    nascent Ile-to-Met AUA reassignment", most advanced in *Oistococcus okinawensis*.
  - peDinoflagellate plastids: "apart from the previously reported Ile-to-Met AUA reassignment" (refs: Matsumoto et al.
    2011 *Mol Phylogenet Evol* 60:68–72, doi:10.1016/j.ympev.2011.04.010; Matsuo et al. 2022 *Front Plant Sci*
    13:918543), plus AGA/AGG → most likely Ala, and UUA/UCA → stop via modified pRF2.
  - 0 hits for Phaeophyceae, brown alga or Xanthophyceae.
- **Lukeš J., Eliáš M., ... 2025. Natural and artificial variations of the standard genetic code. *Curr Biol*.
  doi:10.1016/j.cub.2025.09.071.** I searched the submitted version (CWI repository PDF). Quote: "Later, however,
  plastids of certain algae joined this class, as they exhibit the Ile-to-Met reassignment of AUA [67,68]". Ref 67 is
  Matsumoto 2011 (*Lepidodinium*) and ref 68 is Turmel 2019. The review does not mention brown algae or *Choristocarpus*.
- Turmel 2019 has 21 citing papers in Europe PMC. Only the three above discuss the plastid code, judging by titles plus
  the texts I read.

**Verdict, item 2: FOUND.** The k1_gb result for *C. pacifica* is a blind recovery of Turmel et al. 2019. That paper
names this exact strain (RCC4656) and reports the same trnI(cau)/tilS absence. The lab's earlier README line ("Published
for ... *Chloroparvula*") was right; NOTE.md's "not in lab notes" was not. The only open point is bookkeeping:
NC_042489.1 still carries `transl_table=11`, and NCBI has no plastid table with AUA = Met.

## NCBI genetic code tables (wprintgc.cgi?mode=c, fetched and parsed 14 Nov 2026)

- Tables reading ATA as Met: **2** (vertebrate mito), **3** (yeast mito), **5** (invertebrate mito), **13** (ascidian
  mito) and **21** (trematode mito).
- Every other table (1, 4, 6, 9–12, 14–16, 22–33) reads ATA as Ile, including 11 (plastid), 23 (Thraustochytrium mito)
  and 32 (Balanophoraceae plastid).
- The page has no text matching Phaeo-, brown, Xantho-, yellow, Pedino-, Lepidodinium, dinoflag- or Stramenop-. The
  only "Chloro" hit is table 16, Chlorophyceae mito (TAG = Leu).
- So no NCBI table covers Xanthophyceae, Pedinophyceae, peDinoflagellate or Chloroparvula AUA = Met. These published
  cases are simply untabled.

## Sources: what I read and how

- **Full text, searched and key passages read:** Turmel et al. 2019 (PMC XML); Barcytė et al. 2025 (PMC XML); the
  Ishigeales 2026 *J Phycol* paper (PMC XML); Choi et al. 2024 *Curr Biol* (repository PDF of the published version);
  Lukeš et al. 2025 (repository submitted PDF); the eustigmatophyte (2016, 2026), *Plocamiomonas*, *Leucomyxa*,
  *Sargassum*, *Coccophora* and Starko 2021 papers (PMC XML); Phaeoexplorer (bioRxiv v1 HTML; supplementary
  Table S2 not opened).
- **Full text, read entirely:** Ehara et al. 1998 *J Phycol* (4 pages).
- **Abstract only:** Ehara et al. 1997 *J Mol Evol*; Ehara et al. 2000 *Curr Genet*.
- **Records read in full:** GenBank/RefSeq headers and features for NC_085326.1, ON552980.1, NC_085297.1 and NC_042489.1;
  NCBI Taxonomy for *C. tenellus*.
- **Not read:** Žihala et al. 2020 (cited from lab notes); Matsumoto et al. 2011 and Matsuo et al. 2022 (cited via
  Barcytė); supplementary tables of Turmel 2019 (per-gene AUA counts).
- **Could not open:** the Cell Press HTML for Choi et al. 2024 (HTTP 403). The repository copy was used instead.
