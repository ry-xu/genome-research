# Prior art: mitochondrial AUA = Met in *Choristocarpus tenellus*

Checked 2026-11-14 by a literature subagent. Only public sources were used.

## Verdict: NOT FOUND for *Choristocarpus* or *Discosporangium*. The change itself is known elsewhere in stramenopiles.

- No source reports AUA = Met, or any other non-standard mitochondrial code, in *Choristocarpus* or *Discosporangium*.
  - Every paper that uses the *Choristocarpus* mitogenome (NC_085326 / ON552980) says nothing about a genetic code or
    uses the standard code.
  - NCBI assigns mito code 1 (Standard).
  - The GenBank records carry no transl_table, so they are translated with table 1.
- The same reassignment, AUA Ile → Met, is published in:
  - **Xanthophyceae**, the class closest to Phaeophyceae (Ehara et al. 1997);
  - **Labyrinthulea** and **MAST8b**, both stramenopiles (Žihala et al. 2020).
- Ehara et al. 1998 (*J. Phycol.* 34:1005; co-author H. Kawai) sampled only crown brown algae. These were 5 new cox1
  sequences plus *Pylaiella*, from Ectocarpales, Laminariales and Chordaceae. They concluded that the change "most
  probably occurred in an ancestor of the Xanthophyceae after it had branched off from the Phaeophyceae" and called an
  AUA = Met state shared with Phaeophyceae "highly unlikely". No Discosporangiales were sampled.
- So a *Choristocarpus* AUA = Met would be new for brown algae. If *Discosporangium* reads AUA as Ile, it would be a
  separate origin within Discosporangiales. It is not covered by the Xanthophyceae finding.

## Q1. *Choristocarpus* / *Discosporangium*, the Choi et al. records and Phaeoexplorer

- **The Choi et al. mito paper is unpublished.**
  - Both GenBank records (efetch, read today) still give the reference as "Organelle phylogenomic analysis reveals
    isogamous ancestral polymorphisms in brown algal evolution", JOURNAL "Unpublished".
    - NC_085326.1 was created 01-FEB-2024.
    - ON552980.1 was updated 05-MAR-2025 and submitted 19-MAY-2022.
  - The title was not found in any of these:
    - Crossref (bibliographic query);
    - Europe PMC (including preprints; queries "Choristocarpus", "isogam* AND brown alga* AND organelle 2022-2026");
    - WebSearch (exact phrase "isogamous ancestral polymorphisms").
  - The records have no /transl_table qualifier on any of their 35 CDS, so their translations use table 1 (AUA = Ile).
- **Choi et al. 2024, *Curr. Biol.* 34:740–754, "Ordovician origin and subsequent diversification of the brown algae"**
  (doi:10.1016/j.cub.2023.12.069; same group; *C. tenellus* is in the tree).
  - FULL TEXT (NSF PAR accepted manuscript), searched.
  - It is plastid-only: "The plastid genome data that are produced in this study were publicly deposited…"
  - There is no mitochondrial analysis and no genetic code text.
- **Hwang, Kim, Choi, Kawai & Yoon 2026, *J. Phycol.* 62(15), "Organelle genome evolution in the early-diverging brown
  algal lineage, Ishigeales"** (doi:10.1111/jpy.70180, PMC13280774, epub 26 May 2026).
  - FULL TEXT (Europe PMC XML), searched.
  - Its table includes "Discosporangiales *Choristocarpus tenellus* Mitochondria 32,351 30.60% … NC_085326". The text
    calls the *Choristocarpus* mitogenome an example of "strong structural conservation".
  - Annotation: "initially performed using GeSeq server … and subsequently verified with BLASTx".
  - No hits for "genetic code", "AUA", "isoleucine", "methionine" or "translation table".
- **Phaeoexplorer: Denoeud et al. 2024, "Evolutionary genomics of the emergence of brown algae as key components of
  coastal ecosystems"** (*Cell* 187:6943, doi:10.1016/j.cell.2024.10.049; bioRxiv 10.1101/2024.02.19.579948).
  - The bioRxiv v1 FULL TEXT was searched. The Cell version returned 403 and was not read.
  - It reports 34 sets of plastid and mitochondrial genomes. Methods: NOVOPlasty, "Annotation of protein-coding genes was
    performed with GeSeq v2.03", ARAGORN for tRNAs.
  - The only mitochondrial result is tree discordance in Ectocarpales.
  - No genetic code, AUA or translation table text.
  - Its Table S2 (organelle genome characteristics) was not seen.
- **Other brown algal mitogenome papers, 2020–2026 (FULL TEXT searched for code terms):**
  - None mentions *Choristocarpus* or a code other than the standard one.
  - Two state the standard code outright:
    - *Dictyota coriacea* 2025 (PMC12565469): "MFannot …, using Translation Table 1 (standard)";
    - *Sargassum hemiphyllum* 2024 (PMC11354929): "35 PCGs with the standard genetic code".
  - The others have no code statement:
    - *Cystophora brownii* 2026 (PMC13485213);
    - *Padina usoehtunii* 2024 (PMC11290263);
    - Starko et al. 2021 *GBE* (PMC8290108);
    - Choi et al. 2020 *Sci. Rep.* (PMC7005149).
- **GenBank has no complete *Discosporangium* mitogenome**, only cox1 partials (AB899176.1, LM995414.1).
  - Your *Discosporangium* comparison must therefore come from a WGS assembly. The Phaeoexplorer one is GCA_964200235.2,
    strain MT17-79, Malta (RCC7272).

### Strain note (not literature, but relevant)

- The BioSample of GCA_964200185.2 (SAMEA115803497) is a "virtual sample" that merges ERS17727554 and ERS17727555.
- Its "additional details" field says the strain "was originally thought to be *Discosporangium mesarthrocarpum* strain
  KU-2346 … but … rbcL … 99% identical with that of *Choristocarpus tenellus* (KU-1152) and … the COX1 sequence of this
  strain is identical to that of KU-1152".
- ON552980 is from "Villefranche-sur-Mer", collected 1980-09-01. That is the same locality and date as KU-1152.
- So the Choi record and the Phaeoexplorer assembly are very likely the same culture line and not independent isolates.
  Make sure the *Discosporangium* used as the AUA = Ile sister is MT17-79 and not this mislabelled strain.

## Q2. AUA = Met in other stramenopile and non-animal, non-fungal mitochondria

### Stramenopiles (closest first)

1. **Xanthophyceae.** This is the closest precedent, a related PX-clade class. [Correction by Lyra, 14 Nov, from the fact-check: the sister class of Phaeophyceae is Schizocladiophyceae, not Xanthophyceae.]
   - Ehara M, Hayashi-Ishimaru Y, Inagaki Y, Ohama T 1997, *J. Mol. Evol.* 45:119–124 (doi:10.1007/PL00006210).
     ABSTRACT only: "six manifested the expected deviant genetic code where AUA codes for methionine (AUA/Met), but not
     for isoleucine". This is from cox1 of six xanthophytes; *Monodus* and *Ophiocytium* and the eustigmatophytes read
     AUA as Ile.
   - Ehara et al. 1998, *J. Phycol.* 34:1005–1008. FULL TEXT (PDF on Kawai's Kobe University page).
     - Abstract: "all the algae described above bear the universal genetic code, including the assignment for the AUA
       codon". This covers *Chorda filum*, *Colpomenia bullosa*, *Ectocarpus* sp., *Pseudochorda nagaii*, *Undaria
       pinnatifida* and *Chattonella antiqua*.
     - Xanthophyte AUA/AUG ratios at Met were 14/4 to 17/1.
2. **Labyrinthulea (thraustochytrids, *Aplanochytrium*, LAB14 S2/S4) and the uncultivated MAST8b.**
   - Žihala D, Salamonová J, Eliáš M 2020, *Mol. Phylogenet. Evol.* 152:106908 (doi:10.1016/j.ympev.2020.106908;
     corrigendum 2026, 220:108593, not read). ABSTRACT only.
   - Quotes: "labyrinthulean mitochondria also exhibit a sense-to-sense codon reassignment, manifested as AUA encoding
     methionine instead of isoleucine … this change evolved independently in the uncultivated stramenopile lineage
     MAST8b, together with the reassignment of the AGR codons from arginine to serine."
   - Note: NCBI table 23 (*Thraustochytrium*) still has AUA = Ile.
3. **Negative controls in recent ochrophyte papers:**
   - *Leucomyxa plasmidifera* (Ochrophyta), *Open Biol.* 2024, PMC11528492, FULL TEXT: "no departures from the
     standard codon meaning in the *L. plasmidifera* mitochondrion were indicated by a FACIL analysis".
   - Eustigmatophytes, *GBE* 2026, PMC13449250, FULL TEXT: "no indications for codon reassignments were detected in any of
     the eustig organellar genomes".
   - *Plocamiomonas* (Pelagophyceae), *Sci. Rep.* 2026, PMC13582874, FULL TEXT: uses code 4 (UGA = Trp), "a rare
     feature among Stramenopiles". There is no AUA change, though it has "a single tRNA-Ile bearing an unusual UAU
     anticodon".
   - Diatom and *Cafeteria* mitochondria have UGA = Trp only (Shulgina 2022 Table A.7).

### Non-animal, non-fungal, outside stramenopiles

- **Green algae:**
  - all pedinophytes: Barcytė, Žihala & Eliáš 2025, *PLoS Genet.* 21:e1011901, PMC12574857, FULL TEXT: "All
    pedinophyte mitochondria seem to decode AUA as methionine rather than the standard isoleucine";
  - *Pycnococcus provasolii*: Turmel et al. 2010, with UGA = Trp and UUR = stop (listed in the Codetta supplement; not
    read directly).
- **Plastid, not mitochondrial:** peDinoflagellate plastids have AUA = Met (Barcytė 2025 cites the earlier papers).
- **Red algae (Cyanidiophyceae):** our own unpublished lab call. No publication found earlier
  (`commons/2026-11-10-new-stop-codons.md`).

## Q3. Reviews and screens

- **Shulgina Y 2022, Harvard PhD dissertation "The Evolution of Alternative Genetic Codes"** (DASH). FULL TEXT.
  - Table A.7, "All stramenopile, alveolate, and rhizarian (SAR) clades known to use a mitochondrial alternative genetic
    code", lists for stramenopiles:
    - Thalassiosirales UGA → Trp;
    - **Xanthophyceae AUA Ile → Met [Ehara 1997]**;
    - *Cafeteria*/MAST12 UGA → Trp;
    - Labyrinthulea AUA → Met + UUA → stop (± UAR → Tyr, AGR → stop);
    - MAST8b AUA → Met + AGR → Ser.
  - No Phaeophyceae and no *Choristocarpus*.
- **Shulgina & Eddy 2023, *Bioinformatics* 39:btac802 (Codetta)**, PMC9825746. FULL TEXT plus supplementary files 2 and 4,
  downloaded from Europe PMC.
  - The mito benchmark is 82 genomes chosen because they have *reported* alternative codes.
  - The stramenopiles in it are *Thraustochytrium* AF288091, S2/S4, *Aplanochytrium*, MAST8b and the Xanthophyceae list
    above.
  - No brown alga is included, and *Choristocarpus* (released 2024) postdates it.
  - Its only unexpected AUA = Met calls were "Two green algal mitochondrial genomes".
- **Shulgina & Eddy 2021, *eLife* 10:e71402.** It screens bacteria and archaea only (per its abstract and the
  dissertation). NOT READ further.
- **Sengupta, Yang & Higgs 2007, *J. Mol. Evol.* 64:662** (PMC1894752). FULL TEXT searched.
  - AUA → Met in heterokonts is cited only for xanthophytes: "has also been observed (Ehara et al. 1997) in some species
    of yellow-green algae".
  - *Laminaria* and *Pylaiella* are treated as standard.
- **Knight, Freeland & Landweber 2001; Kollmar & Mühlhausen 2017.** NOT READ. Both predate the first *Choristocarpus*
  mitogenome and so cannot list it.
- **NCBI genetic codes page** (wprintgc.cgi, read). No table covers stramenopile AUA = Met. Table 23 keeps AUA = Ile.

## Q4. NCBI taxonomy

- E-utilities efetch, read today. The Browser HTML returned a reCAPTCHA page.
- Each of these taxa has MitoGeneticCode MGCId **1, Standard**, not 4:
  - genus *Choristocarpus* (116064);
  - *C. tenellus* (116065);
  - Choristocarpaceae (116066);
  - Phaeophyceae (2870);
  - Ochrophyta (2696291);
  - *Discosporangium mesarthrocarpum* (376270).

## What was searched

- **Europe PMC REST, all sources including preprints**, with these queries:
  - Choristocarpus;
  - Choristocarpus + mitochondri*;
  - Discosporangium;
  - Discosporangiales + mito;
  - Phaeoexplorer;
  - "GCA_964200185";
  - (Phaeophyceae OR brown alga*) + AUA/ATA + Met/Ile + mito;
  - stramenopile/heterokont/ochrophyt* + mito + AUA + methionine;
  - brown algae mitogenome 2024–2026;
  - "genetic code" + brown algae + mito 2020–2026;
  - Ehara M as author;
  - Shulgina Y as author.
- **Other sources:**
  - Crossref, for the Choi et al. title;
  - WebSearch, standard and extended, on the Choi title, *Choristocarpus* mitogenome and code, Phaeoexplorer and Codetta;
  - NCBI nuccore, taxonomy and Datasets, for the records, assemblies and BioSamples.

## Could not read

- The Cell version of Denoeud 2024 (403) and its Table S2.
- The Žihala 2020 full text and its 2026 corrigendum.
- The Ehara 1997 full text.
- Knight 2001 and Kollmar & Mühlhausen 2017.

A search for absence cannot rule out a thesis, poster or unindexed supplement. The likeliest place for a missed mention
is the unpublished Choi et al. mito paper itself.
