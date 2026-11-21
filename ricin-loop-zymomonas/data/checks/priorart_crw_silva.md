# Prior art: natural SRL tetraloop variants (GAAA in Zymomonas, GCAA in P. freudenreichii)

Positions use E. coli numbering. The SRL tetraloop is 2659-2662; E. coli reads GAGA.
The two variants are GAAA (G2661A) in *Zymomonas mobilis* and GCAA (A2660C+G2661A) in *Propionibacterium freudenreichii*.
Searched 2026-11-20. Constraint followed: nothing was submitted, nobody was contacted, and no personal details were sent to any service. NCBI E-utilities were called with no email parameter.

## Verdicts

- **A. CRW (Gutell) conservation.**
  - The CRW diagrams show 2659-2662 as conserved, but not as invariant.
  - In the three-domain diagram (592 sequences, Nov 2001), G2659, A2660, G2661 and A2662 are all drawn upper case, meaning "98+% conserved".
  - In the 3-domain + 2-organelle diagram (930 sequences, Oct 2001), 2660 is lower-case "a" (90-98%) and 2661 is an open circle (<80%). Both reflect the mitochondrial GAAA and CAAA loops.
  - The CRW nucleotide-frequency file lists the bacterial minority directly. In the bacterial set ("B", 237 sequences), 2660 C is 0.42% (1 sequence) and 2661 A is 0.84% (2 sequences).
  - The 2001 CRW bacterial alignment contains **both of our organisms**, with exactly these loops:
    - `Z.mobilis` (AF086792) reads GAAA.
    - `P.freudenr` (Y10819) reads GCAA.
  - So both variants have sat, unremarked, in CRW's public data since 2001.
  - RiboVision gives only Shannon-entropy scores over 123 species. It has no variant list, and no E. coli LSU data are archived.
- **B. Papers or preprints on the SRL of Zymomonas or P. freudenreichii, or on ricin or α-sarcin sensitivity of their ribosomes.**
  - Not found.
  - There is one near-miss. St-Onge et al. 2007 names P. freudenreichii in a 23S alignment, but at other "sarcin-ricin motif" sites, not at the 2659-2662 tetraloop.
  - Lee et al. 2001 (Mol Cells) proposes a Z. mobilis 23S secondary structure. It may show the GAAA, but it could not be opened (paywalled).
- **C. Blanchet et al. 2014.**
  - **No organism is named.**
  - Table 1, row 15 reads `UUAGUACGCAAGGACCG` (Freq 1) and comes from "a 23S rRNA bacterial alignment" supplied by E. Westhof. The paper has no supplementary data file.
  - This 17-mer is identical to P. freudenreichii (GenBank Y10819 from 1997; also the CIRM-BIA1 genome). That identification is our inference; the paper does not state it.
  - There is no GAAA row, so Zymomonas is evidently absent from that alignment, or its SRL is.
- **D. Kudrin et al. 2018.**
  - **No.**
  - The paper uses A2660C/G2661A as an engineered SRL mutant. It gives no rationale, natural organism or citation for the choice.
  - The supplement does not mention it either.
- **E. Other reports of natural bacterial SRL tetraloop variants.**
  - Not found.
  - Doris 2015 used the SRL only as a filter (≤4 mismatches) and lists no variants.
  - The 2015-2026 Europe PMC searches found only:
    - engineered mutants (Lancaster 2008, Shi 2012, Kudrin 2018, Aval 2025);
    - statements that the GAGA tetraloop is "highly conserved";
    - organellar variants (human mito GAAA; Koripella 2020, already known).
  - The only bacterial listings are data rows without comment: CRW 2001 (both organisms) and Blanchet 2014 Table 1 (GCAA, unnamed).

## Evidence

### A. CRW conservation diagrams, frequency file and alignment; RiboVision

All files were retrieved from the Wayback Machine as raw `id_` copies. Copies of the diagrams, SRL crops, frequency file, 23S.B.GB, Blanchet Table 1 and the parser (`gbsrl.py`) are in `work/curation/src/srl_priorart/`.

1. **Three-domain diagram:** `https://crw-site.chemistry.gatech.edu/SAE/2B/ConsStruc/Diagrams/cons.23.3.3DOM.3.pdf` (snapshot 20250329151938).
   - Title: "Phylogenetic conservation superimposed onto the Escherichia coli Large Subunit rRNA secondary structure - 3' half", "Three Phylogenetic Domains", "1.cellular organisms", "November 2001", "Number of sequences: 592".
   - Legend, verbatim: "Positions with a nucleotide in more than 95% of the sequences are shown in one of four categories: ACGU - 98+% conserved / acgu - 90-98% conserved / • - 80-90% conserved / ○ - less than 80% conserved".
   - I rendered the SRL region at 600 dpi and read it. The hairpin ...A G U A C **G A G A** G G A... is entirely upper case, so 2659-2662 are each 98+% conserved. The diagram does not mark the loop "invariant"; the scheme has no separate invariant class.
2. **Three-domain + two-organelle diagram:** `.../cons.23.t.3DOM.2ORG.3.pdf` (snapshot 20250329124328).
   - Header: "Three Phylogenetic Domains / Two Organelles", "October 2001", "Number of sequences: 930".
   - The SRL reads ...a G U A C G **a ○ A** G G A...
     - 2659 G is upper case.
     - 2660 is "a" (90-98%).
     - 2661 is an open circle (<80%).
     - 2662 A is upper case.
   - This matches the mitochondrial GAAA and CAAA loops.
   - Lyra's earlier Mears et al. 2002 CONS.1 table (in `organelle-reading/checks/srl/priorart.md`) gives the same picture: 2661 "w" (<80%) and 2660 lower case in 3P2O.
3. **Nucleotide-frequency file:** `https://crw-site.chemistry.gatech.edu/SAE/2A/nt_Frequency/SB/HARRY_NOLLER_23S_FREQ_0.txt` (snapshot 20210119140919).
   - It gives per-position percentages by alignment set: T = total, 3 = three domains, B, A, C (chloroplast), E, M (mitochondria). "Seq" is the number of sequences covering the position. Verbatim rows:

   ```
   Position: 2660            A      G     C     U    -    Seq
   T   97.96 0.00 2.04 0.00 0.00 737
   3   99.76 0.00 0.24 0.00 0.00 417
   B   99.58 0.00 0.42 0.00 0.00 237
   M   94.66 0.00 5.34 0.00 0.00 262
   Position: 2661
   T   32.02 67.98 0.00 0.00 0.00 737
   3    0.48 99.52 0.00 0.00 0.00 417
   B    0.84 99.16 0.00 0.00 0.00 237
   M   89.31 10.69 0.00 0.00 0.00 262
   ```

   - 2659 and 2662 are 99.6-100% G and A in B. Archaea (41), chloroplasts (58) and eukaryotes (139) are 100% GAGA at 2659-2662.
   - In B, 0.42% of 237 is 1 sequence (C at 2660), and 0.84% of 237 is 2 sequences (A at 2661). That is consistent with exactly one GCAA and one GAAA in the bacterial set.
4. **CRW bacterial alignment (2001):** `http://www.rna.ccbb.utexas.edu/SIM/4B/AA_AG/ALIGNMENTS/23S.B.GB` (snapshot 20100627033237; editor header "Fri May 11 12:33:53 2001"; companion file `23S.B.aln`).
   - I stripped gaps and searched each entry for `AGUAC....GGAC`. 222 entries cover the SRL: 220 read GAGA, 1 GAAA and 1 GCAA. Another 105 lack the region (partial sequences or group dividers), and 3 more give only off-target matches (T. maritima M67498, Pirellula X07408, L. carnosum S60371).
   - The only bacterial entries with a non-GAGA tetraloop at the true SRL were:

   ```
   Z.mobilis   SOURCE Zymomonas mobilis.  GenBank AF086792   ...CCUUAGUACGAAAGGACCGGG   (GAAA)
   P.freudenr  SOURCE Propionibacterium freudenreichii. subsp. ITS shermanii  GenBank Y10819   ...CCUUAGUACGCAAGGACCGGG   (GCAA)
   ```

   - AF086792 is Lee, Jin, Kang & Kang, submitted 26-AUG-1998, "Zymomonas mobilis 16S ribosomal RNA gene, complete sequence" (6,250 bp operon). I confirmed GAAA in it via efetch.
   - The CRW alignment carries no annotation of either variant.
5. **Not opened:** the newer CRW total alignment `DAT/3C/Alignment/Files/CRWformat/23S.T.crw` (snapshot 20210416). The Wayback copy is truncated at 1,048,576 bytes and is not a valid zip.
6. **RiboVision** (`apollo.chemistry.gatech.edu/RiboVision`):
   - The live host does not respond.
   - In the Wayback copy of `Documentation/AlignmentData.html` (snapshot 20141117192436), conservation is described as: "a multiple sequence alignment from a representative set of 123 species from each domain of life ... The sequences were aligned with SILVA. The scores reported are Shannon Entropy scores with gaps penalized."
   - These are per-position scores only, with no variant table. Only Drosophila (`DM_LSU_3D`) residue data are archived; there is no E. coli LSU entry in the CDX listing. The ribovision2 site responds but was not queried, since it serves the same entropy-style data.

### B. Literature on the SRL of Zymomonas or P. freudenreichii, and toxin sensitivity

**Europe PMC full-text queries.** Hit counts are given in parentheses. Every hit was checked by title, and full text was grepped where relevant.

| Query | Hits |
|---|---|
| `freudenreichii AND sarcin` | 4 |
| `freudenreichii AND ricin` | 3 |
| `freudenreichii AND "23S rRNA"` | 42 |
| `Zymomonas AND sarcin` | 3 |
| `Zymomonas AND ricin` | 10 |
| `Zymomonas AND "23S rRNA"` | 78 |
| `Propionibacterium AND sarcin` | 6 |
| `Propionibacterium AND ricin` | 26 |
| `(Zymomonas OR freudenreichii OR Acidipropionibacterium) AND ("sarcin" OR "SRL" OR "H95" OR "helix 95")` | 66 |
| `(Zymomonas OR freudenreichii) AND ("ribosome-inactivating" OR "ribosome inactivating" OR restrictocin OR "alpha-sarcin" OR "α-sarcin")` | 13 |

- None of these hits reports either organism's SRL sequence or tests ricin, α-sarcin or restrictocin on its ribosomes.
- The Zymomonas "sarcin" and "ricin" hits are unrelated. For example:
  - PMC5633077 is the human mitoribosome.
  - PMC1364539 is a GM-crop review in which Zymomonas appears only for aroA.
  - PMC8272150 covers antimicrobial-peptide expression.

**Near-miss: St-Onge K, Thibault P, Hamel S, Major F (2007).** "Modeling RNA tertiary structure motifs by graph-grammars." NAR 35:1726, doi:10.1093/nar/gkm069, PMID 17317683, PMC1865062.
- The paper scans a Westhof bacterial 23S alignment of 806 sequences for the bulged-G "sarcin-ricin motif" at several sites (near loops L11, L13, L21 and L23, plus the last site).
- P. freudenreichii is named twice, both times at other sites:
  - "Another possible sequencing error would be at position 29 in the sequence of Propionibacterium freudenreichii, where a G would fit better with all other 805 sequences."
  - "MC-Fold supports a sarcin–ricin at the same position in the sequence of P. freudenreichii, but keeps the gap at position 40."
- Neither concerns 2659-2662. This shows a P. freudenreichii 23S was in the Westhof alignment.

**Not opened: Lee JS, Jin SJ, Kang HS (2001).** "Molecular organization of the ribosomal RNA transcription unit and the phylogenetic study of Zymomonas mobilis ZM4." Mol Cells 11:68-74, PMID 11266123, doi:10.1016/s1016-8478(23)17005-1.
- The abstract says: "The secondary structure models of 16S, 23S and 5S rRNA are proposed."
- The ScienceDirect full text returned 403, so it is unknown whether the GAAA is drawn or discussed.

**Public sequence history.**
- Zymomonas GAAA is in:
  - AF086792 (submitted 26-AUG-1998);
  - AF088897 (cosmid 65G3, Lee & Kang, 31-AUG-1998);
  - AF117351 (rrnB/rrnC, 02-JAN-1999);
  - the ZM4 genome AE008692 (2005, 3 copies).
- P. freudenreichii GCAA is in Y10819 (subsp. shermanii partial 16S, 23S and 5S; Dasen, Smutny, Teuber & Meile, Syst Appl Microbiol 21:251-259, 1998; submitted 27-JAN-1997), reading `GTCCTTAGTACGCAAGGACCGGGACGG`.
- Older P. freudenreichii records lack the SRL region: M85129, AF280074/5 and AF218432/4.

### C. Blanchet et al. 2014

Blanchet MF, St-Onge K, Lisi V, Robitaille J, Hamel S, Major F (2014). "Computational identification of RNA functional determinants by three-dimensional quantitative structure-activity relationships." NAR 42:11261, doi:10.1093/nar/gku816, PMID 25200082, PMC4176186.

**Table 1 "SRL variants".** This is an image; I read it from the PMC figure file. The rows span E. coli 2652-2668.

- **Row 15, verbatim:** `15 | UUAGUACGCAAGGACCG | Freq 1 | (no growth data) | 2D dark gray | Score 1.00`
- **Footnote:** "the middle region eight variants from a 23S rRNA bacterial alignment ... Freq: frequencies of occurrence of each sequence in the alignment."
- **Methods:** "Eight variant sequences were extracted from a 23S rRNA bacterial alignment, IDs #12 to #19."
- **Acknowledgements:** "We also thank Eric Westhof for providing the bacterial 23S rRNA alignment."
- **Organism:** not named in the text, table, footnote or legend. PMC lists no supplementary data. The Europe PMC supplementaryFiles zip holds only figure and table images. The MC-3DQSAR website (major.iric.ca) was not checked.
- **Frequency rows, verbatim:**
  - E. coli CUAGUACGAGAGGACCG 457; #02 2; #05 1;
  - #12 UUAGUACGAGAGGACCG 271; #13 CUAGUACGAGAGGACCA 72; #14 AUAGUACGAGAGGACCU 6;
  - #15 UUAGUACGCAAGGACCG 1; #16 1; #17 1; #18 1; #19 1.
  - The total is 814, close to the 806 in St-Onge 2007.
  - There is **no GAAA row**.
- **Our inference, not stated in the paper:** row 15 is P. freudenreichii.
  - It is identical over all 17 nt to P. freudenreichii Y10819 and to the CIRM-BIA1 genome (GCA_000091725.1, both copies `TTAGTACGCAAGGACCGG`).
  - The same lab's earlier paper (St-Onge 2007) used a Westhof bacterial alignment that included P. freudenreichii.
  - The 2001 CRW bacterial alignment has GCAA only in P. freudenreichii.

### D. Kudrin et al. 2018

Kudrin P, et al. (2018). "The ribosomal A-site finger is crucial for binding and activation of the stringent factor RelA." NAR, doi:10.1093/nar/gky023, PMID 29390134, PMC5829649.

**Full text grep.** Zero hits for: Blanchet, Marchant, phylogen, Propionibacter, Zymomonas, mitochondri. No reason is given for choosing GCAA. Quotes:
- "A2660C/G2661A mutant in the SRL of domain VI, helix H95."
- "In the case of RelA, both the A2660C/G2661A SRL mutation and alpha-sarcin cleavage inhibit activation by the 70S alone (Figure 4B). However, activation by 'starved' complexes is retained…"
- "the A2660C/G2661A mutant grows normally at 32–42°C, is cold-sensitive and displays a minor growth defect at 30°C."
- Mutagenesis primer: 5′-gggctgctcctagtacgcaaggaccggagtggacgc-3′.

**Supplement** (gky023_supp.docx via the Europe PMC supplementaryFiles zip): Supplementary Figure 1 (translocation kinetics) only. It has no GCAA and no natural organism.

### E. Other reports of natural bacterial SRL tetraloop variants

**Doris et al. 2015** (PMC4574749) used the SRL as a filter: "Accessions that did not contain the 14-nt sarcin–ricin loop (SRL) AGUACGAGAGGAAC sequence at least 70% conserved (i.e., ≤4 mismatches) at the appropriate structural position at the 3′ end of 23S–28S rRNA were eliminated." It also writes "The SRL coincides with uCNE10, which is conserved in ≥90% of rRNA sequen[ces]". It lists no individual SRL variants.

**Europe PMC queries, 2015-2026:**

| Query | Hits |
|---|---|
| `"sarcin-ricin loop" AND (natural OR naturally) AND (variant OR variants OR polymorphism) AND PUB_YEAR:[2015 TO 2026]` | 125 |
| `("sarcin-ricin" OR "sarcin/ricin") AND (GAAA OR GCAA) AND PUB_YEAR:[2015 TO 2026]` | 22 |
| `("sarcin-ricin loop" OR "α-sarcin loop" OR "alpha-sarcin loop") AND (phylogen* OR "across bacteria" OR "sequence variation") AND PUB_YEAR:[2018 TO 2026]` | 102 |
| `"sarcin-ricin" AND (G2661A OR A2660C OR "G2661" OR "A2660") AND PUB_YEAR:[2010 TO 2026]` | 48 |
| `SRC:PPR AND ("sarcin-ricin loop") AND (variant OR natural OR GCAA OR GAAA)` | 1 (irrelevant: SARS-CoV-2 RetroMAD1) |
| Earlier: `GCAA AND "sarcin-ricin"` | 22 |
| Earlier: `GAAA AND "sarcin-ricin"` | 42 |
| Earlier: `GCAA AND "sarcin"` | 24 |

I inspected the top 30 titles for each query. I grepped the full text of the relevant ones for GCAA, GAAA, G2661, A2660, "natural variant" and tetraloop:

- **Aval et al. 2025** (PMC11912913): engineered mutants only (G2655C, A2660U, ΔU2653ΔC2667; Δ2653–2667::GAAA or ::CUUU). The "GAAA" there is a Lancaster 2008-style replacement loop, not a natural variant.
- **"Ribosome-inactivating proteins in bacteria"** (Microb Genom 2026, doi:10.1099/mgen.0.001839, PMC13592729) says only: "The sarcin–ricin loop is highly conserved, but a major factor underlying the different activity levels can be the nature of the RIPs' interactions with ribosomal proteins". It has no bacterial SRL variant.
- **Spiroplasma toxin-arsenal paper 2025** (PMC12422007): "Nucleotide identity is highly conserved across the SRLs present in this study and across domains of life with regards to the GAGA tetraloop regions". It covers insect hosts only.
- **Organellar variant already known:** human mitochondrial GAAA (Koripella 2020, PMC7395135).
- **Earlier report** (`work/frozen-core/results/lsu/bac/PRIORART.md`) covered Marchant & Hartley 1994/1995, Lancaster 2008, Shi 2012 and Doris 2015, and found no natural bacterial or archaeal G2661/A2660 variant described in text.

**Net result.** The only bacterial records of natural SRL tetraloop variants found are uncommented data:
- CRW 2001 alignment and frequency tables: Z. mobilis GAAA and P. freudenreichii GCAA, with frequencies 1-2/237;
- Blanchet 2014 Table 1 row 15: GCAA, unnamed, Freq 1.

No text was found that names either organism as carrying a variant SRL tetraloop.

### Not opened / limits

- Lee et al. 2001 Mol Cells full text (403).
- The 23S.T.crw zip (truncated Wayback copy).
- major.iric.ca/MC-3DQSAR.
- St-Onge 2007 Figure 5, which might show whether its "last" site is H95.
- ribovision2 data.
- Searches covered Europe PMC (including the preprint source PPR), NCBI nuccore and Wayback copies of the CRW and RiboVision sites. bioRxiv was not searched separately; Europe PMC indexes bioRxiv preprints as PPR. Google Scholar and paywalled full texts were not searched.
