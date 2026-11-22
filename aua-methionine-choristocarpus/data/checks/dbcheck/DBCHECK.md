# Reference-database check: Choristocarpus mito AUA = Met (draft outbox/2026-11-14-choristocarpus-aua-met/post.md)

Fresh, independent check of the reference databases (not papers). Run on 2026-11-14 from 15:22 to 15:33. No email or tool
parameters were sent, and no accounts were used. Downloads are in `data/`, scripts are `d*.py` / `*.sh`, and every number below
comes from the matching `.out` file (or, for the tRNA finders, from the raw outputs in `trna/`).

## Bottom line
No database records AUA = Met (or any non-standard mitochondrial code) for *Choristocarpus*, *Discosporangium* or any
brown alga. Every record that encodes a code uses the standard code (table 1), and that is consistent with the draft. Real tRNA-finder runs
(tRNAscan-SE 2.0.13, ARAGORN 1.2.41) agree with the draft's similarity-based tRNA statements. Three things need qualifying (see the last section):
- The NCBI code assignment carries no information. NCBI also gives the standard mito code to classes where AUA = Met is published.
- No existing NCBI translation table is "standard + AUA = Met".
- tRNAscan-SE labels *Discosporangium*'s extra CAU tRNA "Ile2".

## 1. NCBI Taxonomy
- **Queries:**
  - esearch db=taxonomy, `<name>[Scientific Name]` for Choristocarpus tenellus, Choristocarpus, Choristocarpaceae, Discosporangium,
    D. mesarthrocarpum, Discosporangiaceae, Discosporangiales and Phaeophyceae. Also `[All Names]` for Discosporangiales and Discosporangiophycidae (fetch1.sh, fetch3.sh, d14).
  - efetch db=taxonomy id=116065,116064,116066,376270,376269,376268,2870 retmode=xml (fetch2.sh, d1).
  - NCBI Datasets v2 `taxonomy/taxon/116065,2870,2833,35131,240563/dataset_report` (d2a).
- **Findings:**
  - MitoGeneticCode = 1 (Standard) for all 7 taxids: *C. tenellus*, *Choristocarpus*, Choristocarpaceae, *D. mesarthrocarpum*,
    *Discosporangium*, Discosporangiaceae and Phaeophyceae (d1). Nuclear GC = 1. The Datasets plastid code = 11.
  - **Comparison taxa:** NCBI also gives mito code 1 to Xanthophyceae (2833), *Sphaerosorus* (240563) and Labyrinthulomycetes
    (35131) (d2a). These are the groups the draft cites as published AUA = Met.
  - **"Discosporangiales":** 0 hits as Scientific Name and 0 as All Names (d14). NCBI places both genera in **Sphacelariales**
    (order 27945) (d1).
  - **Inherited / flagged:** not checked. The efetch XML and the Datasets JSON carry no inheritance flag; that is only in
    nodes.dmp, and I did not download the taxdump (about 65 MB). The species record was last updated 2005/01/19.
- **Verdict:** NCBI assigns the standard mito code throughout. It does the same for known AUA = Met stramenopile groups, so the
  assignment is a default, not a curated judgement.

## 2. GenBank ON552980 / RefSeq NC_085326, Genoscope contig, and tRNA finders
- **Queries:** efetch db=nuccore rettype=gbwithparts/fasta for ON552980, NC_085326, CAXNAD020000071.1 (Genoscope *Choristocarpus*
  KU-1152 "Cten_MITOCHONDRIA") and CAXNAE020000396.1 (*Discosporangium* MT17-79 WGS contig). Controls: NC_013473.1 (*Saccharina*),
  ON552993.1 (*Phaeostrophion*) and NC_007684.1 (*Desmarestia*) (fetch1/2.sh, d3).
- **Record findings (d3, d13):**
  - **CDS and translations:**
    - All 35 CDS in each *Choristocarpus* record (ON552980, NC_085326, CAXNAD020000071.1) lack /transl_table, so they use table 1.
    - 0 CDS have a /translation that differs from table 1.
    - The CDS contain 227 ATA codons (229 in the Genoscope contig). All are translated as Ile.
  - **Text:** no /note on any feature of ON552980 or NC_085326. 0 hits for "genetic code", "AUA", "isoleucine", "methionine" or "code".
    The only COMMENT is NC_085326's "PROVISIONAL REFSEQ … identical to ON552980".
  - **Codons:** all 35 CDS start with ATG; stops are TAA 26, TGA 6, TAG 3 (NC_085326).
  - **tRNA features:**
    - ON552980 and NC_085326 have 22 tRNA features with no /anticodon qualifiers.
    - Gene names include **two trnM-CAU** (2625..2701, 2859..2931) and one **trnI-GAU** (18862..18934). There is no trnI-CAU.
    - The Genoscope contig has 23 tRNA features, with names only in locus_tags.
    - The same submitters label *Phaeostrophion*'s canonical trnI(CAU) (ON552993 22257..22330) "trnM-CAU". GenBank CAU labels
      from this group therefore cannot distinguish Ile-CAU from Met.
- **tRNA finders:**
  - **Tools and runs:**
    - ~/dev/tools/bin/tRNAscan-SE 2.0.13, run in -O (organellar), -G and -B (bacterial isotype models) modes, all with --thread 1.
    - aragorn1.2.41 -t -gc1, default and -ps.
    - Runs are in d4_run_trna.sh and d5_run_trnaB.sh, tabulated in d6, d7 and d8.
  - **Counts:**
    - *Choristocarpus*: 21–23 tRNAs per tool.
    - *Discosporangium*: 23–25.
    - *Saccharina*: 24–25.
  - ***Choristocarpus* (both assemblies):**
    - Every tool finds exactly 2 CAT-anticodon tRNAs, 1 GAT (Ile) and no TAT.
    - Neither CAT tRNA is the brown-algal trnI(CAU):
      - **Slot.** In all 3 controls, trnI(CAU) sits between trnW and trnQ. In *Choristocarpus*, trnW → trnQ has no tRNA between them.
      - **First CAT (2626–2700).** Full-length identity is 63/75 = 0.840 to the *Saccharina* initiator Met, and 51/73 = 0.699 to the
        *Desmarestia* trnI(cau). It is the fMet-like copy.
      - **Second CAT (2859–2931).** It sits between trnL(UAA) and trnH. Best identity is 51/73 = 0.699 to the *Desmarestia* Met-elongator,
        and 40/73 = 0.548 to trnI(cau).
      - **Isotype calls.** tRNAscan -B calls both "Met", with best isotype Pro and Lys respectively.
  - ***Choristocarpus* weak extra hit:** tRNAscan -O alone reports "Ile AAT" (18773–18841, score 20.2) in both assemblies.
    - -G, -B and ARAGORN do not report it.
    - The anticodon sits at 31–33 with a malformed anticodon stem.
    - Its best match to any tRNA in the 6 genomes is 0.640 (d8).
    - It is not a credible tRNA. With an inosine at position 34 it could in principle read AUA, so it is noted here. It was not checked against RNA.
  - ***Discosporangium* contig:**
    - **Slot.** It also has no CAT tRNA between trnW (21715) and trnQ (21921), so the canonical slot is empty, as the draft says.
    - **Copies.** It has 3 CAT tRNAs:
      - 16956, fMet-like: 65/75 = 0.867 to the *Saccharina* initiator.
      - 33176, Met-elongator-like: 52/73 = 0.712.
      - 33023, the "extra divergent" one: 49/73 = 0.671 to both *Saccharina* Met-e and *Desmarestia* trnI(cau).
    - **Isotype call.** tRNAscan-SE -B labels the extra one **"Ile2"**, the CAU-anticodon Ile isotype that reads AUA (isotype scores:
      Ile2 64.5, Met 57.8, best Lys 67.2).
    - **Calibration caveat.** -B labelled 0 of the 3 control trnI(CAU)s "Ile2". It called them Met, with best isotype Arg. So the label
      is not calibrated for brown algal mitochondria.
  - **Controls:** all 3 have trnM-e, trnI(CAU) and trnM-f by synteny and identity. The *Desmarestia* GenBank annotation names them
    trnM(cau)e, trnI(cau) and trnM(cau)f.
- **Verdict:**
  - The records are translated with the standard code, carry no notes, and have no trnI-CAU annotated.
  - Real tRNA finders confirm that *Choristocarpus* has only fMet-type and elongator-type CAU tRNAs, with no homologue of trnI(CAU).
  - They also confirm the *Discosporangium* pattern the draft describes.

## 3. Organellar tRNA / code databases
- **mitotRNAdb:** http://mttrna.bioinf.uni-leipzig.de now redirects to tRNAdb 2026 (tdb.bioinf.uni-leipzig.de). I used its public
  JSON search endpoint, POST /browse (tdb_q1–3.json, d9, d10).
  - *Choristocarpus tenellus*: 13 entries, all from NC_085326 (mt-tRNA_DNA). One CAT entry, labelled Met. No Ile entry at all,
    not even the GAU tRNA, so the set is a filtered subset.
  - *D. mesarthrocarpum*: 0 entries. The genus is absent from the organism list.
  - Clade Phaeophyceae with anticodon CAT: 61 entries from 61 organisms, every one labelled Met. It keeps one CAT entry per
    organism and never labels a brown-algal CAU tRNA as Ile.
  - **Verdict:** nothing about AUA or the code. It cannot distinguish Ile-CAU, so it neither confirms nor contradicts.
- **GtRNAdb (nuclear-focused):** browse.html lists 610 eukaryotic genome links. None is a brown alga, a stramenopile or
  *Choristocarpus*; the only "Phaeo" hits are bacteria (d11). **Verdict:** nothing relevant, as expected for a nuclear tRNA database.
- **NCBI genetic code tables (gc.prt and wprintgc page; d2b–e):**
  - Tables with AUA = Met are 2, 3, 5, 13 and 21.
  - The "Systematic Range" text names vertebrates, yeasts, invertebrates (nematodes, molluscs, arthropods), ascidians and trematodes.
    There are 0 hits for Phaeophyceae, brown alga, stramenopile, Xanthophyceae, yellow-green, Labyrinthul or Ochrophyta.
  - Every AUA = Met table also reads UGA as Trp. No table differs from table 1 only at AUA.
  - *Choristocarpus* CDS end in TGA 6 times, so none of the existing tables fits a "standard + AUA = Met" reading.
  - **Verdict:** no NCBI code lists brown algae (or any stramenopile) under AUA = Met. Re-annotating with an existing table would be wrong.

## 4. Protein databases
- **NCBI Protein:**
  - **Queries:** efetch rettype=gp for the 35 WBP69798–WBP69832 (ON552980) and the 35 YP_011007967–YP_011008001 (NC_085326) (d12).
  - **Sequences:** all 70 equal the table-1 translation of their CDS; 0 equal the AUA→Met translation.
  - **Text:** no transl_table and no code text. The /notes are only automatic CDD domain notes. The YP COMMENT is "PROVISIONAL REFSEQ …
    Method: conceptual translation".
- **UniProt:**
  - **UniProtKB:** organism_id:116065 returns 0 entries, organism_id:376270 returns 0, and a free-text "Choristocarpus" search returns 0.
  - **Proteomes:** 0 for either taxid.
  - **UniParc:** dbid query for the 70 protein ids returns 35 UPIs. None has a UniProtKB accession. All 35 sequences equal the
    table-1 translation.
- **Verdict:** every protein database copy uses Ile at AUA, with no curator notes. There is no UniProtKB entry or proteome.

## Draft sentences the database results contradict or should qualify
1. **"…and NCBI assigns *Choristocarpus* the standard code."** True (MGC 1, d1). As placed, though, it reads as evidence that brown
   algae "have been taken to use the standard code". NCBI gives the same standard mito code to Xanthophyceae, *Sphaerosorus* and
   Labyrinthulomycetes, where AUA = Met is published (d2a). Suggest adding "as it does for the yellow-green algae", or dropping it.
2. **"*Choristocarpus tenellus* is a small filamentous brown alga in Discosporangiales…"** Fine, but NCBI has no taxon named
   Discosporangiales (0 hits, d14) and files both genera under Sphacelariales. The note to Ryan says this; the post does not. A
   reader checking NCBI will not find the order.
3. **"Brown algae normally have an isoleucine tRNA with anticodon CAU… I found no such tRNA in *Choristocarpus*." / "These are
   similarity calls, not tRNA-finder runs."** Now confirmed by tRNAscan-SE (3 modes) and ARAGORN (2 modes) on both assemblies,
   with a synteny check (empty trnW–trnQ slot) and controls (3 of 3 standard brown algae have it) (d6, d7). The last sentence can be
   updated.
   - Add that GenBank's two "trnM-CAU" labels are not evidence either way: the same submitters call *Phaeostrophion*'s trnI(CAU)
     "trnM-CAU".
   - Optional: one low-scoring (20.2) "Ile-AAU" hit from tRNAscan -O alone is not supported by the other runs (d8). A reader who
     reruns tRNAscan -O will see it.
4. **"It is also missing in *Discosporangium*… *Discosporangium* keeps one extra, divergent CAU tRNA that *Choristocarpus* lacks."**
   Consistent with the finder runs: the canonical slot is empty, and there are 3 CAU tRNAs against 2 in *Choristocarpus*. But
   tRNAscan-SE -B classifies that extra *Discosporangium* tRNA as **Ile2**, the AUA-reading Ile isotype (d6). The label is weak
   (0/3 calibration on control trnI(CAU)s), but it is the obvious candidate for how *Discosporangium* reads AUA as Ile, and could be
   said so.
5. **"…translated with the standard code."** Confirmed: 35/35 CDS have no transl_table, 35/35 NCBI proteins and 35/35 UniParc
   sequences are table-1 translations, and 227 internal ATA are read as Ile (d3, d12). Optional: say the protein databases
   (NCBI Protein, UniParc) carry these Ile-at-AUA sequences, and that no UniProtKB entry exists.
6. **Possible addition.** No NCBI translation table encodes "standard + AUA = Met"; every AUA = Met table also makes UGA = Trp
   (d2e), and *Choristocarpus* uses TGA as a stop 6 times. A correct database annotation would need a new table. This is useful if
   the post or a GenBank update request suggests re-annotation.

No database entry contradicts the core claim, and none pre-records it.

## Not checked
- The NCBI MGC inheritance flag (needs nodes.dmp).
- The Phaeoexplorer portal's own annotation.
- ENA's view of CAXNAD/CAXNAE beyond the GenBank mirror.
- BOLD (the cox1 barcodes AB899175/6 and LM995414 were not opened).
- RNA coverage of the weak AAT hit.
