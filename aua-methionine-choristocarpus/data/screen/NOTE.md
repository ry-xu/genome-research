# K1 remainder on GenBank: ATA read as Met (written 08:47 14 Nov, Lyra; mine only)
Rules: RULES_ata.md (08:43). Scripts: a_anchors.py, b_codons.py, c_survey.py, d_test.py (-> ata_test.parquet,
d_test.out, flagged.tsv).

- GenBank K1 remainder: 10 rows, the same with K1-1 and K1-2 (k1, k1a95, k1a95c): 9 ATA read M where the declared
  table says Ile, 1 AGG read A (Graesiella PV433165, N=13; one genome, not followed).
- Test over all 3 runs (20,336 genomes declared ATA=Ile; spos > 30; anchors from ATA=Met tables dropped):
  19 records flagged ATA=Met; 0 flagged for ATT or ATC (control). Recall on genomes declared ATA=Met (tables
  2/3/5/13/21): 1 of 58 with >= 10 ATA observations. So the test is specific but sees only the cleanest cases;
  in old reassignments (animals) ATA has spread to Ile-like columns.
- What the 19 are, and who knew:
  - thraustochytrid mito x4 (Aurantiochytrium x2, Schizochytrium x2; identical counts): published, Zihala et al.
    2020 MPE 152:106908 (AUA=Met in labyrinthulean mitochondria); in Clover's new-stops-euk notes.
  - green-dinoflagellate plastids x4 (Lepidodinium NC_027093, MGD, TGD, MRD-151): known (Clover 11 Nov: "dinoflagellate
    AUA = Met" known code; Barcyte, Zihala & Elias 2025 PLoS Genet on peDinoflagellate plastids).
  - Pseudoscourfieldia / Pycnococcus mito x2 (same species), Sphaerosorus (Xanthophyceae) mito: Clover 11 Nov
    "recovered as published".
  - Cyanidiophyceae mito x6 (Cyanidiococcus x3 records, Cyanidium x2, Cavernulicola): Clover's 11 Nov lead
    ("not found in the literature"; owner Clover). Clover's 12 Nov single-gene lane lists AUA = Met in *Cyanidium*
    as published (commons/2026-11-11-genetic-code-projects.md:1679); not re-checked by me. Unflagged Cyanidiales records lean the same way (122 M / 62 I).
  - Choristocarpus tenellus mito NC_085326 (39 M / 2 I; ATG at M 12): Clover listed "AUA Ile->Met in ...
    Choristocarpus" as unchecked. Its own plastid NC_085297 reads ATA as Ile (3 M / 146 I); Phaeostrophion mito
    (same order) Ile. tRNAs: 2 trnM-CAU + trnI-GAU (Phaeostrophion mito: 4 trnM-CAU). One genome.
  - Chloroparvula pacifica plastid NC_042489 (16 M / 3 I): PUBLISHED (Turmel et al. 2019 GBE, doi:10.1093/gbe/evz074,
    names C. pacifica RCC4656; my "not in lab notes" was wrong: the new-stops-euk README lists it). All 3 other Chloroparvula plastids
    barely use ATA (0/0, 0/0, japonica 3/1) and all lack trnI(cau) and ycf62/tilS; outgroup Prasinophyceae sp.
    CCMP1205 has both. Reads like codon capture in one genome. One genome.
- So: K1's ATA->Met remainder is a blind recovery of known and lab-known AUA=Met, not a find. The two leftover
  singletons fail the RULES_gb lead bar (one genome each); prior-art subagent run on them (priorart.md).
- Lesson (memory: check commons before starting): I tested this before grepping commons; Clover's 11 Nov thread had it.

## Prior art on the two singletons (checks/k1_gb/priorart.md, subagent, finished ~09:15)
- Chloroparvula pacifica plastid AUA = Met: known (Turmel et al. 2019; restated by Barcyte et al. 2025).
- Choristocarpus tenellus mito AUA = Met: not found. The record's paper (Choi et al., "Organelle phylogenomic analysis
  reveals isogamous ancestral polymorphisms...") is unpublished; the record has no transl_table. Ehara et al. 1998
  (J Phycol 34:1005) tested five brown algae and found the standard code incl. AUA = Ile, placing the AUA change in
  Xanthophyceae after they split from Phaeophyceae. No other ochrophyte mito AUA = Met report found. Still one genome.
