# Polyangiaceae TTA-stop clade: catalogue search log

Target clade: f__Polyangiaceae genera g__JADLFU01, g__JANWXX01, g__DASTPP01, g__JAXLYQ01, g__CAYNEK01 (GTDB R232).
Searched 2026-11-11, 10:45 to 12:10 EST. Scripts are in `src/catalogues/`, tables in `results/catalogues/`,
and downloads in `data/catalogues/`.
Result: **9 new candidate members plus 1 redundant copy** (`members.tsv`). These are candidate clade members by
marker AAI. None has been tested for TTA use.

**Update, 2026-11-11 16:44 EST (Basil).** The NCBI screen was resumed at 12:39 and finished at 16:26: all 4,477
selected assemblies were screened, 4,441 scored and 36 too small to score (20–51 kb each, no markers; 0 download
failures). The frozen rule now gives 12 NCBI assemblies: the 3 below plus GCA_060838125.1, GCA_060838225.1,
GCA_060838245.1, GCA_060894045.1 (wastewater, Virginia), GCA_025453615.1 (wetland, California), GCA_027451725.1
(mine-drainage sediment, China), GCA_947450785.1 (lake, Czech Republic), and GCA_004144535.1 plus GCA_004144565.1
(poplar leaves, Norwich). GCA_060838125.1 is the same strain as GCA_060838225.1. GCA_004144565.1 is probably the rest of
GCA_004144535.1, one genome split into two bins (same four SIP samples; 7 shared clade-like markers where ~48 would be
expected; `scratch/factcheck3e/REPORT.md`). The TTA tests on all candidates (Arm B, checks 2 and 3) are in
`RESULT.md` (project root), "Outside GTDB", and `results/clade_members/cat_members.tsv` (with `cat_members_NOTE.md`); the candidates table below is the
12:10 one and is not updated.

## R232 cutoff
- The GTDB release notes give no sequence cutoff date. R232 (R11-RS232) was released 2026-04-15.
- The latest `ncbi_seq_rel_date` in `bac120_metadata.tsv.gz` (R232) is **2025-09-04**, so that date is used as the cutoff.
- An assembly also counts as "outside R232" if it was released before the cutoff but is absent from the R232 metadata.
  Absence is checked by base accession, GCA or paired GCF.
- In the R232 tree (`tree_neighbours.py`), the five genera form one clade: support 100, 8 species
  representatives. Its sister is a 251-tip Polyangiaceae subtree with support 69.

## Screen (the single sequence test)
`screen.py`, `driver.py` and `catalogue_screen.py` implement a competitive marker-AAI test.

1. **Gene calls.** pyrodigal predicts genes on contigs of at least 500 bp. Genomes with less than 100 kb in such contigs are
   not scored.
2. **Markers.** There are 437 core markers. They are proteins of GCA_020849655.1 that are at least 100 aa long and are
   single-copy in at least 85% of a 94-genome panel. Diamond blastp settings: e ≤ 1e-10, query and subject cover ≥ 50%.
   The panel is:
   - the 10 clade genomes;
   - 84 outgroup genomes, one R232 species representative for each non-clade Polyangiaceae genus
     (`outgroup_panel.tsv`).
3. **Scoring.** For each marker, the candidate's best protein is taken. The test reports:
   - `top_hit_clade_frac`: the fraction of markers whose top-scoring panel hit is a clade genome;
   - `clade_aai` and `out_aai`: mean identity to the nearest clade genome and to the nearest outgroup genome. Only panel
     genomes sharing at least half of the candidate's markers (and at least 20) count.
   A genome is never scored against itself.
4. **Call rule.** The rule was fixed from the controls before any candidate was scored:
   - **member** = n_markers ≥ 50, `top_hit_clade_frac` ≥ 0.5, and `clade_aai` > `out_aai`;
   - **borderline** = 0.2 ≤ `top_hit_clade_frac` < 0.5. There were none.
5. **Redundancy check.** `redundancy.py` uses pyskani:
   - ANI of each candidate to the 10 references, where ANI ≥ 99.5% means a redundant copy. The same BioSample also counts
     as redundant.
   - Pairwise ANI among all candidates. No pair of candidates reached ≥ 80% ANI.
   - pyskani reports nothing below about 80% ANI, so it is useless for calling genus membership. That is why the main test
     is AAI-based (`skani_calib.py`).

**Controls** (`controls_*.tsv`):

| Control | Genomes | `top_hit_clade_frac` | Passes call rule |
|---|---|---|---|
| Positive, leave-self-out (FAST mode) | 10 | 0.76 to 1.00 | all 10 |
| Positive, leave-genus-out (all genomes of its genus removed from the panel) | 10 | 0.76 to 0.93 | all 10 |
| Negative, 84 panel outgroups, leave-self-out | 84 | max 0.079 | 0 |
| Negative, held-out reps of 35 other c__Polyangia families and 15 other Myxococcota families (`neg_controls_heldout.tsv`) | 50 | max 0.141 | 0 |

- In the leave-genus-out control, JAXLYQ01 (0.76) is the single-genome genus.
- In the held-out negatives, 3 of 50 had `clade_aai` > `out_aai`, but all had a low fraction.
- The gap between 0.14 (highest negative) and 0.76 (lowest positive) is wide. Every candidate scored 0.80 or more.
- Among all non-candidate genomes screened from every source, the highest fraction was 0.14, an NCBI genome, at 12:10.
  After the full NCBI screen (16:26) it is 0.192 (GCA_061528845.1, 26 markers); among genomes with ≥ 50 markers, 0.189
  (GCA_012031025.1). Later candidates scored as low as 0.53 (GCA_004144565.1, see the update at the top), so
  "every candidate scored 0.80 or more" held only at 12:10.
- **Scale for reading AAI.** In the leave-genus-out control, each reference's marker AAI to the nearest *other* clade
  genus was 63.7–66.5%. Its AAI to the nearest outgroup (always GCA_903893335.1) was 59.7–62.1%.

## Candidates (`members.tsv`) and how to read them

| Candidate | Source | Nearest clade genome (genus) | Marker AAI | Fraction | Reading |
|---|---|---|---|---|---|
| GCA_982529355.1 | NCBI | JANWXX01 | 69.9% | 0.94 | above the between-genus level; a predatory soil bacterium from Finland, grown in co-culture with its prey (not a pure isolate; corrected by Basil 12:17 from the genome note, F1000Research 14:1153 v2) |
| GCA_057137765.1 | NCBI | GCA_963169165.1 (JADLFU01) | 89.7% | 0.99 | ANI 86.9%; paddy-soil MAG from China |
| GCA_035762045.1 | NCBI, = SMAG S400bin.4 | GCA_035762175.1 (DASTPP01) | 85.0% | 0.85 | farm soil, China; released 2024-01-17 but absent from R232, probably failed QC (3,387 contigs, 49–56% complete); only 268 markers |
| spire_mag_00671464 | SPIRE | GCA_964605685.1 (JANWXX01) | 96.4% | 0.93 | ANI 94.7%, close to the species boundary; Yellowstone hot spring; only 50.7% complete |
| MAG.1207 | AnMAGC | JANWXX01 | 79.6% | 0.99 | |
| bin242 | Brown23 | JANWXX01 | 79.4% | 0.97 | |
| bin474 | Brown23 | JAXLYQ01 | 85.7% | 0.91 | |
| bin496 | Brown23 | JANWXX01 | 63.7% | 0.85 | out_aai 60.6% |
| spire_mag_00879468 | SPIRE | JANWXX01 | 63.3% | 0.80 | out_aai 59.0%; wastewater |
| MGYG000501200 | MGnify soil | DASTPP01 | 100% | 0.99 | redundant copy of GCA_035762175.1 |

- MAG.1207 and bin242 do not match each other at ≥ 80% ANI.
- bin496 and spire_mag_00879468 sit at or below the clade's between-genus AAI. They may be new genus-level lineages.
- A lineage branching just *outside* the clade's MRCA would also score high on this test. These two are the ones most in
  need of a GTDB-Tk placement (de novo or R232 reference tree) to confirm they fall inside the clade.

**Per-genome AAI profile** (`candidate_aai_profile.tsv`, from `candidate_profile.py`): marker AAI of each candidate to each
of the 10 clade genomes, with the 10 references alongside for comparison.
- **Caveat on the outgroup column.** That table requires only ≥ 20 shared markers. So the low-quality outgroup
  GCA_001312805.1 (Byssovorax, about 140 markers) shows an inflated AAI (62–65%) for every genome, references included.
  Ignore it; the screen itself excludes it by the half-of-markers rule.
- **Clearly inside the clade**, nearer one genus than to the others:
  - GCA_982529355.1: 68.3–69.9% to the four JANWXX01 genomes, against 63.2–65.6% to the other genera. This looks like a
    sister of JANWXX01, or a deep member of it.
  - GCA_057137765.1: JADLFU01, 86.9–89.7%.
  - spire_mag_00671464: JANWXX01, 82–96%.
  - MAG.1207 and bin242: JANWXX01, 77–80%. Within-JANWXX01 AAI among the references is 80–89%.
  - bin474: JAXLYQ01, 85.7%.
- **Flat profiles:**
  - bin496: 61.0–63.7% to all five genera, 60.6% to GCA_903893335.1.
  - spire_mag_00879468: 60.0–63.3% to all five genera, 59.0% to GCA_903893335.1.
  - These look like the profile of the reference JAXLYQ01 genome: 60.2–63.8% to the other genera, 59.7% to
    GCA_903893335.1. JAXLYQ01 is inside the clade in the R232 tree.
  - So AAI cannot tell "a new deep genus inside the clade" from "sister to the clade" for these two. Treat them as the
    least certain candidates.
- **Tree check, unfinished.** A concatenated tree was started with `marker_tree.py`: 150 markers, mafft, iqtree3 LG+G4
  `-fast`, 1 thread. The machine's load average was about 70, so it was stopped at 12:04 during model estimation.
  - The alignment and checkpoint are in `data/catalogues/tree_check/` (`concat.faa`, `tree.ckp.gz`).
  - To finish it, run with more threads:
    `iqtree3 -s data/catalogues/tree_check/concat.faa -m LG+G4 -T 4 -fast --alrt 1000 --prefix data/catalogues/tree_check/tree`.
  - Then check whether each candidate falls under the MRCA of the 10 references when the tree is rooted on
    GCF_000601485.1. The tree has 26 tips: 10 references, 8 outgroups and 8 candidates; GCA_035762045.1 is not included.

## Sources

| Source (version, date) | Genomes | f__Polyangiaceae | Screened | Hits | Notes |
|---|---|---|---|---|---|
| **NCBI GenBank+RefSeq**, Datasets v2 REST, listed 2026-11-11 | 9,434 Myxococcota (taxid 2818505) + "Deltaproteobacteria bacterium" (2026735) assemblies | NCBI names Polyangia-lineage: P1 = 1,013 post-cutoff | **all 4,477 selected** (finished 16:26; 720 by the 12:04 deadline) | 12 assemblies, 10 genomes (see update) | 5,194 not in R232; 2,984 released after 2025-09-04 |
| **SPIRE v01** (bulk table, README 2023-09; plus website API taxonomy, 2026-11-11) | 1,158,553 MAGs | 462 | 163 (162 blank or non-R232 genus, 1 target genus) | 2 | 299 assigned to an R232 non-clade genus were not downloaded |
| **GEM** (Nayfach 2021; genome_metadata.tsv 2020-11-30; GTDB r89) | 52,515 | 20 (16 blank genus) | 49 (16 blank-genus Polyangiaceae and 33 other o__Polyangiales in r89 families now inside R232 Polyangiaceae, or with a blank family) | 0 | max frac 0.12 |
| **OMD / OMDB v2.0** (API taxonomy R226; data table 2025-05-23) | 274,282 | 44 | 0 | 0 | all 44 are in R232 non-clade genera, so none were downloaded; the metadata directory returns 403 |
| **SMAG** (Zenodo 8223844 v1.0; GTDB r202/r214) | 40,039 | 200 (0 target) | 15 via NCBI (see below) | 1 (S400bin.4 = GCA_035762045.1, counted under NCBI) | details below |
| **MiDAS / Singleton 2021** HQ MAGs (Suppl. Data 3, GTDB r89) | 1,083 | 10 | 0 | 0 | all 10 have GCAs that R232 places in non-clade genera (JAAFHV01, JADJKB01, Xiheicaenibacterium, JADJCA01, JADJRH01, Polyangium) |
| **AnMAGC** anammox catalogue (Water Res 2024, doi 10.1016/j.watres.2024.122356; GTDB r207/r214) | 1,768 | 21 | 4 (target or blank/non-R232 genus) | 1 | MAG.1207 is labelled r214 g__JANWXX01 |
| **Brown et al. 2023** activated-sludge MAGs (Zenodo 10028566) | 876 | 11 (10 blank genus, 1 Labilithrix) | 10 | 3 | |
| **GROWdb** (Borton 2024; GTDB r207) | 2,093 | 3 (all SZUA-265) | 0 | 0 | MAGs are in NCBI under PRJNA946291 |
| **MGnify** genome catalogues (FTP genomes-all_metadata.tsv, 2026-11-11; per-catalogue counts in `mgnify_counts.tsv`) | 518,808 in 19 catalogues | 99 | 3 | 1 (redundant) | details below |

### NCBI details
- **Listing.** `ncbi_list.py` produced `ncbi_myxococcota_assemblies.tsv` (9,434 assemblies; 4,240 in R232, 5,194 not).
- **Selection.** `ncbi_select.py` produced `ncbi_to_screen.tsv` (4,477 assemblies). RefSeq copies whose GenBank twin is
  listed were dropped. The selection has three priorities:

  | Priority | Released | NCBI name | Assemblies |
  |---|---|---|---|
  | P1 | after the cutoff | lineage inside Polyangia | 1,013 |
  | P2 | after the cutoff | uninformative (Myxococcota/Myxococcales/Deltaproteobacteria bacterium, or "uncultured …") | 1,405 |
  | P3 | on or before the cutoff, but not in R232 | same name rules | 2,059 |

  The 717 not-in-R232 assemblies whose NCBI name sits in another named Myxococcota lineage were not selected.
- **Screening order.** Screening ran in this order:
  1. the 15 SMAG blank-genus MAGs that have a GCA and are not in R232;
  2. P1;
  3. the 309 Polyangia-named P3 assemblies;
  4. P2;
  5. the rest of P3.
- **Coverage at the deadline (stopped 12:04 EST): 720 of the 4,477 selected assemblies screened.**
  - All 720 were scored; there were 0 download failures.
  - P1: 705 of 1,013 (70%). P2: 0 of 1,405. P3: 15 of 2,059, the SMAG-linked ones.
  - **Not screened:** 308 P1, the 294 other Polyangia-named P3, all of P2 and the rest of P3.
  - Hits: GCA_982529355.1 and GCA_057137765.1 (P1), and GCA_035762045.1 (SMAG-linked P3).
- **Resume command.** It skips everything already in `screen_ncbi.tsv`:
  `uv run --no-project --with pyrodigal python -I src/catalogues/driver.py data/catalogues/ncbi_all_ordered_v2.txt ncbi 60`
- **Storage.** Genomes are downloaded in batches of 60 and deleted after screening unless they are provisional hits
  (`clade_aai` ≥ 60 or frac ≥ 0.3). Kept genomes are in `data/catalogues/ncbi/fna`.

### SMAG details
- 117 of the 200 Polyangiaceae MAGs match an NCBI GCA.
- 15 of those, all blank-genus in SMAG, are in R232 in non-clade genera or as GCA_035762175.1 itself, the DASTPP01
  reference. They are covered by R232.
- 57 are not in R232 and are in the NCBI screen list. The 15 blank or non-R232-genus ones were screened first.
  - One is a hit: S400bin.4 = GCA_035762045.1 (DASTPP01, 85% AAI).
  - The other 42, which carry an r214 genus outside the clade, were not screened.
- 17 qualifying MAGs have no GCA (blank or non-R232 genus, for example SRR6266545rbin.4: 90% complete, wetland, USA).
  They were **not reachable**:
  - the SMAG website no longer resolves in DNS;
  - the FASTAs exist only inside 38–39 GB tar files on Zenodo.

### MGnify details
- **Genome counts.** The 19 catalogues parsed hold 518,808 genomes. Most host-associated catalogues contain no
  Polyangiaceae.
- **Polyangiaceae found:**

  | Catalogue | Polyangiaceae | Blank genus | Target genus |
  |---|---|---|---|
  | soil v1.0 | 81 | 2 | 1 |
  | marine-sediment v1.0 | 5 | 0 | 0 |
  | maize-rhizosphere v1.0 | 5 | 0 | 0 |
  | tomato-rhizosphere v1.0 | 5 | 0 | 0 |
  | marine v2.0 | 3 | 0 | 0 |

  All other Polyangiaceae in MGnify are in R232 non-clade genera.
- **Screened genomes.** The three genomes below were screened. The FTP holds only `.gff.gz` files, so contigs were
  extracted from the `##FASTA` section by `mgnify_fetch.sh`.

  | Genome | Genus | Result |
  |---|---|---|
  | MGYG000501200 | DASTPP01 | redundant: same BioSample SAMN37960669 as GCA_035762175.1, ANI 100% |
  | MGYG000495370 | blank | no; frac 0.011 |
  | MGYG000495393 | blank | no; frac 0; 100% AAI to panel genome GCA_963782055.1, so it is that genome |

- **Not parsed:** none. All 19 catalogues listed on the MGnify FTP were parsed, human-gut v2.0.2 (289,232 genomes) and mouse-gut v1.0 (112,951) included; neither has any Polyangiaceae.
- The MGnify API `taxon_lineage` filter is ignored server-side, so the FTP metadata tables were used.

## Not searched or not reachable
- **SMAG:** the 17 non-GCA MAGs described above.
- **MiDAS:** 2,650 MQ MAGs are only in a 4.56 GB tarball with no taxonomy table. MiDAS 5 is 16S only.
- **OMD:** the metadata directory returns HTTP 403. Taxonomy came from the site's table API instead.
- **Hot springs and wetlands:** no dedicated catalogue with a GTDB table was found. SPIRE includes JGI Yellowstone hot
  spring MAGs, and one hit comes from there.
- **Login-only sources** were not attempted under the no-accounts rule: IMG/M and JGI GOLD-only bins.
- **SPIRE MAGs with an R232 non-clade genus label (299)** were not downloaded. The same applies to labelled-genus
  Polyangiaceae in OMD, MiDAS, GROWdb and MGnify. This assumes GTDB-Tk would not have put a clade member into a
  non-clade genus.

## Negative results
- **GEM:** 0 of 49.
- **OMD:** 0 target genus.
- **GROWdb:** 0.
- **MiDAS HQ:** 0. All 10 are already in R232, in non-clade genera.
- **MGnify:** 0 new. The only target-genus genome is a duplicate of GCA_035762175.1.
- **AnMAGC:** 3 of the 4 screened genomes are not members.
- **Brown23:** 7 of 10 are not members.
- **SPIRE:** 161 of 163 are not members.
