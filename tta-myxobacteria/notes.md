# Notes: bacteria that end their genes at TTA

The agents' full write-up (method, evidence, checks, open points). The short version is in [`README.md`](README.md).

**Summary.** We screened all 199,923 GTDB R232 species-representative genomes for codons that pile up where related
proteins end. In one clade of Polyangiaceae (phylum Myxococcota; 5 GTDB genera, 8 species), protein-coding genes end
at TTA, exactly where related genes end at an ordinary stop. TTA is read as leucine in every bacterial genetic code
known. Three GTDB genomes pass every check fixed before scoring, including raw reads and a blinded re-check by a
different method. Reads carry the TTA at every testable site in all 12 clade genomes checked (285 of 285). No clade
genome has a detectable tRNA-Leu(UAA), and RF2 carries a pair of changes found together in no other GTDB
representative. One member, 1-FT3.2, is kept in co-culture with its prey. **All evidence is computational. The
genomes are metagenome-assembled except 1-FT3.2, which is assembled from a co-culture. No protein data exist.**

## What was known
- **Sense codons that became stops** are known in organelles, plus one 2026 proposal in archaea (KACI, below).
  - Mitochondria: UUA in all labyrinthuleans, *Thraustochytrium* among them (NCBI table 23; Žihala, Salamonová & Eliáš
    2020); UCA in *Scenedesmus* and other Sphaeropleales; AGA/AGG in vertebrates (read by mtRF1; Krüger et al. 2023,
    Saurer et al. 2023). UUA/UUG were proposed as stops in *Pycnococcus* mitochondria (Turmel, Otis & Lemieux 2010) and
    as putative stops in Marsupiomonadales mitochondria, relying on changes in mtRF1a (Barcytė, Žihala & Eliáš 2025).
  - Plastids: UUA and UCA in the green dinoflagellate *Lepidodinium* and relatives (Barcytė, Žihala & Eliáš 2025,
    *PLoS Genet* 21:e1011901), "the first case of a sense-to-stop codon reassignment ever discovered outside
    mitochondria". Their bacterial-type plastid RF2 has a mutated Ser204 (as in labyrinthulean mtRF2a, proposed to
    mediate UUA recognition) and further changes. This is the closest precedent for an RF2 route. We have not mapped
    their residue numbering onto ours.
- **Archaea.** KACI (Melnykov 2026, *PNAS* 123:e2610659123; preprint bioRxiv 10.64898/2026.04.02.715157) applies a fast
  code-inference method to thousands of assemblies. Its significance statement reports "the first case outside of
  eukaryotic organelles of a sense codon reassigned to serve as translation stop". In the preprint this is CGA in two
  Njordarchaeales (Asgard) genomes, raised because CGA never appears in their protein-domain matches and
  tRNA-Arg(UCG) is missing; the preprint adds that its approach "does not allow direct inference of stop codons". It
  does not look at gene ends. (The *PNAS* body was not readable to us; PMC embargo.) On our screen, one of those two
  genomes (GCA_027068385.1, a GTDB representative) does show genes ending at CGA, but below our bar: 20 families,
  4.4× expected, P 8e-8 against a bar of 1e-10.
- **Bacterial code changes.** Stop to sense (UGA→Trp in *Mycoplasma*, symbionts and Eggerthellaceae; UGA→Gly in
  Gracilibacteria; TAG in [`../tag-bacteria`](../tag-bacteria)) and sense to
  sense (five arginine-codon changes, Shulgina & Eddy 2021, *eLife* 10:e71402; ACA Thr→Asp, Fakih et al. 2026,
  bioRxiv 10.64898/2026.05.03.722492). None is sense to stop.
- **Unassigned codons in bacteria.** *Mycoplasma capricolum* CGG has no tRNA and stalls ribosomes without release in
  a cell-free system (Oba et al. 1991). *Micrococcus luteus* AGA/AUA were proposed unassigned (Kano et al. 1993); the
  later genome uses both and encodes tRNAs that may read them.
- **What existing screens can see.** Codetta outputs a stop as "?", the same symbol as a codon with too little data.
  gTranslate (GTDB, 2026) chooses only among tables 11, 4 and 25. KACI is described above.
- **Searches** (web, Europe PMC including bioRxiv, to 2026-11-11) found nothing on the genetic code of
  Polyangiaceae/Myxococcota. The genome note on 1-FT3.2 (Demina et al., *F1000Research* 14:1153) does not discuss its
  genetic code.
- **Mechanism precedent.** A mutant *E. coli* RF2 with two charge changes (E149K/E170K) releases at UUA in vitro,
  though also at UGG, UAC, UCA and CAA (Uno, Ito & Nakamura 2002, *PNAS* 99:1819).

## Method
1. **Arm A (screen).** Codon use near gene ends in all 199,923 GTDB R232 representatives (`src/gtdb_scan.py`,
   `src/armA.py`, `src/armA_run.py`). 2,052 genomes were flagged (`data/screen/arm_a_flags.tsv.gz`).
2. **Arm B (homology).** Each protein is read on to the next TAA/TAG/TGA and aligned to Swiss-Prot (DIAMOND); we record
   the codon just past the point where the homolog's C-terminus ends (`src/term_codon.py`). Pass rule, frozen before
   scoring: at least 10 families, at least 3× expected, P ≤ 1e-10. Run on 4,599 genomes (`src/armB_select.py`);
   31 (genome, codon) candidates passed (`data/screen/arm_b_candidates.tsv`; this table includes the non-TTA
   candidates, which are not claimed here).
3. **Four frozen checks** for a candidate to become a result (`data/prereg.md`, "From an Arm B pass to a result"):
   1. another assembly of the genus passes;
   2. at least 3× stronger than the codon's synonyms;
   3. raw reads carry the codon at 10 or more event sites (depth ≥ 5, ≥ 90% of reads). Each read pair counts once;
   4. a blinded subagent re-check by tblastn of conserved proteins on raw contigs, with no gene calls. A stricter rule
      (observed/expected ≥ 3 downstream) was added when check 4 was unblinded; it can only remove passes.
4. **After the checks** (post hoc; none changes a verdict): an adversarial review on close-relative orthologs with no
   Swiss-Prot; a blinded 16-genome panel scored by a separate agent with a method frozen before scoring; tRNA genes
   (ARAGORN, tRNAscan-SE 2), release factors on raw contigs, and fixed differences across 205 shared genes; all of
   GTDB R232, not only representatives; a search outside GTDB.

## Results
**The frozen chain** (all TTA, all Polyangiaceae):

| genome | genus | Arm B families (ratio, log10 P) | 1 genus | 2 synonym | 3 reads | 4 re-check | status |
|---|---|---|---|---|---|---|---|
| GCA_963169165.1 | JADLFU01 | 15 (30×, −16.8) | yes | 62× | 10/10 | yes | result |
| GCA_964605685.1 | JANWXX01 | 19 (19×, −18.4) | yes | 30× | 20/20 | yes | result |
| GCA_020849655.1 | JADLFU01 | 19 (38×, −23.0) | yes | 70× | 18/18 | yes | result |
| GCA_025453595.1 | JANWXX01 | 32 (50×, −41.8) | yes | 76× | no read runs | yes | unconfirmed |
| GCA_025057345.1 | JANWXX01 | 21 (16×, −18.0) | yes | 20× | 15/15 | no (3 vs 0.5) | unconfirmed |
| GCA_035762175.1 | DASTPP01 | 16 (11×, −11.4) | no other assembly | 17× | 16/16 | no (4 vs 0.5) | unconfirmed |
| GCA_046351315.1 | JAXLYQ01 | 14 (13×, −10.8) | no other assembly | 18× | 11/11 | no (4 vs 0.5) | unconfirmed |
| GCA_964396175.1 | CAYNEK01 | 5 (10×, −3.8): below the bar | | | | | not a candidate |
| GCA_963169955.1 | JADLFU01 | 20 (39×, −24.5) | yes | 78× | 20/20 | not run | passes 1–3 |
| GCA_033578855.1 | JANWXX01 | 20 (14×, −16.0) | yes | 17× | 14/14 | not run | passes 1–3 |

- Check 4 looks at only about 100 conserved markers per genome, so its expected count sits at a floor of 0.5; it
  was weak for single genomes.
- The last two rows are the only non-representative GTDB genomes of the clade's genera. GCA_963169955.1 comes from
  the same sample and read run as GCA_963169165.1 (a different species; different loci); GCA_033578855.1 is a hybrid
  reassembly of the samples behind GCA_025057345.1. Neither is independent evidence.
- **One clade.** On the GTDB tree the 7 TTA candidates and CAYNEK01 form an 8-tip clade. Arm A flagged exactly 7
  Polyangiaceae genomes, all for TTA, all in this clade.

**Close-relative ortholog ends** (adversarial review; `data/ortholog_ends/`):
- In each clade genome, take the genes that run past a sister ortholog's stop. The next codon is TTA in 21/260,
  36/245, 63/328, 24/524, 32/208, 24/269, 6/196 and 46/309 of them: 252 pooled against about 9.6 expected. In the 8
  sisters it is 0 in 7 and 1 of 195 in the eighth (1 of 2,025).
- By observed/expected, TTA ranks first of the 61 sense codons in every clade genome.
- The sequence after a TTA end evolves like non-coding DNA; the bases around TTA ends look like stop context, not
  leucine context; the signal is not contamination and is spread over about 100 gene families.

**Blinded re-check** (`data/blinded_panel/`): a separate agent scored the 8 clade genomes and their 8 nearest
relatives, shuffled and renamed, for all 61 sense codons (blastx of contigs against Swiss-Prot, no gene calls). Its 8
strongest UUA scores were exactly the 8 clade genomes (1/12,870 by chance). It called 6 outright; it described
JAXLYQ01 as "strong" (not one of its frozen tiers) and CAYNEK01 as "suggestive". All 8 controls showed nothing. No other
codon reached even a suggestive level in any genome.

**TTA trades places with the stops** (`data/shared_ends/`): place a gene end of one clade genome against a sister
ortholog, then read the codon at the same place in the ortholog of each other clade genome. Across genera, TTA ends
are TTA again at 77 of 641 resolved ends (12.0%), a canonical stop at 78.2%, a sense codon at 9.8%. Canonical stops
turn into TTA at 3.9–4.6%. About 90% of TTA ends and of canonical ends are an end in the other genome too. In the
sisters, TTA stands at 3 of 38,417 canonical ends, against 777 of 18,007 in the clade. 43 of 183 TTA-end genes with a
resolved ortholog in another genus share their TTA with it, which fits TTA ends being older than the genera.

**Machinery** (post hoc; `data/machinery/`):
- **tRNA.** tRNA-Leu(UAA) is not detected in any of the 8 clade genomes (59–97% complete). All 8 sisters carry one
  copy (one-sided Fisher P = 7.8e-5). ARAGORN and tRNAscan-SE 2 agree genome by genome; *E. coli* K-12 is the positive
  control. The gap is not completeness: the clade genomes have 30–92 tRNAs (tRNAscan-SE), and the two smallest sister
  sets (20 and 23) still include Leu(UAA). Not specific on its own: 18 Polyangiaceae representatives in a partial
  tRNA table lack it, and only 4 of them are TTA-flagged.
- **RF2.** The recognition motif is SPY in 7 of the 8 clade genomes and SPF in every sister. SPY alone is common (9.7%
  of the 161,077 GTDB representatives with an RF2 call) and not enriched for TTA flags. A second change four residues
  before the motif, *E. coli* L201F, is in all 8 clade genomes. F201 together with SPY occurs in exactly 7 GTDB
  representatives: the 7 SPY clade genomes. Of 205 genes present in all 16 genomes with ≥ 100 aligned sites, RF2 has
  the most clade-derived fixed differences (10 in 262 sites; the next two genes, 8 each; no significance test).
- **RF1** (TAA/TAG) is unchanged (PAT).

**Outside GTDB** (`data/outside_gtdb/`): a competitive marker test (437 clade core proteins, rule fixed from controls
first) over NCBI Myxococcota assemblies not in R232 (4,477), SPIRE, GEM, OMD, SMAG, MiDAS, AnMAGC, Brown et al.
activated-sludge MAGs, GROWdb and all 19 MGnify catalogues. The frozen Arm B rule and check 2 were then run unchanged.
Seven genomes pass: 1-FT3.2 (GCA_982529355.1; 102 families, 25×, log10 P −103.7; reads 104/104), three Brown et al.
activated-sludge bins (no source reads identified), MAG.1207 (anammox; reads 28/28), spire_mag_00879468 (wastewater;
13/13) and GCA_060838245.1 (wastewater plant, Virginia; 16/16). Nine more rank TTA first of 61 codons but fall below
the bar. Two marker trees put all catalogue members inside the clade (SH-aLRT 100); where three of them sit inside it
is not settled.

**1-FT3.2.** A predator isolated from subarctic peatland soil at Pallas, Finland, with *Mucilaginibacter cryoferens*
FT3.2 as prey (University of Helsinki and Natural Resources Institute Finland; Demina, Ihonen, Männistö & Hultman,
*F1000Research* 14:1153). Kept as a co-culture; it grows only slowly and weakly without prey. Genome assembled from
co-culture reads after removing the prey's (CheckM2 81.9% complete, 0.38% contamination). Marker AAI 68–70% to
JANWXX01. RF2 F201 + SPY; no tRNA with a UAA anticodon.

## What could still be wrong
- **Mechanism unknown.** A release-factor stop and an unassigned codon leave the same marks in a genome. The RF2
  changes fit the first; the tRNA loss fits both.
- **JAXLYQ01 keeps SPF.** Its only RF2 found is SPF (genome 59% complete), but it shares L201F. So SPY may be one of
  several RF2 changes, or not the cause.
- **Some TTAs sit inside genes.** 207 interior TTA sites lie in aligned regions of clade genes. In conserved genes
  they come at 1.0–2.8× the rate a canonical stop would show; none of these tests is significant once off-GC contigs
  are removed (P 0.10–0.56). So TTA may not be a strict stop every time.
- **No pure isolate, no proteins.** All but 1-FT3.2 are metagenome-assembled; 1-FT3.2 is a co-culture. The three
  full-chain genomes are 84.6%, 85.6% and 96.9% complete.
- **One genome without reads.** GCA_025453595.1 (wetland) has no public read runs.
- **Possibly not the only case.** 14 Streptosporangiaceae TTA candidates are unconfirmed, not refuted: check 4 had
  almost no power for them, and reads support their TTA in the four checked.
- **Post hoc.** The machinery work, the adversarial review, the blinded panel, the shared-ends test, the
  non-representative members and the catalogue genomes all came after the frozen chain. They test the same claim in
  new ways but were not in the plan.

## What would settle it
C-terminal proteomics would show whether proteins end at TTA or carry leucine there. 1-FT3.2 makes this direct: it
can be grown in co-culture with its prey. The anammox reactor biofilms that carry JADLFU01 are another option. An in
vitro release assay with the clade's RF2 on a UUA codon would test the release-factor route.

## References
- Barcytė D, Žihala D, Eliáš M (2025) Expanded diversity of pedinophytes provides a window into the evolution of the
  genetic code. *PLoS Genet* 21:e1011901. doi:10.1371/journal.pgen.1011901
- Demina T, Ihonen R, Männistö MK, Hultman J (2025) Draft genome sequence of a predatory bacterium from northern
  peatland soil. *F1000Research* 14:1153. doi:10.12688/f1000research.171863.2
- Kano A, Ohama T, Abe R, Osawa S (1993) Unassigned or nonsense codons in *Micrococcus luteus*. *J Mol Biol* 230:51.
- Melnykov AV (2026) Alternative genetic codes in bacteria and archaea identified with a fast k-mer-based algorithm.
  *PNAS* 123:e2610659123. doi:10.1073/pnas.2610659123
- Oba T, Andachi Y, Muto A, Osawa S (1991) CGG: an unassigned or nonsense codon in *Mycoplasma capricolum*. *PNAS*
  88:921.
- Shulgina Y, Eddy SR (2021) A computational screen for alternative genetic codes in over 250,000 genomes reveals
  five new codon reassignments. *eLife* 10:e71402.
- Uno M, Ito K, Nakamura Y (2002) Polypeptide release at sense and noncognate stop codons by localized charge-exchange
  alterations in translational release factors. *PNAS* 99:1819.
- Žihala D, Salamonová J, Eliáš M (2020) Evolution of the genetic code in the mitochondria of Labyrinthulea
  (Stramenopiles). *Mol Phylogenet Evol* 152:106908.
