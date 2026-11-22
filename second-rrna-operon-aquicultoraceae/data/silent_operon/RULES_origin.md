<!-- Copied from the lab record (work/silent-operon/, 22 Nov 2026). Paths such as out/, src/, recheck/ or frozen-core/ refer to the original working folder; the small tables from out/ are in this folder (see ../README.md), the rest is not included. -->
# B_origin rules (written 22 Nov before any result was read)

Question: where did the divergent operon come from, and is something selecting to keep it?
Inputs reused: the 65 order-level genomes and their Prodigal proteins (`frozen-core/results/lsu/bac/cca/aq/fna`,
`kill_div/faa`), the div/norm labels of `kill_div/loci_blast.tsv`, census per-column conservation
(`frozen-core/results/lsu/bac/columns.tsv`, `results/ssu/bac/columns.tsv`; Rfam RF02541 / RF00177 match columns with
E. coli numbering), GTDB R232 tree and SSU files. Alignments: cmalign to the same Rfam models (match columns only).

Controls everywhere: (i) the normal operon of the same genomes; (ii) E. coli K-12 (7 operons, concerted evolution);
(iii) known divergent working operons: Vibrio vulnificus MO6-24/O (rrnI, the specialised-ribosome operon of Song et
al. 2019) and Thermobispora bispora DSM 43833 (two 16S types ~6% apart, both expressed, Wang et al. 1997).

## T1 Origin (gene trees, 23S and 16S separately)
Set: every div and norm copy >= 60% of model length from the 65 genomes; the 10 best non-family hits of each
divergent copy in the census 23S (all copies, `scan23/census_lsu_all`) and GTDB `ssu_all_r232`; one normal copy per
other family of the order; outgroups from other classes. Tree: IQ-TREE 3, GTR+F+R4 (or ModelFinder), UFBoot 1000, 1 thread.
- **Duplication inside the lineage**: divergent copies are monophyletic and their sister (UFBoot >= 90) is a group of
  normal copies from this order (Aquicultoraceae alone = duplication in the family ancestor; several families =
  older duplication with losses).
- **Horizontal transfer**: the divergent clade is nested (UFBoot >= 90) inside copies of another class/phylum, OR a
  non-family sequence is more similar to the divergent copy than the genome's own normal copy is (identity over
  >= 1.2 kb, at least 3 points higher). Then that lineage is the donor candidate.
- **Phage / mobile origin**: any phage or plasmid sequence among the closest hits, or mobile genes in T2 (below).
- **Unresolved**: divergent clade on a long branch whose attachment point has UFBoot < 90, or differs between the 23S
  and 16S trees. Then I say "no donor detectable; consistent with an old duplication but not proven".
- **Vertical inheritance since one ancestor**: the topology among genera of the divergent copies has no conflict
  with UFBoot >= 95 against the normal-copy tree or the GTDB species tree of the same genomes. One or more such
  conflicts = transfer of the operon between genera.
- Control expectation: E. coli's 7 copies form one tight clade; V. vulnificus rrnI and T. bispora's two types sit
  inside their own genus (duplication or concerted evolution, not foreign). If a control lands outside its genus,
  the tree method is not trusted.
- Age: the divergent-vs-normal split is placed relative to the family / order node of the normal copies.

## T2 Neighbourhood (+-10 kb of each 16S-23S operon)
Proteins: Prodigal calls already made. Annotation: DIAMOND vs Swiss-Prot, hmmsearch vs Pfam-A (E <= 1e-5), DIAMOND
vs the lab phage protein set (`new-stops-prok/data/phage/phageref.dmnd`, e <= 1e-5, >= 30% id). tRNAs: aragorn.
Orthogroups: DIAMOND all-v-all among flank proteins, linked at e <= 1e-5 and >= 50% query cover, single linkage.
- **Conserved with the operon**: an orthogroup present within 10 kb of the divergent operon in >= 3 genera.
- **Specific to the divergent operon**: conserved with it AND present at <= 1 normal-operon neighbourhood.
- **Suggestive class**: ribosomal-protein paralog, rRNA/tRNA methyltransferase or other rRNA modifier, RNase,
  toxin-antitoxin (Pfam PIN, RelE/ParE, MazF/PemK, HicA/HicB, HigB, VapB/C, Doc, GNAT-toxin, AbiEii, etc.),
  phage/defence (Cas, RM, Abi, Gabija, retron...), transposase/integrase/recombinase.
- **Island / mobile origin supported** if: a mobile gene (integrase, transposase, phage) is conserved with the
  divergent operon in >= 3 genera, OR the +-10 kb region's GC differs from the genome median of 20-kb windows by
  > 2 SD in >= half the carriers, OR a tRNA gene sits at the island edge with a direct repeat in >= 3 genera.
- **Candidate partner** if a suggestive-class gene is specific to the divergent operon in >= 3 genera.
- Otherwise: "fixed ordinary chromosomal locus, no partner gene seen".
- Controls: the same scoring on the normal operon flanks; E. coli 7 rrn flanks (no mobile genes expected);
  base rate of each class among all genes of the carriers.

## T3 Selection along the molecule
Alignment: cmalign match columns (RF02541 / RF00177) of divergent copies (one per genus) and normal copies, mapped
to E. coli numbering through `columns.tsv`. Regions: 23S domains 0, I-VI and named sites (PTC ring, P loop H80,
A loop H92, H89, SRL 2653-2667, H69 1906-1924, GTPase centre 1051-1108, L1 stalk H76-78, H68, exit-tunnel
uL4/uL22 contacts, 5'/3' end stem H1); 16S 5', central, 3' major, 3' minor domains and h18, h44, 530 loop,
790 loop, h45, anti-SD.
Measures:
1. **Relative rate per region** = substitutions per site in the region / molecule-wide, between divergent copies of
   different genera. Same for normal copies of the same genus pairs, and for functional pairs at matched overall
   divergence (normal copies from different families of the order, and the control pairs from the census chosen to
   match the divergent pairs' identity within 2 points).
2. **Conservation sensitivity** = substitution rate at census-variable columns (f_cons < 0.9) divided by rate at
   census-invariant columns (f_cons >= 0.99), for divergent-vs-own-normal and for matched functional pairs. A neutral
   pseudogene gives ~1-2 (only mutational bias); functional rRNA gives >> 5.
3. **Site intactness**: each named site identical to the normal copy in all divergent copies, or not.
   Expected identical by chance = (per-site identity in the same census-conservation bin)^(sites x copies).
Readings:
- **Functional ribosome under relaxed constraint**: region profile correlates with the matched functional profile
  (Spearman rho >= 0.6) and conservation sensitivity >= half the control's.
- **Neutral decay**: sensitivity < 2 and flat profile (rho < 0.3).
- **Selection for a different, narrow job**: a few regions kept (relative rate <= 0.5 of the divergent average and
  intact beyond chance, P < 0.01) while functional centres (PTC, decoding) drift like the average. The kept regions
  name the candidate job (e.g. SRL/GTPase centre kept = factor binding; processing stems kept = processing; only
  protein-binding helices kept = protein sponge).
- Processing sites: the leader/trailer stems around 16S and 23S (RNase III) — do the bases flanking each divergent
  gene pair into a long stem as in the normal operon (>= 20 bp with <= 3 mismatches)? Yes/no per copy.

Any result I call a lead gets a fresh-subagent recheck with raw inputs before it goes into RESULT.md.
