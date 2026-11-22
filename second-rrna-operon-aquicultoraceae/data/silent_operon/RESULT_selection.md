<!-- Copied from the lab record (work/silent-operon/, 22 Nov 2026). Paths such as out/, src/, recheck/ or frozen-core/ refer to the original working folder; the small tables from out/ are in this folder (see ../README.md), the rest is not included. -->
# T3 Selection along the molecule (B_origin, 22 Nov)

**Finding.** Selection keeps the divergent rRNA folded, but its functional core has been let go. Nucleotides that
touch tRNA/mRNA change 3-6x faster than the rest of the molecule, subunit-bridge nucleotides about 3x faster, and the
excess fades with distance from the tRNA sites. The fold is kept, and the gene has become more G-C-rich, fastest at
paired positions. The divergent 16S has no processing stem in any genus. This is neither a ribosome under weak constraint
nor a decaying pseudogene. It is a structured RNA kept for something other than translation.

Alignment: cmalign match columns (Rfam RF02541/RF00177), census conservation and E. coli numbering from
`frozen-core/results/{lsu,ssu}/bac/columns.tsv`; one divergent copy per population (6 for 23S and for 16S).
Contacts from E. coli 70S with mRNA and A/P tRNAs (PDB 7K00, heavy atoms <= 4 A).

Numbers:
- **Constraint on census-invariant sites is mostly gone, and stays gone after the genera split.** Between divergent copies
  of different genera (identity 0.82-0.87), invariant columns (f >= 0.99) change at 0.048 per site; functional 23S pairs at
  matched identity (0.80-0.89, n = 90) change them at 0.0014 (max 0.011). Share of changes falling at invariant columns
  relative to their share of columns: divergent 0.35 (23S), 0.47 (16S); matched functional 0.009 / 0.018; putative
  decaying Arsenophonus copy 0.67-0.69. "Conservation sensitivity" (variable-site rate / invariant-site rate): divergent
  5.2 (23S) and 3.4 (16S), functional ~250, decaying 1.7. By the pre-written rules this is neither "functional under
  relaxed constraint" (needs >= half the control's) nor "neutral decay" (< 2).
- **What is kept and what drifts** (deviation from census consensus at conserved columns f >= 0.95, observed / expected
  from each copy's own rate at the same conservation level and pairing status, 6 copies pooled):
  - tRNA/mRNA-contact nucleotides: 23S 94 / 17.3 (5.4x), 16S 54 / 17.2 (3.1x). Column level, 19 of 33 such 23S columns
    differ in some divergent copy, against 79 of 514 contact-free columns. 8 of the 33 carry the same change in all copies (ancestral).
  - 16S-23S bridge nucleotides: 23S 2.8x, 16S 3.7x.
  - r-protein-contact nucleotides: 23S 154 / 263 (0.59x), 16S 220 / 289 (0.76x); contact-free 1.03x / 1.06x.
  - Controls: the Arsenophonus decaying copies are ~1x in every class (ligand 1.1-1.5, protein 1.05); the Tenuifilum
    divergent operon keeps its ligand sites (0.76-0.84x) and does not keep protein sites (0.77-0.99x). So the pattern is
    specific to this operon.
  - Exploratory, not pre-registered: the best-kept protein footprints in 23S are those of uL24, bL21, bL36, bL17, bL20,
    uL14, bL28, uL4, uL13, bL34, uL15 (0-0.4x); the most drifted are bL27, bL32, uL22, uL16, bL33, bL35 (1.2-2x). In 16S, S15, S19, S2 and S13 are kept
    (0-0.3x); S3, S5, S10 and S8 drift (1.4-1.8x). Several kept ones are early-assembly binders, and the drifted ones
    are late binders near the PTC, tunnel and decoding site, but the column counts are small.
- **Named sites** (obs/exp, pooled): PTC 123/38.7, tunnel 47/9.2, H69 31/7.9, P loop 18/4.5, A loop 20/5.5, SRL 22/8.4,
  GTPase centre 46/23.1, decoding h44 24/8.2, h34 71/32, h28 neck 53/21 all drift *more* than the molecule (P < 1e-3).
  No functional site is kept intact beyond chance. The 530 loop (13/18.3) and the h23/1338 P-site (14/10.1) drift at the average rate.
- **Domain profile:** 23S domain III is fastest (1.7x the molecule rate vs 1.2 in functional pairs) and V slowest (0.79 vs 0.82); Spearman rho
  with the functional profile 0.60 over 6 domains. 16S 3' minor domain (h44/h45) is 2.4x the molecule rate vs 0.40 in functional pairs.
- **Fold kept and stems strengthened.** Fraction of census base pairs canonical: divergent 0.93-0.95 (23S) vs normal 0.95,
  Arsenophonus decaying 0.80. Among canonical pairs, the G-C share rose in every population: 23S 0.645-0.672 vs normal 0.608-0.613,
  16S 0.677-0.713 vs 0.641-0.645, while G-U fell. The decaying Arsenophonus copy goes the other way (0.39 vs 0.56). Where two
  divergent copies differ at both bases of a pair, 87% (23S) / 79% (16S) of the pairs stay canonical. These pair figures
  use a covariance-model alignment, which favours pairing; the kill agent's model-free value for div-vs-normal was 0.83-0.86.
- **Processing stems** (best complementary stem between 5' and 3' flanks vs 20 shuffles): every normal 16S in the family
  has one (z 3.9-9.5, score 45-56). 5 of 6 divergent 16S have none (z -1.0 to 1.2, score 22-33); CADDYS01 GCA_902810715.1
  has a weaker one (score 41) only with 450-nt windows. Divergent 23S keep a stem in 5 of 6 genera (score 37-51 vs 50-69 normal;
  none in Aquicultor_A). E. coli and Vibrio stems found (z 6.6-20.7); the test misses Tenuifilum even for its normal copies.
  This agrees with A_processing T1, which used different code on 2 genomes, and extends it to 6 genera.

**Mechanism.** The copy is under selection, but not for translation. The parts that grip tRNA, mRNA and the other subunit
have been let go, or changed early on (24% of the 23S contact columns carry a change shared by all genera). What is
kept is the overall fold and the periphery of the molecule, with paired positions gaining G-C. That the protein-binding
faces are kept on purpose is suggestive in the 16S only; in the 23S it is confounded with distance. The 16S seems to lack
the usual RNase III stem, though that signal is modest. This fits an rRNA that still folds, and so can still bind
ribosomal proteins, but cannot form a working ribosome. Two examples are a decoy or sponge that sets the free
r-protein level (r-protein translation is feedback-controlled by rRNA mimics) and a dead-end assembly intermediate.
These are hypotheses, not tested. It fits the RNA data: transcribed, the 16S seen only as precursor, the
23S not accumulating. It does not fit a second working ribosome.

**Recheck** (fresh subagent, own alignments, contact table and scorers; scripts in `recheck/`):
- tRNA/mRNA excess **holds**: mafft 5.3x / 3.0x, cmalign posterior >= 0.9 6.9x / 3.8x. Counting each column once, 20 of 30
  contacted 23S columns deviate vs 4.3 expected (16S 13 vs 3.8), permutation p <= 0.0005. So the honest independent count
  is ~20 columns, not 94 events.
- **Caveat that changes the reading:** the effect is a gradient with distance from tRNA/mRNA (23S: 5.9x within 4 A,
  2.7x at 4-12 A, 1.3x at 12-20 A, 0.57x beyond 30 A). The 23S protein-contact deficit mostly disappears once distance is
  controlled (beyond 12 A: protein 0.66x vs uncontacted 0.69x). It survives in the 16S (0.74x vs 0.95x). So "r-protein sites kept"
  is **weakened** to: the whole functional core drifts, and the periphery, where most protein contacts are, is kept.
  The per-protein list above inherits this confound.
- GC **holds as a measurement** under mafft too (23S 0.646-0.676 vs 0.607-0.618). It is **not stem-specific**: loops gained
  GC as well. Per site, AT->GC changes run ~2.5x faster than GC->AT at paired columns vs ~0.6-0.7x at unpaired ones,
  which still points to selection on stems. Two of the six populations lack a full own normal copy.
- Processing stems **weakened**: divergent 16S have no stem in all 6 and in JAMCWL01 (z -1.4 to 1.0, any window or 3'
  shift). But normal 16S stems are weak with its scorer (z 0.9-4.8; A. secundus ~1), and the Tenuifilum divergent copy also
  lacks one. Its scorer finds a 23S stem in Aquicultor_A as well (z 2.8), so the divergent 23S keep one in 6 of 6.

Files: `src/t3_profile.py`, `src/t3_summary.py`, `src/t3_sites.py`, `src/t3_contacts.py`, `src/t3_classes.py`,
`src/t3_gc.py`, `src/t3_stems.py`; outputs `out/t3_*`, `out/contacts_7k00.tsv`.
