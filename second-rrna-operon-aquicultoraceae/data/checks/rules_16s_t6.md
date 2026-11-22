<!-- Copied from the lab record. Paths that point to files in this repository were rewritten to repo-relative ones; any other path refers to the original working folder and that file is not included (see data/README.md). -->

# T6 rules: are the divergent 16S's RNA bridges compatible with the normal 23S? (Fomalhaut, frozen before any A. secundus base at a contact position is read)

**Why.** Deneb: the divergent 23S of *Ca.* Aquicultor secundus is absent from RNA (0 : 129). If the divergent 16S does
accumulate (T5), it can only work by docking on the normal 23S, the layout of engineered orthogonal ribosomes (orthogonal
16S + wild-type 23S). Engineered orthogonal 16S rRNAs keep the intersubunit bridges; a natural divergent 16S that changed
its bridge bases together with the divergent 23S partners would be co-adapted to a 23S that is not there. This test is
sequence-only and is worth reading whatever T5 says.

**Inputs.**
- 16S: `out/t5/sec16.fa` (GCA_009993245.1 div = JAACTV010000069, norm = JAACTV010000080; full-CM spans, T5-A1).
- 23S: Deneb's `frozen-core/results/lsu/bac/cca/expr/loci.fa`, loci `GCA_009993245.1|JAACTV010000069.1|21177|23501|-|hmm8-2122`
  (div) and `GCA_009993245.1|JAACTV010000080.1|12238|14589|+|hmm3-2289` (norm), oriented to the E. coli sense strand by best
  alignment score to E. coli 23S.
- E. coli numbering: `data/contacts_16s/ecoli16S_num.tsv`, `data/contacts_16s/ecoli23S_num.tsv` (PDB 6OG7 chains, contiguous 2-1540 and 1-2903).
- Contacts: `data/contacts_16s/contacts_16s23s.tsv` from `src/orthogonal_ribosome/t6_contacts.py` (6OG7, non-rotated 70S, heavy atoms <= 4.0 A): 53 pairs,
  28 16S residues (set B16), 32 23S residues (set P23). RNA-RNA only; bridges made by 50S proteins are not included.

**Alignment.** mafft `--localpair --maxiterate 1000` on {E. coli, norm, div}, one run per rRNA. A residue is assessable if
E. coli, norm and div all have a base in its column (23S: the div locus is partial, so P23 residues outside it are not
assessable). Changed = div base != norm base.

**Statistics.**
- Bridge changes: number of assessable B16 residues changed. Background f16: changed share of assessable 16S columns
  outside B16; f23 likewise for 23S outside P23.
- Co-adaptation: over the contact pairs (i, j) with i a changed B16 residue and j assessable, the number with j also
  changed (div 23S != norm 23S); one-sided binomial against p = f23.

**Verdict (primary).**
- Gate: >= 24 of the 28 B16 residues assessable. Otherwise UNASSESSABLE.
- **COMPATIBLE:** <= 1 bridge change.
- **CO-ADAPTED TO THE DIVERGENT 23S:** >= 2 bridge changes, and the co-adaptation share >= 0.5 with binomial p <= 0.01.
- **CHANGED, PARTNERS NOT:** >= 2 bridge changes and the co-adaptation test fails. Compatibility then not settled by
  sequence.
- Reported, not in the verdict: per-residue table (E. coli number, E. coli/norm/div bases, partners and their bases);
  bridge changes vs f16 (binomial both tails; bridges sit in conserved regions, so depletion is expected anyway); the
  reverse view (changed P23 residues whose 16S partners are unchanged); the E. coli position of the div 16S insertions.

**Control PC-T6a (alignment and covariation logic), run and passed before the verdict is read.** Intra-16S base pairs
from `cmalign` (RF00177, `frozen-core/data/rfam/RF00177.cm`) SS_cons on norm and div; pairs where both sides are aligned
in both. Among pairs with at least one side changed, the share with both sides changed must be >= 3 x the share expected
if sides changed independently (f per side from the paired columns), one-sided binomial p <= 0.01; and >= 80 % of the
both-changed pairs must still pair in div (Watson-Crick or G-U). If PC-T6a fails, T6 is not read until explained.

**Caveats fixed now.** One E. coli state (non-rotated, RF2-bound). Bridge geometry can tolerate some base changes
(A-minor, ribose contacts), so a change is not proof of incompatibility and identity is not proof of assembly. Protein
bridges (uL2, uL14, bL19, uS13/uL5) not tested. One genome (the 8 Crystal Geyser genomes carry near-identical copies).

## T6 result note and T6-B (frozen 18:50, after T6 and an unfrozen look at A. secundus only; no other carrier's 16S read)
T6 verdict as frozen: CHANGED, PARTNERS NOT (data/contacts_16s/t6.txt). Unfrozen look (A. secundus div only, cmalign with E. coli):
the divergent 16S also changes universally conserved decoding/P-site residues (1054, 1492, 1493, 1400, 1498, 790, 966,
967, 1518). Question T6-B: is that rewritten state shared by the divergent 16S of the other carrier genera (maintained since
before the genera split), or particular to A. secundus (decay)?

- Units: the 5 other carriers of T1 (GCA_036504445.1 DASXVE01; GCA_902810715.1, GCA_027392615.1 CADDYS01;
  GCA_014859275.1 Kmv33; GCA_035643475.1 DASQQM01): 4 other genera. 16S = the locus span in out/t5/ref16.fa/class16.tsv.
  A copy is div if its blastn (-task blastn) top bitscore to A. secundus div > to A. secundus norm, else norm. A genome is
  used if it has >= 1 div and >= 1 norm (first of each by locus order). All copies cmaligned together (RF00177) with E. coli.
- R = A. secundus rewritten positions: the 14 changed bridge residues + changed named residues 1054, 1492, 1400, 1498, 790,
  966, 967, 1518 (22 positions). G = every other E. coli position where A. secundus div and norm both have a base and differ.
- s_R(genome) = share of R positions (other div has a base there) where the other div base = A. secundus div base; s_G likewise.
- Gate (alignment across genera): at R positions the other genomes' norm bases equal the A. secundus norm base in >= 90%.
- **SHARED REWRITE:** s_R >= 0.8 and s_R >= s_G in >= 3 of 4 other genera (CADDYS01: both genomes must pass).
- **A. SECUNDUS-SPECIFIC:** s_R <= 0.4 in >= 3 of 4 other genera. Otherwise **MIXED**.
- Reported: per-position table; s_R vs s_G (one-sided binomial, s_R > s_G); div-div identity over all aligned columns.
- **T6-B1 (~18:51, after classing copies by blastn bitscore, before any base was read).** (1) The non-Crystal-Geyser spans
  are HMM-only and stop at model 1465-1501, which would cut the 3' positions; spans are recomputed with T5-A1 `refine()`
  (full CM, hit +-400 nt). (2) Only GCA_027392615.1 and GCA_036504445.1 have an assembled normal 16S; the norm copy serves
  only the gate, so the gate uses those two normal copies, and every genome with a div copy is used for s_R.
- **T6-B result (~18:51):** gate 30/30; s_R 0.67-0.76 vs s_G 0.72-0.77 in every genus -> **MIXED** as frozen. Descriptive:
  the other divergent copies keep the E. coli state at only 5-14 % of R; 1492U, 1495G, 1498G, 1517A, 1518G, 966U, 967U,
  702U, 783G, 784C are in all 6 divergent copies (5 genera); 1493 is U in A. secundus and 3 others, C in DASQQM01, A (canonical) in DASXVE01 (data/contacts_16s/t6b.txt).

## T6-C (frozen ~18:51, before computing it): has the divergent 16S kept its structure since the carrier genera split?
Shared states alone can be relics of one ancestral copy decaying since. Between A. secundus div and each other genome's div
(T6-B alignment, out/t6/tmp/*/all.sto), over SS_cons pairs with all four bases present: pairs with >= 1 side changed,
both-changed share vs independence (f per side), and the share of changed pairs still WC/GU in the other copy.
Control: the same between A. secundus norm and the DASXVE01 and CADDYS01 (GCA_027392615.1) norm copies.
- **MAINTAINED:** in >= 3 of 4 genera, changed pairs still WC/GU >= 0.8 and both-changed share >= 3 x independence (p <= 0.01).
- **DECAYING:** changed pairs still WC/GU < 0.6 in >= 3 of 4 genera. Otherwise **UNSETTLED**.
- The norm control is read for comparison (alignment's own pairing bias); if the norm pairs-kept share is < 0.8 the
  verdict is UNSETTLED. CM alignment favours pairing, so MAINTAINED is read against the norm control, not alone.
- **T6-C result (18:52):** **MAINTAINED** (3/4 genera; DASXVE01 changed pairs kept 0.74). Div vs div: changed pairs still
  WC/GU 0.74-0.86, both-changed 5.1-7.1 x independence (p <= 3e-18); normal control 0.88/0.95. Reported null (not frozen
  beforehand): A. secundus div with 15 % random substitutions, same cmalign path (5 seeds): kept 0.24-0.40, 1.4-2.2 x
  (data/contacts_16s/t6c_neutral.txt), so the CM's own pairing bias does not produce the observed values.
- **Audit (19:02, fresh subagent):** no code or numbering bugs; contacts, mapping and counts reproduce. Wording
  corrected in RESULT.md: T6-B shows the paralog (not the specific rewrite) predates the split; T6-C passes by one pair
  (Kmv33 66/82) and div copies change ~5x faster than normal ones; neutral nulls post hoc (uniform and 55 % transitions,
  src/orthogonal_ribosome/t6c_neutral_gen.py); 783/784 bases differ between mafft and cmalign (both call them changed); H69 indel caveat.
