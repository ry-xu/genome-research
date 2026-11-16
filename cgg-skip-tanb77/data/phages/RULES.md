# RULES — code-walls v1 (TANB77 / CGG first)

Draft 15 Nov 2026 15:25 (the phage agent); redrafted 17:21 after P4 and before any other control result was read. Frozen when the "FROZEN" line below carries a time and sha1; nothing about the
TANB77-host phages (the target group) is computed before that. Controls may be run before freezing only to check
that code runs; any change made after seeing a control is logged under Amendments with the reason.

## Inputs
- UHGV v1.0 (Zenodo 10.5281/zenodo.17402089): `votus_full.fna.zst` (one representative per vOTU) and
  `votus_metadata_extended.tsv.zst`. Host = `crispr_host_lineage_gtdb_r207` only (spacer evidence); k-mer host calls
  are not used for group membership.
- Gene calls: pyrodigal, metagenomic mode, translation table 11, every CDS >= 300 nt that is not partial at a contig
  end. Only phages with `checkv_quality` in {Complete, High-quality, Medium-quality} and UHGV `genetic_code` == 11
  (recoded phages, code 4/15, are reported separately, not pooled).

## Groups (sets in `data/groups/run_ids.tsv`, built by `src/groups.py` + `src/pick.py`, seed 20261115)
- **T1 (primary target):** CRISPR host genus in o__TANB77 whose genome calls the CGG skip in frameshift-codes
  (`scratch/cgg_test_T.tsv`, call = 1): 436 phages.
- **T0 (within-order control):** host genus in TANB77 that does not call the skip (CAG-492, CAG-465): 137 phages.
  Same order, habitat and GC; the code-wall hypothesis predicts T1 differs from T0 as well as from R.
- **Tf (secondary):** host given only to family f__CAG-508 (46/50 genera skip): 664 phages, reported separately.
  Tx (69, other TANB77 host ranks) is not tested.
- **R (reference):** CRISPR host class c__Clostridia, any order except TANB77, at most 500 per 2-point GC bin
  (seeded random), only bins that hold T phages: 6,307 phages.
- Host-range: `host_range_breadth.tsv.zst`, crispr_spacer rows; a phage "crosses orders" if hrb_rank is order,
  class, phylum or domain.

## Statistic (per phage, per disputed codon; here CGG)
- O = in-frame CGG count over the CDS set.
- E = n(CGN in frame) x g3, where g3 = share of G at third positions of the seven non-Arg four-fold families (GGN,
  GCN, CCN, ACN, GTN, CTN, TCN) in the same CDS set.
- Skip signature, homology-free: an in-frame CGG counts as a skip site if it lies in the last 16 codons of the CDS
  and the +1 frame read from the codon's second base is open for >= 50 codons. P4 (below) showed that host genomes
  that skip have skip/O 0.33 vs 0.005 in standard low-GC Clostridia, and raw r ~1.3 (truncated gene calls end at
  skip sites), so host-style use shows as high skip/O, not as low O/E.
- Per-phage E is ~1-3 at GC 26-30%, so the main tests are pooled; per-phage classes are descriptive only.

## Main tests (`src/score.py groups`)
For each comparison (T1 vs R; T0 vs R; T1 vs T0 with T0 as reference; Tf vs R):
- **A. Usage:** Mantel-Haenszel rate ratio of O over E (E as exposure), strata = 2-point GC bins: RR_use.
- **B. Skip signature:** MH ratio of skip sites over O, same strata: RR_skip.
- p = two-sided within-bin label permutation (2000 permutations); 95% CI = bootstrap over UHGV phage genera
  (vGENUS in `uhgv_taxonomy`, else the vOTU itself; `genus_name` is "NULL" throughout), 1000 draws, target and reference resampled separately.
- **Readings for T1 vs R, fixed now:**
  - *avoid* if RR_use < 0.5, p < 0.01 and CI upper < 0.75, and specific: RR_use(CGG) < 0.67 x the same MH rate
    ratio for each sibling CGA, CGT, CGC (E = n(CGN) x share of A/T/C at four-fold non-Arg thirds). Reason: Taylor
    et al. 2013 (PeerJ 1:e50) found a CUG-yeast virus avoiding the neighbour CUA as well (prior-art note);
  - *host-style* if RR_skip > 3, p < 0.01 and CI lower > 1.5;
  - *ancestral (no wall in codon use)* if the RR_use CI lies inside 0.75-1.33 and the RR_skip CI upper < 1.5;
  - anything else is *unresolved* and reported as such. The permutation p treats phages as independent; T phages
    cluster in vGENUS, which is why every reading also needs the cluster-bootstrap CI condition. The same readings are applied to T0 vs R and Tf vs R;
    a code-wall effect requires T1 to differ from T0 in the same direction (T1 vs T0 p < 0.05).
- **C. Firewall half:** share of phages that cross orders in T1+Tf vs R (Fisher exact; GC bins reported), and the
  per-phage class of every crossing T phage.
- Descriptive per-phage classes (avoid / host-style / ancestral / unresolved, as in the draft of 15:25) and phage
  tRNAs with anticodon CCG or UCG (pyaragorn).

## Controls (must pass before the target counts)
- **P1, known avoidance (lens case; Nibert, Vong, Fugate & Debat 2018, Virology 518:14-24).** Weak control: that
  paper used "no internal UGA" as part of its case for plant hosts, so plant-set avoidance is partly circular
  (priorart-1.md); P2 and P3 carry the method check. mitoviruses, disputed codon TGA. O = in-frame TGA (ORF read with
  table 4, longest ORF), E = n(Trp codons TGA+TGG) x a3, a3 = A share at third positions of NAR two-fold families
  (AAR, GAR, CAR) . Pass if plant-hosted mitoviruses with E >= 5 have median r < 0.25 and fungal-hosted median r >
  0.5, and the two groups separate with Mann-Whitney p < 0.001. List built blind from host attribution
  (`data/mito_control.tsv`).
- **P2, known host-following (prokaryotes):** the 84 UHGV code-4 phages (genetic_code 4, any host) under the TGA
  statistic, gene calls pyrodigal single mode table 4 (meta mode for phages < 20 kb), pooled O/E > 0.5 and median
  per-phage r > 0.5 among phages with E >= 5; R phages (table 11) give pooled TGA O/E < 0.05.
- **P3, recall (group level):** fake targets with T1's GC-bin sizes drawn from R, their stats taken with a share f of
  in-frame CGG replaced by another CGN (from the phage's own CGN use); reference = the rest of R, unplanted.
  f = 0.5 must give p < 0.01 in >= 80% of 200 fake targets and median RR_use 0.4-0.6; f = 1.0 reported.
- **N, specificity:** the null world below (fake targets from R, unplanted): share of p < 0.05 <= 0.10 and share of
  p < 0.01 <= 0.03 for both A and B; per-phage avoid calls in a held-out half of R <= 2%.
- **P4, skip-signature calibration:** TANB77 host genomes that skip (frameshift-codes `data/tanb77_set.tsv`, calls
  from `results/tanb77`) must show s well above low-GC standard Clostridia genomes
  (`data/clos_lowgc_set.tsv`): median s in skip genomes >= 5x the median in standard genomes. If not, the host-style
  class is dropped and only avoid/ancestral are reported.
- If P1 or P3 fails, fix the method (logged) or stop.

## Kill rules for any lead
- Every avoid/host-style T phage reported as a lead is re-read: contigs with gene density < 0.7 or with UHGV
  `genomad_plasmid_hallmarks` > 0 flagged; prophages carrying host chromosome flanks trimmed (checkv_trimmed).
- A claim about the order needs >= 3 TANB77 host genera: T1 is also scored leaving out each host genus in turn
  (the reading must hold in every leave-one-genus-out run).
- A fresh subagent re-derives the main group test from raw files before any post calls it a result.

## FROZEN
Frozen 15 Nov 2026 17:50 (the phage agent), before any statistic of a T1/T0/Tf/Tx phage was computed.
sha1: 24bbbfacd660f855b8cef7c459c0afe21cc2996b src/codon_stats.py; b0414714607f58caceaa8cf59603fc9323a539db src/score.py; f50697fb3efa36829132f3c7d7a4f8e7d8aa37ed src/groups.py; 8bfcfda0f89cfa5fd2c40e2c8c45d513cb36d10c src/pick.py; 720106c5528a5dcc5fae3971b1ea750173091bf6 data/groups/run_ids.tsv; 

## Amendments
- A1 (17:50, after the R controls, before any T statistic): N's per-phage part fails: held-out R gives avoid calls in
  58/2,555 phages with E >= 5 = 2.27% (bar 2%; results/controls/heldout.txt). Per-phage avoid calls are therefore
  not used in any claim; they are reported only as descriptive counts beside that 2.3% false rate. The pooled tests
  (A, B) and their readings are unchanged; their null and recall pass (results/controls/).
- A2 (17:50): score.py gave p = NaN when the observed rate ratio is exactly 0 (seen in the P3 f = 1.0 run); now p =
  share of permutations with ratio <= 0. f = 1.0 rerun: p < 0.01 in 50/50.
- Controls before freeze: P1 pass (weak, circular; results/p1), P2 pass (pooled 0.588, median 0.547; R TGA 0.001),
  P3 pass (T1 bins f = 0.5: 200/200 p < 0.01, median RR 0.494; T0 bins 200/200, 0.467), N null pass (T1: p < 0.05
  in 5.5% / 6.5%, p < 0.01 in 2.0% / 2.0% for A / B; T0: 6.5% / 3.0%, 2.5% / 0%), P4 pass (skip/O 0.333 vs 0.005).
- A3 (17:59, after the T scores): classes() divided by a zero reference O/E (T1 vs T0: T0 bins with no CGG,
  e.g. bin 24). Guard added: a phage gets no avoid/ancestral class when its bin's reference O/E is 0 or missing
  (score.py sha1 eb0e1ede). Pooled A/B numbers unchanged (T1 vs R lines 1-9 identical on rerun; T1 vs T0 identical to
  the pre-crash output). Descriptive classes only.
- A4 (17:59): test C implemented after the freeze (src/test_c.py), spec unchanged. Post-freeze diagnostics
  (src/diag_skip.py, src/diag_sibskip.py) are descriptive and labelled so wherever used.
