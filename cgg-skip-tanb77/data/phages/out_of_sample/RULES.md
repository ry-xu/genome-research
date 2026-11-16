# RULES: referee check of the phage agent 18:27 (TANB77-host phages carry the CGG +1 skip signature)

Written by the referee, 15 Nov 2026, before any skip number for the splits below was computed. Only set sizes were
counted (49 out-of-sample phages, none in code-walls/data/phages.fna). Inputs: UHGV v1.0 metadata
(work/virus-reading/data/uhgv_dl/votus_metadata_extended.tsv.zst, columns cut to meta_cols.tsv), kill1's per-site
table work/code-walls/kill1/out/meta.cgg_sites.tsv (pyrodigal meta, table-11 CDS only), and votus_full.fna.zst
(Zenodo 10.5281/zenodo.17402089) for the out-of-sample set. Skip site, as frozen by the phage agent: in-frame CGG in the last
16 sense codons of a complete CDS >= 300 nt whose +1 frame read from the CGG's 2nd base has no stop in 50 codons.
GC bin = even floor of gc_percent. "Low GC" = bins 24-30.

The question: does the signature follow the host when the host comes from a predictor independent of CRISPR spacers
(UHGV k-mer host, column kmers_host_lineage_gtdb_r207)?

## H1-H3 (zero CPU, kill1 site table joined to the k-mer host)
- H1: R phages (CRISPR host another Clostridia order) split by k-mer host: order TANB77 vs anything else. skip/O each,
  low GC and all bins.
- H2: T1 phages split by k-mer host: TANB77 / other / none. skip/O each, low GC.
- H3: T1 phages at GC >= 32: share with k-mer host in TANB77, vs the same share at low GC (tests "host calls wrong").
All descriptive; reported with counts.

## H4 (out-of-sample, the decisive one)
- KT: UHGV phages with no CRISPR host, k-mer host in o__TANB77, CheckV Complete/High/Medium, genetic_code 11,
  >= 5 kb (49 phages).
- KR: same filters, k-mer host in c__Clostridia but not o__TANB77, no CRISPR host; per GC bin up to 5x the KT count
  of that bin (seed 15), bins restricted to those holding KT phages.
- Gene calls: pyrodigal 3.x meta mode, table-11 CDS only (as kill1); skip sites by my own code (src/skip.py), which
  must first reproduce kill1's T1 and R low-GC counts (490/834, 115/5835) on code-walls/data/phages.fna within 1%.
- Readings, low-GC bins (24-30):
  - REPRODUCES if KT skip/O >= 0.2, KT has >= 30 in-frame CGG, and MH ratio KT vs KR (strata GC bins) >= 5;
  - FAILS if KT skip/O < 0.1 with >= 30 CGG;
  - otherwise UNRESOLVED. Also reported: all bins; KT restricted to phages whose UHGV vGENUS holds no T1/T0/Tf phage
    (lineage-novel), descriptive.

## H5 (zero CPU): within viral genus
vGENUS (from uhgv_taxonomy) holding >= 1 T (T1+T0+Tf) and >= 1 R phage with in-frame CGG: pooled skip/O per side and MH
ratio with vGENUS strata. Descriptive; kill1 did the family-level version.

## Machine
1 thread, wait while load >= 10, download <= 1 GB (votus_full.fna.zst, resumable), filtered to the needed ids then deleted.

## Amendment H6 (written 19:14 15 Nov, after the Pfam-check agent 19:09; before any cluster count below was computed)
Effective n. The Pfam-check agent showed 964/1,543 T1+Tf skip CGGs fall in recurring conserved sites. Test whether the host
contrast survives with sites collapsed.
- Windows: every in-frame CGG at low GC (bins 24-30) in T1, Tf and R from kill1's meta.cgg_sites.tsv, plus KT and
  KR from my h4_skip run (same rule): 15 nt + CGG + 15 nt in gene orientation. Single-linkage clusters at Hamming
  <= 6 over 33 nt (the Pfam-check agent's rule), all sets pooled.
- H6a: unit = (cluster, side), side T = T1+Tf, side R = R. A unit is "skip" if >= half its CGGs meet the signature.
  Report skip units / units per side and Fisher exact. **SURVIVES** if T share >= 5x R share and p < 1e-3; else
  **WEAKENED**. Also reported: the same with units restricted to clusters in one viral genus (singletons).
- H6b (site-level host switch): clusters holding CGGs from both T and R phages. Per side: skip CGGs / CGGs, and
  the number of such clusters where T is skip-majority and R is not, and the reverse. Descriptive (n unknown beforehand).
- H6c: KT skip sites: how many fall in a cluster that holds a T1/Tf CGG (i.e. not new sites). Descriptive.
