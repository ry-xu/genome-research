# REFEREE: The phage agent 18:27, TANB77-host phages carry the CGG +1 skip signature

The referee, 15 Nov 2026. RULES.md (sha1 75ff0ede1de0) frozen 18:45, before any split was scored. My own code (src/).

## Verdict: CONFIRMED out of sample, with an independent host predictor. Two caveats carried.

- **Self-check.** src/skip.py (own code, pyrodigal 3.7.1 meta, table-11 CDS) on code-walls/data/phages.fna, low GC:
  T1 490/834, R 115/5835. That is kill1's numbers exactly (selfcheck.tsv).
- **H4 (decisive, out_h4.txt).** 49 UHGV phages with no CRISPR host whose k-mer host (kmers_host_lineage_gtdb_r207)
  is o__TANB77 (KT), against 230 GC-matched phages with no CRISPR host and a k-mer host in another Clostridia
  order (KR, seed 15). None of these phages is in the phage agent's sets; sequences come from votus_full.fna.zst.
  - At low GC (24-30): KT 20/50 = 0.400 (24 phages; sites in 11 of the 17 that carry CGG, at most 4 per phage).
    KR 1/146 = 0.007. MH 59.5. **REPRODUCES** under the frozen bar (>= 0.2, >= 30 CGG, MH >= 5).
  - Lineage-novel KT (vGENUS holds no T1/T0/Tf phage), low GC: 7/28 = 0.25 (12 phages). Descriptive.
  - At all GC: KT 24/366 = 0.066, KR 17/1786 = 0.010, MH 7.7. The high-GC KT phages show no signal (4/316).
- **H1-H3 (out_h135.txt, kill1 site table + k-mer host).**
  - T1 at GC >= 32: 2/42 have a TANB77 k-mer host, against 304/394 at low GC.
  - The TANB77 host genomes at GC >= 32 that skip are CAG-508 CALYAK01, CAMCHB01 and WQVL01 (32.6-35.3%, call 1 in
    frameshift-codes scratch/cgg_test_T.tsv). No UHGV phage has a CRISPR host in those genera.
  - So the signal-free high-GC "TANB77 phages" have hosts in low-GC genera, and only 2 of them have a TANB77 k-mer
    host. Wrong or promiscuous host calls fit better than a GC-limited trait, but neither is proven.
  - Where CRISPR and k-mer hosts disagree at low GC, the signature follows the CRISPR host both ways:
    - R phages with a TANB77 k-mer host: 3/141 = 0.021;
    - T1 phages with a non-TANB77 k-mer host: 26/39 = 0.667.
- **H5.** Only 4 viral genera hold both T and R phages with CGG. Only one is low-GC (vGENUS-06445: T 3/5, R 1/4).
  - So at genus level the signature cannot be separated from phage lineage: low-GC TANB77-phage genera are host-pure.
  - Kill1's family-level contrast (outside vFAM-00025: 344/646 vs 11/355) is the best lineage control available.
  - The 12 lineage-novel KT phages (0.25) add a little.

## Caveats for any outbox text
1. The phage signature cannot be separated from phage genus, because no low-GC genus spans hosts. "Follows the
   host" rests on the family-level contrast plus the independent k-mer predictor. It is not a within-genus host switch.
2. The high-GC absence is unexplained either way. The KT phages (k-mer host TANB77) at GC >= 32 lack it too.
3. The host-level skip itself is the 11 Nov outbox claim. another agent's junction-peptide search is still pending.

## Outbox paragraph check (18:55; work/code-walls/notes/outbox_phage_para.md, drafted 18:57)
Checked each number against its source file. CLEARED once the following fixes are made.
1. "at 24-30% GC": the bins are 24-30, which means GC 24.00-31.99%. Say "below 32% GC" (or "24-32%").
2. "other codons 0.7-3.6%": diag_sibskip.py tests 9 related codons only (CGA CGT CGC AGG AGA GGG CCG GGA AAG):
   T1 0.0069-0.0360. Say "nine related codons (the other Arg codons, GGG, CCG, GGA, AAG)".
3. Read-through 31/40: add the controls, so the reader has a baseline. These are 3/40 at T1 gene ends on other
   codons and 8/40 at R CGG sites (kill1 Task 3, R-phage ORFs).
4. Out-of-sample sentence: "whose DNA composition matches TANB77 genomes" is wrong. It is UHGV's k-mer host
   prediction (kmers_host_lineage_gtdb_r207), which does not use spacers. Also add the scope: "below 32% GC" and
   "phages whose k-mer host is another Clostridia order, matched for GC".
5. High-GC limit: "their composition points to other orders" overstates it. Of the 42, 21 have a k-mer host
   outside TANB77, 19 have none and 2 are TANB77 (my H3). And the k-mer-TANB77 phages above 32% GC also lack the
   signature (4/316). That cuts against the wrong-host reading, so it must be stated next to it.
   "1-4 spacer links" checks out: 33/42 (crispr_total_connections; 24 have 1).
Checked and fine: 436; 524/896; 116/5,874; RR_use 1.07; Gly 27/38; the Pfam-check agent 239/600 vs 18/600; KT 20/50 vs 1/146;
"all 5 families with >= 30 CGG on both sides split by host" (results/posthoc/vfam.txt: vFAM-00026, -00003, -00001,
-00015, -00006; T 0.43-0.64 vs R 0-1/45); vFAM-00025 R 2/8.
- 18:56: revised paragraph re-read; all five fixes in. The high-GC split 18 other order / 3 Clostridia with no order / 19 none / 2 TANB77 is checked against meta_cols.tsv (my H3 "21 outside" = 18 + 3). **CLEARED.** The post is held by the main analysis agent (agents/oregano/STATE.md line 20; not a live session); posting is Ryan's call.
- 18:57: The phage agent wrote outbox/2026-11-11-cgg-skip/addendum-phages.md. I diffed its body against the cleared paragraph: identical. post.md is untouched.

## H6 (frozen amendment, RULES sha 816ffcf4df30), run 19:15: SURVIVES
- After the Pfam-check agent 19:09 (effective n). Script src/h6.py, output out_h6.txt (3.9 s, 0.46 GB).
- Windows: T 2,458, R 5,834, KT 50, KR 146. Clusters (all pooled): 4,481.
- H6a, (cluster, side) units: T 353/682 = 0.518 vs R 57/3,662 = 0.0156; ratio 33.3; Fisher p 9.0e-259 -> SURVIVES.
  - Single-vGENUS clusters: 218/443 vs 49/3,057; ratio 30.7; p 3.4e-169.
- H6b: 6 mixed T/R clusters, each with 1 R CGG; T-skip-only 0, R-skip-only 0. No test of a host switch at one site is possible.
- H6c: 12 of the 20 KT skip CGGs fall in clusters holding T1/Tf CGGs. Outside them, KT 8/31 vs KR 1/146 (Fisher p 3.1e-6). Out-of-sample holds for phages; the sites are partly shared with the training set.
- Descriptive T1-only (src/h6_t1.py, T1+R clustered alone) reproduces the phage agent's sites_t1.txt exactly: 202/384, one-genus 137/281 (CGGs 159/311), 27 sites in ≥5 genera hold 300/832.
  - kill1's re-called genes give 834 T1 CGGs (490 skip) vs main 896 (524).
- Addendum sentence: cleared to the phage agent with 3 fixes ("Many CGGs", the 832 vs 896 clause, and the H6c clause for the out-of-sample bullet).
- 19:16: read addendum-phages.md after the phage agent's edit (sha 73b6943220bb). The Shared sites bullet and the out-of-sample clause match out_h6.txt / out_h6_t1.txt; all 3 fixes are in. CLEARED.
