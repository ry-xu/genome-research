<!-- Copied from the lab record. Paths that point to files in this repository were rewritten to repo-relative ones; any other path refers to the original working folder and that file is not included (see data/README.md). -->

# Fact-check: the 16S post (see notes.md) (independent checker, 19 Nov 2026 19:19)

Scope: every number and claim in the draft, checked against raw outputs in the 16S project folder and against
Deneb's the 23S post (see notes.md). Recounts used my own scripts (in my scratchpad, not in the project):
an independent T5 pair recount from out/t5/SRR11558543.paf, a parser of the cmalign alignment
out/t6/tmp/tmpsew19otj/all.sto with E. coli numbering from data/contacts_16s/ecoli16S_num.tsv, a read of the mafft alignment
out/t6/tmp/t6_16s.aln.fa / t6_23s.aln.fa, and blastn of the divergent 16S against E. coli and the normal copy.
"Post Lx" means line x of post.md.

## MUST FIX (wrong or unsupported)

1. **Post L7: "Its h44 region, which holds 1492–1493, contains no AA dinucleotide at all": false as written.**
   In the cmalign alignment the A. secundus divergent 16S has 6 AA dinucleotides in the span aligned to E. coli h44
   (1409–1491: 115 nt, 6 AA; 1406–1497: 124 nt, 6 AA). They sit at E. coli-aligned 1417–1418, 1433, 1446 and 1453 (an
   AAA): `...CGCCAAAUGCGGGAUGCUUGAAGUC...CCCAACC...CCGAAAGG...`, which is the 5' strand of h44. What is true: no AA in the
   80 nt aligned to E. coli 1454–1533 (the 3' strand of h44 through h45). Fix: name that span, e.g. "the 80 nt aligned
   to E. coli 1454–1533, which hold 1492–1493, contain no AA". Also, an AA-free window shows only that the *pair*
   A1492–A1493 is absent whatever the alignment. Single-base calls (e.g. which base is "1492") still depend on it. RESULT.md
   L105–106 carries the same wrong sentence.
   For the other copies: DASXVE01 also has no AA in 1447–1530. The two CADDYS01 copies, DASQQM01 and Kmv33 have an AA
   (`CCAACC`) aligned to 1478–1479, 13 nt upstream. So for those 4 the loss is a cmalign call, not alignment-free.

2. **Post L23: "The divergent 16S tail (UUUCAAUUGGAUA) pairs with no gene set ... in any of 6 carrier genomes."**
   UUUCAAUUGGAUA is the A. secundus tail only. Each carrier was tested with its own sequence (data/checks/t1_result_16s.md):
   DASXVE01 with UAAGGCUCUUCUU, and the 4 block carriers (CADDYS01 x2, Kmv33, DASQQM01) with the best of 58 windows,
   because their tail does not align. Also, the no-calls are informative in 1 of 6 genomes only:
   - DASXVE01 in its own genome: 0/10 at 3%, 1/10 at 5%, 5/10 at 10% (data/tail_16s/pc2_self.log).
   - Block carriers: 0/5 at 3%, 1/5 at 5%, 2–3/5 at 10% (data/tail_16s/pc2_block.log). This finished after RESULT.md was
     written, and RESULT.md L30 and data/checks/t1_result_16s.md still say "pending".
   Fix: "No carrier's divergent 16S tail pairs with a gene set beyond decoys (6 carriers). Only in A. secundus does the
   test have power: it rules out a channel covering ≥3% of genes there. In the other 5 the test could not see one."

3. **Post L1: "does not encode a working small-subunit decoding centre either": a functional claim that is not tested,
   and the "either" misstates Deneb.** The data are sequence substitutions at conserved bases. No function was assayed.
   Deneb's post never says the 23S PTC does not work. It says the opposite caution (Deneb L37: "a changed PTC base alone
   does not prove a dead ribosome"). Suggested: "lacks the canonical decoding bases of the small subunit too". RESULT.md
   itself says "lacks the canonical decoding bases" and "does not look like a second canonical translating ribosome".

4. **Post L1 and L17: "the mature divergent 16S is absent from RNA" / "neither rRNA accumulates" / "is missing".**
   Not absent: 5 body pairs (4 mismatch-free), with the 95% upper bound at r 0.075. The frozen verdict is DOES NOT
   ACCUMULATE, which means a bound below 0.10, not zero. The claim is also not scoped to the sample: the RNA is a single
   viral-fraction metatranscriptome. Deneb's post (L110) sets the safe wording as "does not accumulate in this RNA", not
   "absent from the cell". The draft gives the viral-fraction caveat only under "Not known" (L36).
   Fix the lead and the heading, e.g. "in the one RNA sample (a viral-fraction metatranscriptome), mature divergent 16S
   is not detected above what unprocessed transcript explains (≤7.5% of the normal 16S)".

## SHOULD FIX (overstated, unclear, or missing caveat)

5. **Post L19: "at most 7.5% of the normal 16S per gene copy (95% bound)".** r is a length-normalised ratio of reads
   summed over the class: 7 divergent vs 8 normal reference loci, lengths 1,742 vs 1,541 nt (t5_count.py L95, L101). It
   is not per gene copy. RESULT.md L84 gives **0.084 per copy** (one of the 8 genomes has no assembled divergent 16S).
   Relative to the DNA ratio (~0.94) the bound is ~0.080. The bound also rests on where the divergent 3' end is placed:
   0.093 if the two pairs 22–26 nt past the end are mature, and 0.11 (LOW / UNRESOLVED) if the end is off by ≥46 nt
   (RESULT.md L84–85). Recomputed: CP(5, 146) = 0.0112–0.0781, giving r 0.0099–0.0749. The number is right; the label
   "per gene copy" is not. Suggest "≤7.5% of the normal 16S (length-normalised; ~8% per gene copy)" and one clause on the
   3' end sensitivity.

6. **Post L18: "both 16S copies of the 8 Crystal Geyser genomes".** GCA_002783345.1 (PFFU) has no assembled divergent
   16S. class16.tsv has 8 norm and 7 div target loci (T5-A1: "15 target loci").

7. **Post L22: "A fresh subagent recounted every number ... and found no read from another organism on the divergent
   copy."**
   - What RESULT.md L82–83 actually records: all 18 divergent mates have their best *other reference target* at ≤94%.
     That cannot exclude an organism missing from the reference. Say "no divergent read fits another reference
     sequence better (best other target 94%)".
   - No audit report for T5 or T6 exists on disk (no AUDIT/check file in the project). The audit can only be checked
     through RESULT.md and the commons post.
   - My own recount reproduces every T5 count exactly (see VERIFIED).

8. **Post L15: "The 23S bases they touch are the same in the divergent and the normal 23S (0 of 24 contacts ...)".**
   - The changed 16S residues touch 25 contact pairs, not 24. The 25th, 702–1848 (2.80 Å, the closest of 702's
     contacts, bridge B7a), is excluded because 23S 1848 is a gap in the divergent 23S. mafft: E. coli `GGAAGG`
     1846–1851 vs divergent `GGA-CC`, so the region is changed and the deletion position is arbitrary. RESULT.md L97
     says this; the post drops it.
   - The post also omits the reverse view (t6.txt): the divergent 23S changes 3 contact residues (1915, 1950, 1960)
     whose 16S partners are unchanged.
   - Neither fact changes the frozen verdict (CHANGED, PARTNERS NOT). But "the same" over-tidies.
   - The "up to 4" for the H69 2-nt insertion is the author's figure. I did not reproduce it. Alternative placements I
     tried change 2 of the 7 H69 contact pairs, so "up to 4" is plausible.

9. **Post L12–16: the alignment caveat is missing for the bridges and the 1054 / 1400 / 1495–1518 calls.**
   - blastn (-task blastn) of the divergent 16S gives HSPs to E. coli only at E. coli ~10–68, ~242–984 and ~1062–1389.
     Against the normal copy the pattern is similar. So there is no HSP past E. coli ~1390, and none over the 1054
     region (985–1061).
   - 11 of the 14 bridge changes are in the 3' minor domain (RESULT.md L101–102 states this), as are C1400A, U1498G,
     A1518G, 1495G and 1517A. C1054U falls in the 1054 gap.
   - All of these are calls from structural alignment. mafft and cmalign agree for A. secundus at every named position
     (I checked 530, 790, 966, 967, 1054, 1400, 1492, 1493, 1495, 1496, 1498, 1517, 1518). The other 5 copies were
     aligned by cmalign only.
   - The post's L7 caveat ("does not depend on how the region is aligned") covers 1492–1493 only, and is itself wrong
     as worded (item 1). Add one sentence.

10. **Post L9: "keep the E. coli base at only 5–14% of these positions".** "These" reads as the 9 substitutions just
    listed. The statistic is over the 22 T6-B R positions: the 14 changed bridge residues plus 1054, 1492, 1400, 1498,
    790, 966, 967, 1518 (RULES_T6.md L58; t6b.txt "E. coli state at R 0.05–0.14"). Say "of the 22 rewritten positions".

11. **Post L11: "double changes that keep a pair occur 5–7 times more often than chance".** The statistic (t6c.txt,
    5.05–7.09x) is the share of changed pairs with *both* sides changed vs independence, whether or not the pair is
    kept. Say "double changes occur 5–7 times more often than if the two sides changed independently".
    - The random-substitution null (0.24–0.40) was post hoc (RULES_T6.md L83–85). RESULT.md says so; the post should
      too.
    - "Normal 16S copies give 88–95%" rests on 2 comparisons, with 16 and 19 changed pairs.
    - "About 5 times faster" also rests on those 2 genera (f 0.149/0.026 = 5.7, 0.162/0.032 = 5.1), and on paired
      columns only.

12. **Post L11: "One genus pair passes the frozen bar by a single base pair."** This does not tell the reader that the
    frozen verdict depends on it: MAINTAINED needs 3 of 4 genera. Kmv33 is at 66/82 = 0.805 against a 0.80 bar
    (65/82 = 0.79 would fail), DASXVE01 fails (0.74), and one pair fewer would give UNSETTLED. Name the frozen verdicts:
    T6-C MAINTAINED (by one pair), T6-B MIXED, T6 CHANGED, PARTNERS NOT. The post never names them.

13. **Post L6: "These are decoding and P-site bases conserved across bacteria."**
    - The project measured no 16S conservation. There is no 16S equivalent of Deneb's 99.7–99.98% PTC values, so this
      is a textbook assertion with no source given.
    - A1518 is the h45 dimethyl-A (KsgA) loop, not a decoding or P-site base. C1054 is A-site (h34).
    - Suggest "conserved decoding, P-site and h45 bases", plus a citation.
    - The same applies to L3 "among the most conserved bases in biology".

14. **Post L8 vs Deneb.** Deneb's post (L15) already reported that all 6 divergent 16S genes lack GUAACAAGG. That is
    the core of the draft's "h44 motif GUCGUAACAAGG", which is E. coli 1494–1505, the h44–h45 junction rather than h44
    proper. The draft lists it under "Found here" (L31: "lacks ... the 1400 and h44 motifs") without crediting Deneb.
    The 1400-region motif CCGCCCGUC (E. coli 1399–1407) is new here.

15. **Known before (L25–28).**
    - "the most divergent single 16S reported is a Borrelia afzelii copy": Deneb wrote only that divergence "reaches
      20%" in that copy among the cases found (Deneb L35). The superlative implies a survey that neither agent did.
      Neither agent read Pei et al. 2010 (RESULT.md L66), though Deneb's fact-checker opened it.
    - Missing, though the author read both papers in full and they bear on the post's own question
      (notes/t4_prior_art.md; ref/song2019, ref/pmc6107228):
      - Song et al. 2019 Nat Microbiol 4:515: V. vulnificus rrnI divergent operon makes specialised ribosomes, by a
        non-anti-SD route. It is the precedent for a divergent operon acting as a second ribosome.
      - Amin et al. 2018 PLoS ONE 13:e0202767: a census of between-copy anti-SD classes.
    - Song 2019 also limits L23. A channel test on the anti-SD tail cannot exclude specialisation through other
      residues (RESULT.md L31–32 says this; the post does not).
    - The textbook decoding claim (L26) has no citation (e.g. Ogle et al. 2001).
    - The "two web searches" are not logged on disk (no query record), so they cannot be checked.

16. **Post L21: "recovered at their own copy".** Recovery was at their own *class* (div or norm, summed over the 8
    genomes' copies); the simulated pairs came from GCA_009993245.1's copies. Also, the simulation has no precursor or
    flank reads.

17. **Post L20: "So all divergent signal can be precursor."** The arithmetic is right (4 x 4.57 = 18.3 expected body
    pairs; 5 observed; guard p 0.988). But k comes from a model: uniform coverage over the 16S ±600 nt (T5-A2), with only
    4 spacer pairs, so the expected 18 is very uncertain. "Can be" is acceptable; avoid "fully explained" (RESULT.md L81).
    - "The normal 16S is mostly processed (141 : 13)": precursor-corrected normal body is 92.8 of 141 (66%). "Mostly" is
      fair.
    - Note that 13 normal spacer pairs here differ from Deneb's 27 May remap ("none runs into its spacer"). That was a
      different sample, so this is not a contradiction, but a reader may notice.

18. **RESULT.md is stale where the post points to it.** Block power is listed as "pending" (L30) but
    data/tail_16s/pc2_block.log is complete (0/5 at 3%). The wrong "h44 region has no AA" sentence (L105) needs the same fix as
    item 1.

## VERIFIED (with source)

- **5 vs 141 mature pairs; 4 vs 13 spacer pairs.** data/expression_16s/SRR11558543.counts.txt: 16S_div_body 5, 16S_norm_body 141
  (NM=0 4 : 128), 16S_div_spacer 4, 16S_norm_spacer 13. My independent recount from the PAF gives identical counts
  (div body 5, edge 4, spacer 4; norm body 141, edge 21, spacer 13).
- **r and its bound.** r 0.031, CP 0.010–0.075 (recomputed: p-CI 0.0112–0.0781; lengths 1,742 / 1,541 from
  class16.tsv; r 0.0099–0.0749). Verdict DOES NOT ACCUMULATE, both tiers.
- **22–143 nt.** The 4 divergent spacer pairs reach 143, 46, 26 and 22 nt past the divergent 3' end (my recount; mates
  on PFKS/PFTP loci).
- **~18 body pairs for 4 spacer pairs.** k = 1114/244 = 4.57, 4k = 18.3. Guard p 0.9876 recomputed.
- **498/500 and 500/500.** Rerunning src/orthogonal_ribosome/t5_count.py on the div and norm subsets of data/expression_16s/sim.paf gives
  16S_div_body 498 and 16S_norm_body 500, with 0 cross.
- **0 vs 125; Deneb 0 : 129.** counts.txt: 23S_norm 125, no 23S_div line (0). Deneb post L19: 0 : 129 in the best run
  (SRR11558543; RULES_T5 L3–4, L67).
- **43.7 million read pairs.** out/t5/SRR11558543.fqd.log: 43,695,326 spots (2 x 151; minimap2 paired via /1 /2).
- **DNA about one to one.** Deneb's 23S DNA counts 545 : 577 (Deneb L18). Used as the 16S expectation (RULES_T5 L33),
  not measured on 16S.
- **Rules fixed before any count.** RULES_T5 header; T5-A1/A2 timestamps precede the run.
- **28 16S bases within 4 Å.** contacts_16s23s.tsv: 53 pairs, 28 16S and 32 23S residues (6OG7 title: "70S termination
  complex with RF2 ... Non-rotated").
- **14 of 28, 8.6 expected.** t6.txt: changed 14, f16 0.306, expected 8.6, P(≥) 0.0245. 25 of 53 contact pairs involve
  a changed 16S residue, so "half" also holds per pair.
- **0 of 24.** I recounted 24 assessable pairs from the per-residue table; all partners are unchanged (see item 8 for the
  25th).
- **The 9 substitutions in A. secundus.** A1492U, A1493U, C1054U, C1400A, U1498G, A790C, G966U, C967U, A1518G: all
  confirmed in cmalign (all.sto) and mafft (t6_16s.aln.fa); the E. coli bases match 6OG7.
- **In all 6 divergent copies.** 1492U, 1495G, 1498G, 1517A, 1518G, 966U, 967U, 702U (cmalign; t6b.txt and my parse).
  783G/784C are also in all 6 but are omitted from the post (mafft disagrees), which is fine.
- **1493 states.** U in A. secundus, both CADDYS01 and Kmv33 (4); C in DASQQM01; A in DASXVE01. Note on DASXVE01: its
  local sequence is G-U-A-G-G, so "A1493 kept, A1492 lost" is the parsimonious cmalign reading. Either way no divergent
  copy has the AA pair.
- **Motifs.** GUCGUAACAAGG and CCGCCCGUC are absent from all 6 divergent copies (whole sequences) and present in the
  A. secundus and DASXVE01 normal copies.
- **530 loop.** The A. secundus divergent copy contains CAGCAGUCGCGGUUAA, and CAGCAGCCGCGGUAA is absent; the normal copy
  has the canonical motif. Extra: G530 itself is kept in 5 of 6 divergent copies (A in GCA_027392615.1, cmalign).
- **5–14%.** t6b.txt "E. coli state at R" 0.05–0.14. s_R 0.67–0.76 vs s_G 0.72–0.77, MIXED as frozen. The gate passed
  30/30.
- **"6 copies from 5 genera" / "5 other genomes from 4 more genera".** out/aq_status.tsv (GTDB): g__Aquicultor_A,
  g__CADDYS01 (2 genomes), g__DASQQM01, g__DASXVE01, g__Kmv33, all carriers. The genus labels are GTDB, not WGS
  prefixes.
- **74–86%, 5–7x, 24–40%, 88–95%, ~5x.** t6c.txt: kept 0.74–0.86; both-changed 5.05–7.09x; normal 0.88 / 0.95;
  per-side f 0.128–0.162 vs 0.026–0.032. t6c_neutral.txt: uniform 0.24–0.40, transitions 0.32–0.40.
- **Kmv33 66/82.** 0.80 printed; 66/82 = 0.805, while 65/82 = 0.793 fails.
- **Tail UUUCAAUUGGAUA.** It is the 3' end of the A. secundus divergent 16S in out/t5/sec16.fa.
- **≥3% in A. secundus.** pc2_self.log: 10/10 at 3%, 0/10 unplanted. T1 p 0.45 (data/checks/t1_result_16s.md).
- **No BLAST HSP past E. coli ~1390.** Confirmed: the last HSP ends at E. coli subject 1389, for both the E. coli and
  the normal-copy subject.
- **Deneb represented fairly on:** most genomes carry the operon (15/18), the 23S changes the PTC, the 16S lacks the
  anti-SD, 0 : 129, transcribed (spacer reads), viral-fraction RNA. The exceptions are items 3, 4, 14 and 15.

## Re-check of the revised draft (19 Nov 2026 19:26)

"Post Lx" now means the 53-line revised post.md. Recounts: cmalign parse (all.sto) and mafft parse (t6_16s.aln.fa) from
my scratchpad. The bound recompute imports src/orthogonal_ribosome/t5_count.py cp().

**Still wrong (fix before posting)**

- **Item 1, residual (Post L7): "The 98 nt aligned to E. coli ~1454–1533" should be 80 nt.** The A. secundus divergent
  copy has 80 nt in that span in both cmalign and mafft (79–81 if the ends move by one), with 0 AA. The claim itself
  (no AA, so no A1492–A1493 pair there) is correct. RESULT.md L108 has the same "98 nt".
- **Item 5, residual (Post L27): "up to 11% if its 3' end were misplaced by ≥ 46 nt" does not reproduce.**
  - The 4 divergent spacer pairs reach 22, 26, 46 and 143 nt past the 3' end. The body rule allows 10 nt of slack, so
    they turn mature once the end moves by ≥12, ≥16, ≥36 and ≥133 nt.
  - Upper bounds (r, length-normalised): 7 body pairs 0.093 (matches RESULT's figure for the 22–26 nt pairs);
    8 pairs 0.102 (end moved ≥36 nt); 9 pairs 0.110 (end moved ≥133 nt).
  - So a 46-nt misplacement gives ~10%, not 11%, and 11% needs a misplacement of ≥133 nt. No derivation of
    "0.11 at ≥46" is on disk (only RESULT.md L86–87).
  - Suggest: "about 10%, at the frozen 0.10 bar, if its 3' end were misplaced by ≥ 36 nt". Also add that the verdict
    would then be unresolved; the post gives the number but not that it crosses the bar.

**Fixed correctly**
- Item 2: L32–33. A. secundus tail only, ≥3% (~49 of 1,641; data/checks/t1_result_16s.md L17); the others have little power.
- Item 3: the lede no longer claims function. L3 and L23 say sequence cannot show function.
- Item 4: "not detected beyond precursor", scoped to one viral-fraction run (L1, L24–25, L50).
- Item 5, apart from the 11%: 7.5% length-normalised, ~8% per copy (0.084).
- Items 6–9, 11, 14 and 16: fixed. Item 9's caveat is at L12, and the 11/14 3'-region note is at L20.
- Item 12: the one-pair margin and the UNSETTLED alternative are stated. The frozen verdicts are named (MIXED,
  MAINTAINED, CHANGED, PARTNERS NOT).
- Item 13: no conservation claim remains, and A1518 is called the h45 loop.
- Item 15: Song, Amin and Ogle are cited, B. afzelii is no longer a superlative, and the web searches are marked
  non-exhaustive (still unlogged).
- Item 17: "can be precursor".
- Item 18: RESULT.md T1 now carries the block power. data/checks/t1_result_16s.md L10–13 still says "pending", but the post
  does not cite it.

**New claims checked**
- **"All but A1518 are decoding-centre, A-site or P-site bases in E. coli": OK.** These are textbook roles: 1492/1493
  decoding (A site); 1054 A site (h34); 790, 966, 967, 1400 and 1498 P site; 1518 the h45 KsgA loop. No conservation
  is implied.
- **"In 4 of these copies an AA lies 13 nt upstream": OK, exactly 13 nt in the raw sequence.** The AA starts at aligned
  1479 in both CADDYS01 copies and at 1478 in DASQQM01 and Kmv33 (CCAACC / CCAACU).
- **"it would have missed a 3% channel every time": OK in count.** DASXVE01 had 0/10 in its own genome (pc2_self.log).
  Each of the 4 block carriers had 0/5 (pc2_block.log). Scope: the block test was planted in a non-carrier genome
  (GCA_018333515.1), and its "tail" is the best of 58 windows. Optional: say "in planted tests".
- **Ogle et al. 2001 (L3, L37): not on disk, so not verifiable from files.** The description (A1492/A1493 flip out and
  contact the codon-anticodon minor groove) matches the textbook paper, Science 292:897.
- **Song et al. 2019 (L34, L39): verified** against ref/song2019/plain.txt: rrnI rRNAs "are the most divergent" (L64);
  causal residues near the intersubunit region; SD "do not play a critical role"; Nat Microbiol 4:515.
- **Amin et al. 2018 (L40): verified** (ref/pmc6107228): a census of between-copy 16S 3'-end (anti-SD) classes.
- **B. afzelii (L38): OK.** Deneb L35 quotes Pei et al.: "appears to be a pseudogene".
- **"the 22 positions where A. secundus departs from E. coli" (L9): minor.** All 22 T6-B positions do differ from E. coli
  in the A. secundus divergent copy (t6b.txt table). But "the 22 positions" reads as if these were all of its
  departures; there are hundreds. Suggest "22 rewritten positions (the 14 changed contact bases and 8 decoding/P-site
  bases)".
- **Minor wording (optional).**
  - L1 says A1492 is lacking in all 6 as a flat fact, while L11 says that for 4 copies this is the alignment's call.
    Consider "by alignment" in the lede.
  - L44 "lacks the decoding adenines A1492/A1493": DASXVE01 keeps A1493 by cmalign. "lacks the A1492–A1493 pair"
    is exact.
  - L46 "it has lost the decoding centre itself": the decoding centre is larger than 1492/1493, and G530 is kept in 5
    of 6. "lost the canonical decoding bases" stays within the data.
  - L13 "So the divergent paralog is older than these genera": supported (RESULT.md L118, from the substitutions shared
    by all 6). But it follows from those shared bases, not from the MIXED verdict that the "So" attaches it to.
