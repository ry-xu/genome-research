<!-- Copied from the lab record. Paths that point to files in this repository were rewritten to repo-relative ones; any other path refers to the original working folder and that file is not included (see data/README.md). -->

# Detection limits for T1 (16:14 19 Nov; RULES A3/A5)

PC2 as frozen (random tails, GC 5-8, planted at -14..-7 in non-carrier GCA_018333515.1) FAILS: 3/20 calls at 3%, 0/20 at 1%,
5/20 at 5%, 11/20 at 10%, 13/20 at 20% (out/pc2_result.jsonl). Random-tail power is poor because decoy 13-mers vary
widely in pairing propensity (G-rich ones pair with 10-20% of windows), so a typical random tail sits mid-null.

Tail-specific power for the two defined tails T1 tests (src/orthogonal_ribosome/pc2_tail.py, 10 replicates, seeds 5000+rep; out/pc2_tail_result.jsonl):

| tail | 0% | 3% | 5% | 10% | median planted exclusive |
|---|---|---|---|---|---|
| A. secundus UUUCAAUUGGAUA | 0/10 (z -0.77) | 10/10 (z 5.7) | 10/10 (z 9.1) | 10/10 (z 14.9) | 0.76-0.78 |
| DASXVE01 UAAGGCUCUUCUU | 0/10 (z 2.04) | 10/10 (z 3.5) | 4/10 (z 3.8; was 8/10 before A6) | 8/10 (z 4.1; was 9/10) | 0.81-0.84 |

**Corrected 16:20 (RULES A6, strict p <= 0.01 gate from the audit).** A. secundus is unchanged and its detection limit is 3%
(max p 0.006 at 3%, 0.001 at 5% and 10%; z 5.7+). DASXVE01: p stays at 0.007-0.012 at 3%, 5% and 10% alike (out/pc2_tail_result.jsonl). About
ten GC-matched decoys rise with the plant, so calls are coin flips at the p boundary. The A5 rule gives 3%, but it is
not a real detection limit, and a T1 non-call for the DASXVE01 tail is uninformative.
(Superseded: "Detection limit (smallest share called in >= 8/10) = 3% planted genes for both.") The 0% rows are one identical unplanted
genome (decoy seed fixed). Note the DASXVE01 tail already reaches z 2.0 unplanted in the non-carrier: its AGGCUCUU
pairs with purine-rich, SD-like windows that the canonical-weak gate does not remove, so its margin above
background is thin and the in-family non-carrier control matters for it.
Block carriers: src/orthogonal_ribosome/pc2_block.py (A5), pending.

## In-genome power (A7, 16:37; src/orthogonal_ribosome/pc2_self.py, out/pc2_self_result.jsonl)
| tail, planted in its own carrier | 0% | 3% | 5% | 10% |
|---|---|---|---|---|
| A. secundus in GCA_009993245.1 | 0/10 (= the T1 value, p 0.45) | 10/10 (max p 0.002, min z 4.5) | 10/10 | 10/10 (z 11) |
| DASXVE01 in GCA_036504445.1 | 0/10 (p 0.085) | 0/10 (min p 0.017) | 1/10 | 5/10 |

Reading: A. secundus no call = no channel covering >= 3% of genes (~49 of 1,641), in both genomes. DASXVE01: no detection
limit within 10%; its no-call is uninformative.
