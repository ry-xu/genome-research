<!-- Copied from the lab record. Paths that point to files in this repository were rewritten to repo-relative ones; any other path refers to the original working folder and that file is not included (see data/README.md). -->

# T1 channel test: result (16:42 19 Nov; RULES.md v1 + A1-A7; data/tail_16s/carriers.jsonl, data/tail_16s/t1_tails.log, data/tail_16s/t1_blocks.log)

Controls: PC1 pass, PC2 frozen FAIL (A3, power read per tail instead), NC1 pass (0/120 false calls), independent code
audit clean (A6 fixes applied before NC1/T1).

| carrier | genus | tested | E (exclusive share) | decoy mean | p | z | pos ratio | call | detection limit |
|---|---|---|---|---|---|---|---|---|---|
| GCA_009993245.1 | Aquicultor_A (A. secundus) | tail UUUCAAUUGGAUA | 0.0037 | 0.0042 | 0.45 | -0.13 | 1.40 | no | 3% (in own genome and non-carrier) |
| GCA_036504445.1 | DASXVE01 | tail UAAGGCUCUUCUU | 0.0213 | 0.0104 | 0.085 | 1.60 | 1.03 | no | none within 10% (uninformative) |
| GCA_902810715.1 | CADDYS01 | block, best GAGGUGUGGGCGU | 0.0524 | 0.0412 | 0.082 | 1.43 | 0.90 | no | 0/5 at 3 %, 2/5 at 10 % |
| GCA_014859275.1 | Kmv33 | block, best GGUGUGGGCUUAC | 0.0591 | 0.0444 | 0.047 | 1.66 | 1.00 | no | 0/5 at 3 %, 3/5 at 10 % |
| GCA_035643475.1 | DASQQM01 | block, best GAGGUGUGGGCGA | 0.0394 | 0.0333 | 0.18 | 0.98 | 1.10 | no | 0/5 at 3 %, 2/5 at 10 % |
| GCA_027392615.1 | CADDYS01 | block, best GGUGUGGGCGUCC | 0.0496 | 0.0394 | 0.090 | 1.45 | 0.93 | no | 0/5 at 3 %, 3/5 at 10 % |

**Verdict: no channel.** No carrier is called. The channel call needed calls in >= 2 genera. The specificity controls
and the ortholog-overlap step therefore do not run (frozen: only for called tails).
- Scope: for the A. secundus tail, no private SD channel covers >= 3% of genes (~49 of 1,641). A smaller regulon
  (a handful of genes) is below what this test can see. For DASXVE01 the test has no power. For the four block carriers
  (divergent copy's tail unaligned), see block power.
- In all four block carriers the best window is the same G-rich element (GGUGUGGGC), with no positional enrichment
  (ratio 0.90-1.10). This fits the G-U-wobble background that makes the decoy null heavy-tailed (A3), not a channel.
