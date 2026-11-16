# Follow-up S: independent alignment-free re-count (stranded runs)

**Status (2026-11-12 10:40): DONE. Decision: transcripts keep the C.**
- The 8 runs of `../followS_runs.txt` were counted pooled, spots 1..20,000,000 of each (160,000,000 spots, 320,000,000
  reads). The exact-spot check passed for every run.
- Set A: 6,731 observations at 251 sites, 3 lacks, **D_A = 0.00045 [0, 0.0018]**.
- Set C: 14,577 observations at 2,084 sites, 0 lacks, **D_C = 0 [0, 0]**.
- The test is powered, and hi_A = 0.0018 < 0.10, so the rule's decision is "keep".
- **All 3 lacks are at one Set A site** (MGYG000299934_7:8068, in ffh_1), which has 0 keeps. It looks like a
  site-specific genome difference, not an RNA-wide effect (item 1 below). Without it, Set A has 0 lacks in 6,728
  observations at 250 sites.
- **With the old re-check's SRR37726778 counts added** (reported only): Set A has 8,312 observations at 271 sites, 3
  lacks, D_A = 0.00036 [0, 0.0014]. Set C has 19,045 observations at 2,250 sites, D_C = 0 [0, 0]. The decision is
  again "keep".

I did not read `src/rna_cgg.py`, `src/rna_stream.py`, anything under `data/rna_cgg/`, any `*.obs.tsv` or
`*.totals.txt`, or the README Log entries after "## Follow-up S". The main count was not seen.

## Method

- `count_run.py`, `units.py` and `aggregate.py` are copied from `../recheck/`. `CHANGES.md` lists every change.
  - The changes are paths only, plus two reported extras in `aggregate.py`: `per_study`, and a `--plus-old` flag that
    adds the old SRR37726778 units.
- Matching, exact-spot check, bootstrap and `decide()` are unchanged.
- Sites are `../recheck/sites_recheck.tsv`, unchanged: ok sites only, 2,779 in Set A and 30,943 in Set C.
- Python is `../recheck/.venv`. Nothing was rebuilt or installed.
- Probes are exact keep/lack probes (20 nt + site run + 20 nt), matched with Aho-Corasick on both orientations. There
  is one observation per spot, and a spot whose mates disagree is discordant and dropped.
- The modified `aggregate.py` was first run on a copy of `../recheck/test/split`. It reproduced the old test's
  `per_site_counts.tsv` byte for byte, and its `sets` and `per_run`.

Command (run_in_background; log in `counts/run.log`):
```
cd recheck_S
../recheck/.venv/bin/python -I units.py 2500000 | xargs -P 4 -L 1 ../recheck/.venv/bin/python -I count_run.py >> counts/run.log 2>&1
../recheck/.venv/bin/python -I aggregate.py              > aggregate.log                     # summary.json, per_site_counts.tsv
../recheck/.venv/bin/python -I aggregate.py --plus-old   > aggregate_plus_SRR37726778.log    # *_plus_SRR37726778.*
../recheck/.venv/bin/python -I report_S.py               > report_S.log
```
- **Units.** There were 64 units of 2,500,000 spots, 8 per run, at most 4 at a time.
  - The stream ran from 09:36:23 to 10:37:26. xargs and all 64 fastq-dumps exited 0.
  - Each unit held 2,500,000 spots and 5,000,000 reads (all paired). Its first and last spot ids equal its bounds.
  - A unit took 70–2,389 s.
- **Exact-spot check: passed for all 8 runs.** The 8 SRR37726778 units from `../recheck/counts/` also passed under
  `--plus-old`.
- **Memory.** I checked once, 30 s in, with `footprint`: about 233 MB per `count_run.py` and 0.29–0.85 GB per
  fastq-dump. That is well under 8 GB.
- I killed no process. Only counts are stored.

## Results: pooled over the 8 runs (the decision)

Powered: Set A has 6,731 observations (≥ 200) at 251 sites (≥ 30), and Set C has 14,577 (≥ 1,000). The keep condition
holds: hi_A = 0.0018 < 0.10.

| set | ok sites | sites observed | observations | keeps | lacks | discordant | D | 95% bootstrap (10,000, seed 1) |
|---|---|---|---|---|---|---|---|---|
| A | 2,779 | 251 | 6,731 | 6,728 | 3 | 0 | **0.00045** | [0, 0.0018] |
| C | 30,943 | 2,084 | 14,577 | 14,577 | 0 | 0 | **0** | [0, 0] |

**Extra, not in the rule: exact bounds** (Clopper–Pearson, one-sided 97.5% upper bound):

| | per observation | per observed site (as trials) |
|---|---|---|
| Set A | 3 of 6,731: 0.0013 | 1 of 251 sites with a lack: 0.022 |
| Set C | 0 of 14,577: 0.00025 | 0 of 2,084: 0.0018 |

**Other breakdowns:**
- Set C walks with a CGG: 1,942 observations at 225 of 2,904 sites, 0 lacks.
- Set C walks without a CGG: 12,635 observations at 1,859 sites, 0 lacks.
- Set A by run length:
  - 1 nt: 5,492 observations at 215 sites, with all 3 lacks;
  - 2 nt: 980 at 33 sites, 0 lacks;
  - ≥ 3 nt: 259 at 3 sites, 0 lacks.
- Set C at runs of ≥ 3 nt: 503 observations at 31 sites, 0 lacks.
- Set C by codon: CGA 2,401, CGC 2,219, CGT 9,957, all with 0 lacks.
- Sites with ≥ 10 observations: 70 in Set A and 223 in Set C, none with a lack.

## D per run and per study, and read-1 orientation

Read-1 sense share is the share of observations whose read 1 runs along the gene's sense strand, from `keeps_m1plus`
and `lacks_m1plus`. It is computed as in the old RECHECK.md. In a dUTP library it is near 0; for DNA or an unstranded
library it is near 0.5.

| run | study | Set A obs (sites) | A lacks | D_A [95%] | Set C obs (sites) | D_C | read-1 sense, A | read-1 sense, C |
|---|---|---|---|---|---|---|---|---|
| SRR37726779 | PRJNA1441222 | 2,528 (106) | 0 | 0 [0, 0] | 3,644 (942) | 0 | 3/2,528 = 0.0012 | 97/3,644 = 0.027 |
| SRR37726777 | PRJNA1441222 | 1,146 (82) | 1 | 0.0009 [0, 0.0035] | 2,459 (653) | 0 | 3/1,146 = 0.0026 | 32/2,459 = 0.013 |
| SRR37726776 | PRJNA1441222 | 395 (46) | 1 | 0.0025 [0, 0.0097] | 1,327 (389) | 0 | 0/395 = 0 | 21/1,327 = 0.016 |
| SRR37726773 | PRJNA1441222 | 760 (61) | 0 | 0 [0, 0] | 1,128 (388) | 0 | 6/760 = 0.0079 | 59/1,128 = 0.052 |
| SRR37726775 | PRJNA1441222 | 583 (73) | 0 | 0 [0, 0] | 1,000 (369) | 0 | 0/583 = 0 | 44/1,000 = 0.044 |
| SRR29669627 | PRJNA1129603 | 667 (54) | 1 | 0.0015 [0, 0.0076] | 2,086 (271) | 0 | 3/667 = 0.0045 | 69/2,086 = 0.033 |
| SRR29669633 | PRJNA1129603 | 314 (53) | 0 | 0 [0, 0] | 1,222 (303) | 0 | 7/314 = 0.022 | 119/1,222 = 0.097 |
| SRR29669628 | PRJNA1129603 | 338 (43) | 0 | 0 [0, 0] | 1,711 (219) | 0 | 4/338 = 0.012 | 42/1,711 = 0.025 |
| all 8 | | 6,731 (251) | 3 | 0.00045 [0, 0.0018] | 14,577 (2,084) | 0 | 26/6,731 = 0.0039 | 483/14,577 = 0.033 |
| SRR37726778 (old counts; plus line only) | PRJNA1441222 | 1,581 (82) | 0 | 0 [0, 0] | 4,468 (670) | 0 | 15/1,581 = 0.0095 | 0.0096 |

Every per-run D_C interval is [0, 0]. The SRR37726778 numbers are identical to the old RECHECK.md.

| study | runs | Set A obs (sites) | A lacks | D_A [95%] | Set C obs (sites) | D_C |
|---|---|---|---|---|---|---|
| PRJNA1441222 | 5 new | 5,412 (165) | 2 | 0.00037 [0, 0.0015] | 9,558 (1,563) | 0 [0, 0] |
| PRJNA1129603 | 3 | 1,319 (106) | 1 | 0.00076 [0, 0.0032] | 5,019 (628) | 0 [0, 0] |
| PRJNA1441222 with SRR37726778 (plus line) | 6 | 6,993 (188) | 2 | 0.00029 [0, 0.0011] | 14,026 (1,731) | 0 [0, 0] |

All 8 runs are stranded with read 1 antisense, as in dUTP libraries. Set A's read-1 sense share is 0–0.022 per run, and
Set C's is 0.013–0.097.

## Second pooled line: 8 runs + old SRR37726778 (reported only)

| set | sites observed | observations | keeps | lacks | D | 95% bootstrap |
|---|---|---|---|---|---|---|
| A | 271 | 8,312 | 8,309 | 3 | 0.00036 | [0, 0.0014] |
| C | 2,250 | 19,045 | 19,045 | 0 | 0 | [0, 0] |

- The test is powered, and the decision is "keep".
- Exact bounds:
  - Set A: 3 of 8,312 gives 0.00105; at site level, 1 of 271 sites gives 0.020.
  - Set C: 0 of 19,045 gives 0.00019.

## Set A by genome (8 runs)

Set A was observed in 17 of 22 primary genomes, and Set C in all 22.

| genome | Set A obs | sites | share of A | lacks | Set C obs |
|---|---|---|---|---|---|
| MGYG000302241 | 3,964 | 40 | 58.9% | 0 | 2,031 |
| MGYG000302328 | 1,380 | 56 | 20.5% | 0 | 4,296 |
| MGYG000301895 | 942 | 40 | 14.0% | 0 | 4,083 |
| MGYG000299987 | 202 | 25 | 3.0% | 0 | 594 |
| MGYG000301353 | 76 | 32 | 1.1% | 0 | 1,217 |
| MGYG000301495 | 49 | 7 | 0.7% | 0 | 31 |
| MGYG000302837 | 32 | 12 | 0.5% | 0 | 1,686 |
| MGYG000299934 | 7 | 3 | 0.1% | **3** | 30 |
| 9 others | 79 | | 1.2% | 0 | |

The largest Set A sites:
- MGYG000302241_5:5610 (−), with 1,715 keeps: 25% of all Set A observations;
- MGYG000301895_50:11189 (−), with 367;
- MGYG000302241_3:87744 (+), with 339;
- MGYG000302241_2:108007 (−), with 277.

With SRR37726778 added, MGYG000302328 rises to 2,720 observations (32.7%) and MGYG000302241 falls to 50.1%.

## Things that look wrong or need care

1. **The 3 lacks are at one site, which has no keeps.**
   - The site is `MGYG000299934:MGYG000299934_7:8068:+` (CGG, run length 1, CDS MGYG000299934_00105 = ffh_1).
     ffh_1 (7,732–8,121) and ffh_2 (8,264–8,992) are the two halves of a split Ffh, so this CGG is a skip site in a
     split gene.
   - Its counts are 0 keeps and 3 lacks. Each lack is one spot, from three different runs in both studies:
     SRR37726777, SRR37726776 and SRR29669627.
   - In all 3, read 1 is antisense to the gene (`lacks_m1plus` = 0 on a + gene). That is the same RNA-like orientation
     as the keeps.
   - The 10 old main runs had 0 observations at this site.
   - The lack sequence (`...GGTTCA|GGAAAAACAACAACAGCAGG`) read in the standard code gives ...S G K T T T A. That is
     exactly the protein that the skip reading of the reference gives at this site.
   - So these reads look like a strain or relative whose DNA has a continuous ffh without the C. The rule's own caveat
     ("a genomic deletion of the C in the sampled strains would look the same") applies here, at one site.
   - The lack probe occurs nowhere in the 24 catalogue genomes. Uniqueness was only checked against them, so reads from
     an uncatalogued strain can still match.
   - This does not change the decision: D_A = 0.00045, hi_A = 0.0018. If the DNA check can look at one site, this is
     the one.
2. **Concentration.**
   - One site (MGYG000302241_5:5610) gives 25% of Set A observations, and one genome (MGYG000302241) gives 59%.
   - The genome mix differs from the old main set. MGYG000301895 had 0 Set A observations in the 10 main runs, but
     has 942 here. MGYG000301353 was the second-largest there, and has 76 here.
   - The bootstrap over sites accounts for clustering, and even so hi_A is 0.0018.
3. **Set C has more read-1-sense observations than Set A** (0.033 against 0.0039). They are not spread evenly, as
   DNA would be:
   - 20 Set C sites carry 237 of the 483.
   - They are mainly in MGYG000302837 (219 of 1,686 = 0.13; e.g. MGYG000302837_01154, 37 of 40 sense) and
     MGYG000299987 (80 of 594 = 0.135).
   - This looks like antisense or overlapping opposite-strand transcription at a few loci.
   - All of these are keeps, so D_C is not affected.
   - SRR29669633's Set C share (0.097) is the highest. Its Set A share is 0.022.
4. **Degenerate intervals.** With 0 lacks, D_C's bootstrap interval is [0, 0]. With all lacks at one site, D_A's upper
   end depends on how often that site is drawn. The exact bounds above are given for that reason; they are not part of
   the rule.
5. **New samples.** The ENA table gives only a generic sample title for the PRJNA1441222 runs ("Metagenome or
   environmental sample from gut metagenome"). I could not check that SRR37726773–79 are animals or samples other
   than SRR37726778. They are separate runs from the same NovaSeq submission. The PRJNA1129603 runs are YK1-1, YK1-2
   and YK1-6.
6. **Head of each run.** Each run has 33–106 M spots, and only the first 20 M were counted, as Follow-up S asks.
7. **Exact-match sensitivity.** Observations are fewer than an aligner's, because a read must cover the 41–45 nt
   window with no error. This affects keeps and lacks alike (old RECHECK.md, item 7).

## Files

- `CHANGES.md`: the diff summary of the three copied scripts.
- `count_run.py`, `units.py`, `aggregate.py`: the copied scripts.
- `report_S.py`: new; it prints `report_S.log` from the two summaries.
- `counts/`: 64 units, each with `.tsv`, `.summary` and `.err`, plus `run.log`.
- `summary.json`, `per_site_counts.tsv`, `aggregate.log`: the 8-run decision.
- `summary_plus_SRR37726778.json`, `per_site_counts_plus_SRR37726778.tsv`, `aggregate_plus_SRR37726778.log`: the
  second pooled line.
- `report_S.log`: all the numbers above, including per-genome and sites with ≥ 10 observations.
