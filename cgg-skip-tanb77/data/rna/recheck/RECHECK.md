# RNA CGG test: independent re-check (alignment-free)

**Status (2026-11-12 09:13): DONE under Amendments 3 and 4. Decision: transcripts keep the C.**
- 10 runs from the rewritten `main_runs.txt`, first min(20 M, run) spots each, 198,467,619 spots in all, were streamed
  08:43–09:10. The exact-spot check passed for every run.
- Set A: 5,489 observations at 615 sites, 0 lacks, **D_A = 0 [0, 0]**. Set C: 51,787 observations at 10,675 sites, 0
  lacks, **D_C = 0 [0, 0]**. Powered; hi_A = 0 < 0.10, so the rule's decision is "keep" (section 6).
- With 0 lacks the percentile bootstrap interval collapses to [0, 0]; section 6 adds exact binomial bounds, which are
  not part of the rule.
- Earlier: at 08:07 the coordinator stopped the old run list (Amendment 2: the old screen score measured rRNA); its
  counts are sealed and unused (section 4). Amendments 3 and 4 are implemented here in my own code.

Method, per `../recheck_brief.md`: I built my own sites, then matched probes exactly, with no aligner.
- Keep probe = 20 nt left + site run + 20 nt right. Lack probe = the same with the run one base shorter.
- A read in either orientation is matched against the probes. One observation per spot. A spot whose mates disagree is
  counted as discordant and dropped.
- A site is dropped if its probes are not unique across the 24 genomes, on both strands.

I did not read `src/rna_cgg.py`, `src/rna_stream.py`, anything in `data/rna_cgg/`, or any `*.obs.tsv`.

## 1. Sites compared with `sites.tsv`

`build_sites.py` does the site walk. It imports `walk`, `gff`, `rc` and `read_fasta` from `src/proteomics_db.py`.
- It walks every catalogue CDS in skip mode. The walk goes past the CDS end, up to the contig end.
- It takes Set A (the C of in-frame CGG) and Set C (the C of in-frame CGA/CGC/CGT, in the walk's own frame) from the 22
  primary genomes.
- Positions are pooled per genome as the union over all walks.
- A site is dropped if its C lies within 30 nt of a contig end.
- The site run is the maximal homopolymer on the plus strand that contains the C.

**Amendment 3** is implemented in `unique_probes.py`, independently, from the rule's definition; `sites.tsv` is not used
to filter.
- A site is dropped (status `ncrna`) if its window [run_start0 − 20, run_end0 + 20) overlaps any GFF feature of type
  `ncRNA`. The window is 0-based half-open on the plus strand; GFF start..end is taken as [start − 1, end).
- The catalogue GFFs hold 344 `ncRNA` features: tRNA, tmRNA, RNase P, rRNA, riboswitches and others.
- **It drops 5 Set C sites and 0 Set A sites.** All 5 are in MGYG000300925_01647, on contig MGYG000300925_19, at C
  positions 1281, 1494, 1512, 1566 and 1605 (1-based, + strand). All 5 come from walks with no CGG. Four of them had
  passed the uniqueness filter; the fifth had already failed it.

`compare_sites.py` compares the result with `sites.tsv` (log: `test/compare_sites.log`). Note that `sites.tsv` positions
are 0-based.

| set | rebuilt | sites.tsv (Amendment 3) | shared | only rebuilt | only sites.tsv | codon / run / CDS differences |
|---|---|---|---|---|---|---|
| A | 2,871 | 2,871 | 2,871 | 0 | 0 | 0 |
| C | 32,496 | 32,496 | 32,496 | 0 | 0 | 0 |
| A0 (not used) | 238 | 238 | 238 | 0 | 0 | 0 |

**The site sets are identical.** Before the ncRNA drop they are also identical to `sites_v1.tsv`: A 2,871, C 32,501.
- Set C's `walk_has_cgg` split is 3,066 sites with a CGG in the walk and 29,430 without. This is the same as the rule's
  Amendment 3.
- Run lengths in Set A are 2,549 / 295 / 26 / 1 for lengths 1 / 2 / 3 / 4.
- I did not build Set B, since the brief asks for Sets A and C only.

## 2. Probe uniqueness (`unique_probes.py` → `sites_recheck.tsv`)

Every keep probe and lack probe was searched exactly (Aho-Corasick) in all contigs of the 24 genomes, on both strands. A
site is kept ("ok") only if four things hold:
- it is not dropped by the ncRNA check;
- its keep probe occurs at exactly one place, its own;
- its lack probe occurs nowhere;
- neither probe contains a non-ACGT base.

| set | candidates | ncRNA (Amendment 3) | keep probe at ≥ 2 places | lack probe present elsewhere | **ok (counted)** |
|---|---|---|---|---|---|
| A | 2,871 | 0 | 90 | 2 | **2,779** |
| C | 32,501 | 5 | 1,553 | 0 | **30,943** |

Where the other copies lie:
- **A:** 81 of the 90 keep_multi sites have their other copies in another primary genome, 8 in the same genome only and
  1 in a non-primary genome only. Both lack_present sites have the lack version in another primary genome.
- **C:** most keep_multi sites have their other copies in another primary genome. Before Amendment 3 the split was 1,472
  in another primary genome, 70 in the same genome only and 12 in a non-primary genome only, out of 1,554.

Of the ok Set C sites, 2,904 come from a walk with a CGG and 28,039 do not.

The original filters reads by MAPQ ≥ 20 rather than filtering sites, so its usable sites will differ from mine by about
these numbers. All 22 genomes keep ok sites in both sets:
- Set A: from 31 per genome (MGYG000302139) to 278 (MGYG000302432);
- Set C: from 452 to 6,142 (MGYG000301353).

## 3. Counter and tests

`count_run.py <run> --from A --to B` streams `fastq-dump --stdout --split-spot --skip-technical -N A -X B <run>` and
matches each read against the probes and their reverse complements. `--max-spots N` is the same as `--from 1 --to N`.
- Spots are grouped by read name. In --split-spot output the mates of a spot are consecutive.
- For each spot and site it records one of keep, lack or discordant. It also records whether read 1 runs along the
  plus strand, which is used for the strandedness check.
- It writes only counts, to `counts/<run>.<A>-<B>.tsv`, plus a `.summary` with run, from, to, spots, reads, the first
  and last spot numbers seen in the read names, the exit code and seconds. The .tsv is written only if fastq-dump
  exits 0.

`units.py [CHUNK]` reads `../main_runs.txt` ("dir 1 SPOTS 2 RUN").
- It cuts spots 1..SPOTS of each run into ranges of at most CHUNK spots (default 2,500,000).
- It prints one line per range not yet done, so the same command resumes after an interruption.

`aggregate.py` takes its runs from `../main_runs.txt`.
- **Exact-spot check (Amendment 4).** For each run:
  - every unit exited 0;
  - the ranges start at 1, have no gap or overlap, and end at SPOTS;
  - each range holds exactly B − A + 1 spots, its first read is spot A and its last is spot B;
  - the summed spots equal SPOTS.

  Any failure stops the aggregate.
- It pools the counts, computes D = lacks/(keeps+lacks) per set, and takes a 95% percentile bootstrap over sites with ≥ 1
  observation (10,000 resamples, seed 1).
- It applies the rule's Decision. `decide()` was unit-tested on 9 hand-made cases, all correct.
- It writes `per_site_counts.tsv` (every ok site, zeros included) and `summary.json`.
- `summary.json` also holds extras that are reported but not used in the decision:
  - D per run;
  - Set C split by whether the walk has a CGG (`C_walks_with_cgg`);
  - read-1 strandedness per run;
  - D by run length, by Set C codon and by genome;
  - Set A sites with ≥ 10 observations.

**Known-answer test** (`known_answer.py`, log `test/known_answer.log`; rerun on the Amendment 3 sites): **PASS**.
- Setup: 12,000 simulated spots (22,000 reads) at 1,000 random ok sites from each set. They include:
  - paired spots from both strands, keep and lack versions;
  - single-end reads;
  - discordant spots;
  - substitutions at 0.5% per base.
- Expected counts were computed from coordinates, not from strings.
- Result: 0 of 3,584 touched sites differed.

| set | keeps, lacks, discordant, keeps read-1-plus, lacks read-1-plus: expected | got |
|---|---|---|
| A | 2,155, 1,618, 566, 1,073, 817 | 2,155, 1,618, 566, 1,073, 817 |
| C | 4,891, 1,582, 598, 2,417, 797 | 4,891, 1,582, 598, 2,417, 797 |

**Split test** (Amendment 4 ranges; in the same log): **PASS**.
- The same 12,000 spots were counted as 3 ranges (1–4,000, 4,001–8,000 and 8,001–12,000). `aggregate.py` joined them to
  per-site counts identical to the unsplit run.
- Two broken copies were refused (aggregate exit 1): one with the middle range missing, one with a summary claiming one
  spot fewer.

**End-to-end test** (`test/e2e/`, before Amendments 3–4): `--max-spots 5000` on SRR11426316 and SRR31853479.
- Each run streamed 5,000 spots (10,000 reads) and fastq-dump exited 0, in 5 s per run.
- There were 0 observations, as expected for a 5,000-spot head. The aggregate ran through.
- The `-N`/`-X` range path was first streamed in the main run. There, every one of the 80 units began and ended
  exactly on its bounds, as the aggregate's per-range check confirmed.

**Memory:** peak `phys_footprint` was 233 MB per `count_run.py` (from `footprint -p`); maximum RSS on 1 M spots was
255 MB. With 6 at a time that is about 2 GB, with fastq-dump on top.

**Speed:**
- The matcher alone handles about 90,000 spots/s: 1 M local synthetic 2 × 150-nt spots took 11 s, including the 4 s
  automaton build.
- The stream is the limit. On the stopped old list (4 at a time) it was about 6,000 spots/s per fastq-dump. On the new
  main runs, 6 at a time while the coordinator's 8 streams also ran, it was about 24,100 spots/s per stream over the
  79 full units. Each 2.5 M-spot unit took 82–251 s.

## 4. The old run list (stopped, sealed, not reported)

Before the stop, 8 of the 10 old runs had streamed whole, with spot counts equal to the old `main_runs.txt`.
SRR11426327 and SRR11426338 were killed mid-stream; I killed only my own saved PIDs.
- Their per-site counts, and a flank diagnostic of the first 1 M spots of SRR11426316, were moved to `sealed_oldlist/`.
- They are not summarized here and not used by `aggregate.py`.
- I report no keep/lack split or D from them.

What I report from them is only abundance with no site and no class. `kmer_census.py` (output in `diag/`) classes each
read by exact 31-mers of the genomes' annotated rRNA, of the primary genomes' CDS (sampled every 128 nt) and of other
ncRNA:

| run | reads classed | rRNA only | CDS only | CDS + rRNA | none |
|---|---|---|---|---|---|
| SRR11426316 (first 1 M spots) | 2,000,000 | 646,745 | 11,458 | 4,444 | 1,337,347 |
| SRR10686660 (whole run) | 144,400 | 82,270 | 0 | 0 | 62,130 |
| SRR31853479 (whole run) | 164,606 | 64,490 | 0 | 0 | 100,116 |

This agrees with Amendment 2: the old screen picked rRNA-rich runs with little or no TANB77 mRNA.

## 5. How it was run

The first command ran with run_in_background (log: `counts/run.log`):

```
cd results/rna_cgg/recheck
mkdir -p counts
.venv/bin/python -I units.py 2500000 | xargs -P 6 -L 1 .venv/bin/python -I count_run.py >> counts/run.log 2>&1
.venv/bin/python -I aggregate.py --totals-only   # blinded: totals and power only
.venv/bin/python -I aggregate.py                 # powered, so the full run: summary.json, per_site_counts.tsv
.venv/bin/python -I report.py                    # section 6 numbers (report.log)
```

- Each run is cut into ranges of at most 2.5 M spots, 8 per 20 M-spot run, with at most 6 count processes (and 6
  fastq-dumps) at a time. Together the ranges are exactly `fastq-dump -N 1 -X SPOTS`.
- If the first command stops early (xargs exit ≠ 0, or the memory reaper), rerun it unchanged. `units.py` leaves out the
  ranges already done.
- If Amendment 2's "10 more runs" step had applied, the extra runs would be more lines in the same file, picked up by
  the same commands. It did not apply: Set A had more than 200 observations.
- No run of the new list is in `sealed_oldlist/`.
- **Wall time.** I had expected 1.5–2.5 h from the old streams' rate. The actual time was 26 min: 80 units from
  08:43:20 to 09:09:18, xargs exit 0, every unit exit 0. The ERR15873922 tail re-count took a further 34 s (item 12 in
  section 7). The aggregate took 3 s.
- **ERR15873922.** The list first gave 20,000,000 spots, but the run has 18,467,619. Its last unit, 17,500,001–
  20,000,000, streamed 967,619 spots (fastq-dump exit 0), so it would have failed the exact-spot check. The coordinator
  corrected `main_runs.txt` to 18,467,619 at 08:53. I moved that unit to `counts_superseded/` and re-counted 17,500,001–
  18,467,619. The two count files are byte-identical.

## 6. Results

**Runs** (`main_runs.txt` as rewritten, with ERR15873922 corrected): ERR15738679, ERR15873938, ERR15873922,
ERR15873927, ERR15873933, ERR15873934, ERR15873936, ERR15873944, ERR15873926 and SRR37726778.
- Each run's first 20,000,000 spots were streamed, except ERR15873922, which was streamed whole (18,467,619 spots).
- Every spot had 2 reads.
- **Exact-spot check: passed for all 10 runs** (80 units). In each unit the first and last spot numbers in the read
  names equal its bounds, the ranges are contiguous from 1, and the sums equal SPOTS.

**Decision: transcripts keep the C.** The test is powered: Set A has 5,489 observations (≥ 200) at 615 sites (≥ 30), and
Set C has 51,787 (≥ 1,000). The keep condition is hi_A < 0.10, and hi_A = 0.

| set | ok sites | sites observed | observations | keeps | lacks | discordant | D | 95% bootstrap (10,000) |
|---|---|---|---|---|---|---|---|---|
| A | 2,779 | 615 | 5,489 | 5,489 | 0 | 0 | **0** | [0, 0] |
| C | 30,943 | 10,675 | 51,787 | 51,787 | 0 | 0 | **0** | [0, 0] |

No observation in any run, set, codon or run-length class was a lack. With no lacks, every bootstrap resample gives
D = 0, so the rule's interval is [0, 0].

**Extra, not in the rule: exact binomial bounds** (Clopper–Pearson 97.5% one-sided upper bound for 0 events in n):

| | n observations | bound per observation | n observed sites | bound on the share of sites with any lack |
|---|---|---|---|---|
| Set A | 5,489 | 0.00067 | 615 | 0.0060 |
| Set C | 51,787 | 0.00007 | 10,675 | 0.0003 |
| SRR37726778, Set A | 1,581 | 0.0023 | 82 | 0.044 |

The site-level bound treats each observed site as one trial, which allows for clustering by site. All of these bounds
are far below the rule's 0.10.

### Reported extras (Amendment 3)

**D per run and strandedness per run.** Read-1 sense share is the share of observations whose read 1 runs along the
gene's sense strand.

| run | spots | Set A obs (sites) | D_A | Set C obs (sites) | D_C | read-1 sense share, A | read-1 sense share, C |
|---|---|---|---|---|---|---|---|
| ERR15738679 | 20,000,000 | 337 (60) | 0 | 6,707 (1,093) | 0 | 0.226 | 0.226 |
| ERR15873938 | 20,000,000 | 522 (208) | 0 | 4,828 (2,490) | 0 | 0.500 | 0.492 |
| ERR15873922 | 18,467,619 | 352 (123) | 0 | 4,794 (2,246) | 0 | 0.551 | 0.491 |
| ERR15873927 | 20,000,000 | 656 (167) | 0 | 8,947 (4,277) | 0 | 0.509 | 0.497 |
| ERR15873933 | 20,000,000 | 307 (131) | 0 | 4,458 (2,117) | 0 | 0.456 | 0.502 |
| ERR15873934 | 20,000,000 | 251 (125) | 0 | 3,990 (2,043) | 0 | 0.434 | 0.483 |
| ERR15873936 | 20,000,000 | 234 (113) | 0 | 3,881 (1,840) | 0 | 0.577 | 0.484 |
| ERR15873944 | 20,000,000 | 674 (184) | 0 | 5,364 (2,522) | 0 | 0.518 | 0.511 |
| ERR15873926 | 20,000,000 | 575 (166) | 0 | 4,350 (1,755) | 0 | 0.522 | 0.501 |
| SRR37726778 | 20,000,000 | 1,581 (82) | 0 | 4,468 (670) | 0 | 0.0095 | 0.0096 |
| all | 198,467,619 | 5,489 (615) | 0 | 51,787 (10,675) | 0 | 0.349 | 0.419 |

- Every per-run bootstrap interval is [0, 0].
- **Strandedness falls into three groups:**
  - Eight runs are near 0.5. These are unstranded libraries, or DNA; this measure cannot tell the two apart.
  - ERR15738679 is at 0.23, in both sets alike.
  - SRR37726778 is stranded, with read 1 antisense (0.0095), as in dUTP libraries.
- **SRR37726778 alone** has 1,581 Set A observations at 82 sites and 4,468 Set C observations, with 0 lacks. This one
  clearly RNA-derived run meets the power thresholds by itself and gives the same decision. Its Set A observations are
  concentrated, at about 19 per site.

**D_C for walks with a CGG** (`C_walks_with_cgg`): 6,046 observations at 755 of 2,904 sites, 0 lacks, D = 0 [0, 0].
Walks without a CGG: 45,741 observations at 9,920 of 28,039 sites, 0 lacks, D = 0 [0, 0].

**Set A sites with ≥ 10 observations:** 119 sites, all with 0 lacks.
- The largest are MGYG000302328_5:842 (−) with 335 keeps, MGYG000301353_101:2743 (−) with 162,
  MGYG000302241_5:5610 (−) with 162, and MGYG000302328_12:34372 (+) with 150.
- The full list is `sites_ge10_A` in `summary.json`, and also in `report.log`.
- Set C: 1,173 sites have ≥ 10 observations, and none has a lack.

**Other breakdowns** (all with 0 lacks):
- Set A by site run length: 1 nt, 4,618 observations at 535 sites; 2 nt, 855 at 76; ≥ 3 nt, 16 at 4.
- Set C by site run length: 1 nt, 47,112; 2 nt, 4,162; ≥ 3 nt, 513 at 98 sites.
- Set C by codon: CGA 13,168; CGC 6,684; CGT 31,935.
- Per genome:
  - Set A was observed in 20 of 22 genomes; MGYG000301895 and MGYG000301929 had none.
  - Four genomes give 78% of Set A observations: MGYG000302328 (1,345), MGYG000301353 (1,284), MGYG000302241 (876) and
    MGYG000302837 (761).
  - Set C was observed in all 22.

## 7. Things that look wrong or need care

Items 1, 2 and 5 are now handled in the rule by Amendment 3.

1. **Set C is mostly not "the same genes"** (wording fixed in Amendment 3). 28,039 of the 30,943 ok Set C sites come
   from walks that meet no CGG at all; only 2,904 share a walk with a CGG. D_C therefore measures the genomes' background
   rate. D_C for walks with a CGG is reported separately (`C_walks_with_cgg`).
2. **CDS called inside rRNA** (sites dropped by Amendment 3). Prodigal calls MGYG000300925_01647 (1260–1646 on
   MGYG000300925_19) inside the 23S rRNA (5–2872), and reads there are rRNA. In the stopped old streams, most Set C
   observations were at these sites. That is totals by location only; no keep/lack split was looked at.
3. **Few runs from one study.** Seven of the ten old runs came from one study. In the new list, eight runs are
   ERR15873922–ERR15873944, consecutive accessions that probably come from one submission. D is 0 in every run, so here
   the pooled D does not hide disagreement between runs.
4. **Unique probes versus MAPQ.** 2.9% of Set A sites and 4.5% of Set C sites are lost in my build because a probe is
   shared with another primary genome. Two Set A sites have a lack version present in another primary genome. Reads
   from that genome would be called lacks by exact matching, so I drop those sites. The original's MAPQ ≥ 20 filter
   handles this per read, so its observation counts can be a little higher than mine.
5. **DNA in the library** (strandedness is a reported extra under Amendment 3). Genomic DNA reads give keeps at both
   sets and would pull D toward 0. In a stranded library, the share of observations whose read 1 is sense to the gene
   should sit near 0 or 1 per run. A share near 0.5 points to DNA reads, or to an unstranded protocol.
6. **Genomic deletion.** A "lacks" result can also be a strain whose genome lacks the C, as the rule says. Per-genome
   and per-site D (`per_genome`, `sites_ge10_A`) can show whether lacks come from a few sites or from many.
7. **Exact-match sensitivity.** A read must cover the whole 41–45 nt window with no error. Errors and read ends lose
   keeps and lacks almost equally; the lack probe is 1 nt shorter, which is negligible. So D is unbiased but
   observations are fewer than with an aligner.
8. **First 20 M spots (Amendment 4).** The first spots of a run are not a random sample, for example by flow-cell tile
   or lane order. That should not change keep or lack at a site, but it is one more reason to read D per run.
9. **No lacks at all, so the bootstrap interval is degenerate.** The rule's percentile bootstrap gives [0, 0], which
   says nothing about uncertainty. The decision does not depend on this: the exact bounds in section 6 are ≤ 0.006 for
   Set A even at site level.
   - Set C also has 0 lacks, including 513 observations at runs of ≥ 3 nt. So exact matching picked up no
     sequencing-error deletions at these loci either.
   - The lack path itself is covered by the known-answer test, which uses the same probe table and the same automaton;
     only the input source differs.
   - The real reads gave 57,276 keep observations, so the read text is plain upper-case ACGT, as the matcher needs.
   - Discordant = 0 follows from this: a discordant spot needs one mate to match a lack probe.
10. **Strandedness: most runs could be unstranded or DNA.** Eight runs sit at 0.43–0.58 read-1 sense, which genomic DNA
    would also give, and DNA keeps the C by construction.
    - SRR37726778 is stranded (0.0095). By itself it gives 1,581 Set A observations at 82 sites with 0 lacks, enough for
      the rule's power thresholds.
    - So the "keep" result does not rest on runs that might be DNA. The SRR37726778 subset is my check, not the rule's.
11. **Concentration.** 78% of Set A observations come from 4 of the 22 genomes, and 119 Set A sites have ≥ 10
    observations each. The answer is about these well-expressed genes and genomes, not evenly about all 2,779 sites.
    Even so, 615 Set A sites were observed, and none had a lack.
12. **ERR15873922 run size.** The first list gave 20,000,000 spots for a run of 18,467,619. The exact-spot check would
    have refused it; the coordinator corrected the list. The superseded unit is kept in `counts_superseded/`, and its
    counts are byte-identical to the re-count.

## Files

- `build_sites.py` → `candidates.tsv`: all walk candidates, including A0 and C0 from the non-primary genomes.
- `unique_probes.py` → `sites_recheck.tsv`: Amendment 3 ncRNA check and probe uniqueness; status per site, and only
  `ok` sites are counted.
- `compare_sites.py`: compares the sites with `sites.tsv` and `sites_v1.tsv`.
- `units.py`, `count_run.py`, `aggregate.py`, `report.py`: range planning, the counter, the pooled analysis and the
  section 6 printout.
- `counts/`: per-unit counts, summaries and fastq-dump stderr; `counts/run.log` is the stream log.
- `counts_superseded/`: the ERR15873922 17,500,001–20,000,000 unit (item 12).
- `totals.json` (blinded totals), `summary.json`, `per_site_counts.tsv` (all 33,722 ok sites, zeros included),
  `aggregate.log`, `report.log`.
- `known_answer.py`, `test/`: the known-answer, split and end-to-end tests.
- `kmer_census.py`, `flank_diag.py`, `diag/`: abundance diagnostics, with no site classes in `diag/`.
- `sealed_oldlist/`: counts from the old list; not to be read for the test.
- `.venv`: uv venv with pyahocorasick and numpy.
