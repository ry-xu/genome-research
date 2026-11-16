# Changes from `../recheck/` (Follow-up S copy)

`count_run.py`, `units.py` and `aggregate.py` were copied from `../recheck/` at 09:34 on 2026-11-12. These are all the
changes (`diff ../recheck/X X`). The counting, matching, exact-spot check, bootstrap and `decide()` are unchanged.
Sites still come from `../recheck/sites_recheck.tsv`, and Python is still `../recheck/.venv` (nothing rebuilt or
installed).

## count_run.py
- Docstring: one added line saying this is the Follow-up S copy.
- Added `S = ".../recheck_S/"`. `R` (still `.../recheck/`) is used only for `sites_recheck.tsv`.
- The default `--out` is `S + "counts"` (was `R + "counts"`).

## units.py
- Docstring: one added line, and `../main_runs.txt` changed to `../followS_runs.txt`.
- `R` is `.../recheck_S/` (was `.../recheck/`). It is used for the run list `R + "../followS_runs.txt"` (was
  `../main_runs.txt`) and for the done-check in `R + "counts/"`.

## aggregate.py
- Docstring: added lines on the Follow-up S copy and on `--plus-old`.
- Added constants `S = ".../recheck_S/"`, `OLD = "SRR37726778"` and
  `STUDY = {"SRR377267": "PRJNA1441222", "SRR296696": "PRJNA1129603"}` (run-accession prefix to study).
  `R` stays `.../recheck/` and is used for `sites_recheck.tsv` and for the old SRR37726778 units.
- Default counts dir is `S + "counts/"` and default output dir is `S` (were `R + "counts/"` and `R`).
- Run list: `S + "../followS_runs.txt"` (was `R + "../main_runs.txt"`). Only the messages that named
  `main_runs.txt` were changed to match.
- **New `--plus-old` flag** (task item 4, reported only):
  - It adds the old re-check's `../recheck/counts/SRR37726778.*.summary` units.
  - It asserts that `recheck_S/counts/` holds no SRR37726778 unit.
  - It takes SRR37726778's expected SPOTS from its line in `../main_runs.txt` (20,000,000).
  - The exact-spot check applies to these units as to the others.
  - The outputs get a `_plus_SRR37726778` suffix: `summary_plus_SRR37726778.json` and
    `per_site_counts_plus_SRR37726778.tsv`.
- **New extra `per_study`** (and `study_of_run`), reported only. Per-run counts are summed per study, then the same
  `setstats` is applied as for `per_run`.
  - It runs after `per_run`, on the same RNG stream. So the pooled A/C intervals (computed first, seed 1) are unchanged
    by it.
  - Runs with no listed prefix would go under "other".

## Check of the copy
- I ran `aggregate.py --counts <copy of ../recheck/test/split> --partial` in the scratchpad.
- `per_site_counts.tsv` was byte-identical to `../recheck/test/split/per_site_counts.tsv`.
- `sets` and `per_run` in `summary.json` were identical to the old test output.
- `per_study` matched `per_run` in counts; only the bootstrap draws differ, because it uses later draws from the RNG.

## New file (not a copy)
- `report_S.py`: prints `report_S.log` from `summary.json` and `summary_plus_SRR37726778.json`. It is like
  `../recheck/report.py`, and adds per study, per-genome sorting and Clopper–Pearson bounds when a set has 0 lacks.
