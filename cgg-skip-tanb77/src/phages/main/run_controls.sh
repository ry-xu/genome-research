#!/bin/bash
# N null world, P3 recall (f=0.5, 1.0) and N held-out on R; target bin sizes from T1 (run_ids.tsv only)
set -e
S="uv run --no-project --with numpy python -I src/score.py data/groups/run_ids.tsv"
$S results/main/ctrl11_plant0.tsv null --target T1 --nnull 200 --perm 1000 --seed 11 2>&1 | grep -v Warn > results/controls/null_T1.txt
$S results/main/ctrl11_plant0.tsv recall --target T1 --plantstats results/main/ctrl11_plant0.5.tsv --nnull 200 --perm 1000 --seed 12 2>&1 | grep -v Warn > results/controls/recall_f0.5_T1.txt
$S results/main/ctrl11_plant0.tsv recall --target T1 --plantstats results/main/ctrl11_plant1.0.tsv --nnull 50 --perm 1000 --seed 13 2>&1 | grep -v Warn > results/controls/recall_f1.0_T1.txt
$S results/main/ctrl11_plant0.tsv heldout --seed 14 2>&1 | grep -v Warn > results/controls/heldout.txt
