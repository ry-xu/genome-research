#!/bin/zsh
# usage: run_codetta.sh <codetta_dir> <genome.fna>...  (runs sequentially)
C=$1; shift
for g in "$@"; do
  echo "START $g $(date)"
  HMMER_NCPU=4 uv run --no-project --python 3.12 --with 'numpy<2' --with 'scipy<1.14' python -E -s $C/codetta.py $g -m > $g.codetta.log 2>&1
  echo "END $g rc=$? $(date)"
done
