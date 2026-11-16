"""Stream reads of SRA runs (fastq-dump; nothing stored but hits) and map them to the 24 TANB77 pig-gut
representatives with minimap2 (-x sr). Adapted from another agent's work/tag-proteomics/src/rna_stream.py (read, not modified).
Keeps only primary alignments with MAPQ >= 20: <outroot>/<run>/part_<first spot>.sam, plus part_<first spot>.n =
reads streamed in that part.

usage: python -I src/rna_stream.py <outroot> <first_spot> <n_spots> <streams> <run> [run ...]
"""
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REF = Path(os.environ.get("RNA_REF", ROOT / "data" / "rna_cgg" / "ref" / "tanb77_24.mmi"))  # chicken: RNA_REF
MAPQ = 20  # only primary alignments with MAPQ >= 20 are kept (all that rna_cgg.py ever counts); the rest is rRNA
# multimappers, 99% of hits in screen runs


def stream(outroot, run, first, n, streams):
    out = Path(outroot) / run
    out.mkdir(parents=True, exist_ok=True)
    step = -(-n // streams)
    procs = []
    for a in range(first, first + n, step):
        b = min(a + step, first + n) - 1
        part = out / f"part_{a}"
        if (part.parent / f"{part.name}.n").exists():
            continue
        cmd = (f"set -o pipefail; fastq-dump --stdout --split-spot --skip-technical -N {a} -X {b} {run} 2>{part}.err"
               f" | tee >(awk 'END{{print NR/4}}' > {part}.tmp_n)"
               f" | minimap2 -x sr -a -t 1 --sam-hit-only {REF} - 2>/dev/null"
               f" | awk -F'\\t' '$1 ~ /^@/ || ($5 >= {MAPQ} && int($2/4)%2 == 0 && int($2/256)%2 == 0"
               f" && int($2/2048)%2 == 0)' > {part}.sam"
               f" && sleep 1 && mv {part}.tmp_n {part}.n")
        procs.append(subprocess.Popen(["bash", "-c", cmd]))
    bad = sum(p.wait() != 0 for p in procs)
    reads = sum(float(p.read_text() or 0) for p in out.glob("part_*.n"))
    print(run, "parts failed", bad, "reads streamed", int(reads), flush=True)


def main():
    outroot = sys.argv[1]
    first, n, streams = map(int, sys.argv[2:5])
    for run in sys.argv[5:]:
        stream(outroot, run, first, n, streams)


if __name__ == "__main__":
    main()
