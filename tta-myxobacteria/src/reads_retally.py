"""Re-tally check 3 from the kept SAMs, counting read pairs (one counted alignment per read name per site).

Audit 2 (2026-11-11 12:24, scratch/audit_reads_cat/) found that overlapping mates of one pair both counted at a
site: minimap2 runs with --frag=no, so mates align as separate records sharing a read name. src/reads_check.py
tally() now keeps one per name; this reruns that tally on the SAMs of genomes finished before the fix.
Sites: the genome's own data/reads/<acc>/sites.tsv (as written by the run); runs: those with a .done file;
preset from the run's report line (results/reads/<acc>.txt). Nothing is streamed again.

usage: uv run --no-project python -I src/reads_retally.py <acc> [<acc> ...]
       -> results/reads_pairs/<acc>.tsv and .txt (same format as results/reads/)
"""
import csv
import fcntl
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
import reads_check as rc  # noqa: E402

rc.OUTD = ROOT / "results" / "reads_pairs"


def main(acc):
    d = ROOT / "data" / "reads" / acc
    lock = open(d / ".lock", "w")
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        sys.exit(f"{acc}: a reads_check is running on this genome")
    head = (ROOT / "results" / "reads" / f"{acc}.txt").read_text().splitlines()[0]
    m = re.match(r"(\S+) (\w{3}); runs (\S+); minimap2 -x (\S+);", head)
    assert m and m.group(1) == acc, head
    x, runs, preset = m.group(2), m.group(3).split(","), m.group(4)
    sites = list(csv.DictReader(open(d / "sites.tsv"), delimiter="\t"))
    assert all(r == "" or (d / f"{r}.done").exists() for r in runs), f"{acc}: a reported run has no .done"
    sams = [p for run in runs for p in sorted(d.glob(f"{run}.part*.sam"))]
    assert sams, f"{acc}: no SAMs"
    counts = rc.tally(sams, sites, preset)
    rc.report(acc, x, sites, counts, runs, preset)


if __name__ == "__main__":
    for a in sys.argv[1:]:
        main(a)
