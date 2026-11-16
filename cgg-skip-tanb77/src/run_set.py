"""Scan a genome set in batches, one DIAMOND call at a time; resumable (skips genomes whose tables exist).
usage: uv run python -I src/run_set.py <set.tsv> <outdir> [batch=16] [lane=0/1]
lane i/n takes every n-th genome of the ordered list, so lanes never share a genome (and their logs are separate).
Genomes are grouped by GTDB's NCBI translation table (4, 25, else 11) so each batch runs with its own code. Coverage
grows in tiers: the best-scoring rep of each class, then of each order, then every other family rep. A genome not on disk
is checked again before each batch and logged as missing at the end."""
import csv
import os
import subprocess
import sys
import time
from pathlib import Path

from codon_atlas.cache import genome_path

ROOT = Path(__file__).resolve().parents[1]


def main():
    rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
    outdir = Path(sys.argv[2])
    outdir.mkdir(parents=True, exist_ok=True)
    batch = int(sys.argv[3]) if len(sys.argv) > 3 else 16
    tier = {}
    for lvl, key in ((0, ("domain", "phylum", "class")), (1, ("domain", "phylum", "class", "order"))):
        best = {}
        for r in rows:
            k = tuple(r[c] for c in key)
            if k not in best or float(r["score"]) > float(best[k]["score"]):
                best[k] = r
        for r in best.values():
            tier.setdefault(r["gca"], lvl)
    rows.sort(key=lambda r: tier.get(r["gca"], 2))
    lane, nlanes = map(int, (sys.argv[4] if len(sys.argv) > 4 else "0/1").split("/"))
    rows = rows[lane::nlanes]
    log = open(outdir / f"run.{lane}.log", "a")
    todo = []
    for r in rows:
        if (outdir / f"{r['gca']}.fs.tsv.gz").exists():
            continue
        code = r["ttable"] if r["ttable"] in ("4", "25", "16") else "11"
        todo.append((code, r["gca"], r["local"] or str(genome_path(r["gca"]))))
    log.write(f"{time.strftime('%H:%M:%S')}\tstart\t{len(todo)} to do\n")
    # keep the tier priority, but batch within a code
    while todo:
        ready = [t for t in todo if os.path.exists(t[2])]
        if not ready:
            break
        code = ready[0][0]
        b = [t for t in ready if t[0] == code][:batch]
        todo = [t for t in todo if t not in b]
        t0 = time.time()
        p = subprocess.run([sys.executable, "-I", str(ROOT / "src/fs_scan.py"), str(outdir)]
                           + [f"{n}={pa}" for _, n, pa in b], env=os.environ | {"GENCODE": code},
                           capture_output=True, text=True)
        log.write(f"{time.strftime('%H:%M:%S')}\tbatch\t{len(b)}\tcode {code}\t{time.time() - t0:.0f}s\trc {p.returncode}"
                  f"\t{p.stderr.strip()[-200:]}\n")
        log.flush()
    for t in todo:
        log.write(f"{time.strftime('%H:%M:%S')}\tmissing\t{t[1]}\n")


if __name__ == "__main__":
    main()
