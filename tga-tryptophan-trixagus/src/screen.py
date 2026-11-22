"""Read-through screen of QC-failed genomes: stop-codon-census batch_pooled.py (prep/search reused unchanged) with
this project's output paths and the 'qc' sampling rates. The path overrides are module-level so spawned Pool workers,
which import this file as __mp_main__ before unpickling tasks, see them too.

usage: KEEP_GENOMES=0 uv run --no-project --with pyrodigal,requests python src/screen.py <acc_list> [workers=3] [threads=3] [pool=50]
       acc_list: one accession per line, optional 2nd column = rates key (default qc; "full" = every gene)
"""
import sys
import time
from multiprocessing import Pool
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent.parent / "stop-codon-census" / "src"))
import readthrough as rt  # noqa: E402
import batch  # noqa: E402
import batch_pooled as bp  # noqa: E402

import os  # noqa: E402
rt.OUT = HERE.parent / "results" / os.environ.get("OUTDIR", "rt")  # stage 2 (all genes): OUTDIR=rt_full, key full
rt.GEN = HERE.parent / "data" / "genomes"
batch.RATES["qc"] = {"TAA": 0.2, "TAG": 0.3, "TGA": 0.3}  # same dict object bp.RATES points to

if __name__ == "__main__":
    lines = [l.split() for l in open(sys.argv[1]) if l.strip()]
    workers = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    threads = int(sys.argv[3]) if len(sys.argv) > 3 else 3
    pool = int(sys.argv[4]) if len(sys.argv) > 4 else 50
    tasks = [(x[0].replace("GB_", "").replace("RS_", ""), x[1] if len(x) > 1 else "qc") for x in lines]
    tasks = [t for t in tasks if not (rt.OUT / f"{t[0]}.stats.json").exists()]

    def log(msg):
        print(f"{time.strftime('%H:%M:%S')}\t{msg}", flush=True)

    log(f"# {len(tasks)} genomes to do ({len(lines) - len(tasks)} already done); OUT={rt.OUT}")
    buf, n = [], 0
    with Pool(workers) as p:
        for acc, res, dt in p.imap_unordered(bp.prep, tasks):
            n += 1
            if isinstance(res, str):
                log(f"{acc}\t{res}")
                continue
            buf.append((acc, res))
            if len(buf) >= pool:
                bp.search(buf, threads, log)
                log(f"# {n}/{len(tasks)}")
                buf = []
        if buf:
            bp.search(buf, threads, log)
    log(f"done: {n}/{len(tasks)}")
