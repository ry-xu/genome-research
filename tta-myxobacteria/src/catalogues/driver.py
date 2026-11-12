"""Stream-screen genomes: download a batch from NCBI -> pyrodigal -> marker AAI screen -> keep only provisional hits.

usage: uv run --no-project --with pyrodigal python -I src/catalogues/driver.py ACC_LIST TAG [batch=60]
  data/catalogues/<TAG>/fna/       genomes (deleted after screening unless provisional hit)
  data/catalogues/<TAG>/faa/       proteins (deleted after screening unless provisional hit)
  results/catalogues/screen_<TAG>.tsv   one row per screened genome
Provisional hit (kept on disk): clade_aai >= 60 or top_hit_clade_frac >= 0.3 (generous; final call is made later).
"""
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import screen  # noqa: E402

UV = ["uv", "run", "--no-project", "python", "-I"]


def main():
    accs = [l.split()[0] for l in open(sys.argv[1]) if l.strip()]
    tag = sys.argv[2]
    bs = int(sys.argv[3]) if len(sys.argv) > 3 else 60
    base = ROOT / "data/catalogues" / tag
    fna, faa = base / "fna", base / "faa"
    fna.mkdir(parents=True, exist_ok=True)
    faa.mkdir(parents=True, exist_ok=True)
    out = ROOT / f"results/catalogues/screen_{tag}.tsv"
    done = set()
    if out.exists():
        done = {l.split("\t")[0] for l in open(out)}
    todo = [a for a in accs if a not in done]
    print("todo", len(todo), flush=True)
    for k in range(0, len(todo), bs):
        b = todo[k:k + bs]
        lst = base / "batch.txt"
        lst.write_text("\n".join(b) + "\n")
        t0 = time.time()
        subprocess.run(UV + [str(ROOT / "src/catalogues/fetch_ncbi.py"), str(lst), str(fna), "30"], check=False)
        t1 = time.time()
        files = [str(fna / f"{a}.fna.gz") for a in b if (fna / f"{a}.fna.gz").exists()]
        screen.predict(str(faa), "FAST", files)
        t2 = time.time()
        rows = screen.run(str(ROOT / "data/catalogues/markerdb"), str(faa), str(out)) if any(faa.glob("*.faa")) else []
        got = {r["genome"] for r in rows}
        with open(out, "a") as o:  # genomes that failed download or had too little sequence
            for a in b:
                if a not in got:
                    o.write(f"{a}\tNA\n")
        keep = 0
        for r in rows:
            hit = float(r["clade_aai"]) >= 60 or float(r["top_hit_clade_frac"]) >= 0.3
            keep += hit
            if not hit:
                (fna / f"{r['genome']}.fna.gz").unlink(missing_ok=True)
            # faa always moved out of the batch dir so the next batch only screens new genomes
            p = faa / f"{r['genome']}.faa"
            if hit:
                (base / "hits_faa").mkdir(exist_ok=True)
                p.rename(base / "hits_faa" / p.name)
            else:
                p.unlink(missing_ok=True)
        for p in faa.glob("*"):
            p.unlink()
        print(time.strftime("%H:%M:%S"), f"{k + len(b)}/{len(todo)} dl {t1 - t0:.0f}s pred {t2 - t1:.0f}s "
              f"screen {time.time() - t2:.0f}s kept {keep}", flush=True)


if __name__ == "__main__":
    main()
