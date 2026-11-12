"""Download NCBI assemblies (genomic FASTA only) via the Datasets v2 REST API, in batches, as <acc>.fna.gz.

usage: uv run --no-project python -I src/catalogues/fetch_ncbi.py ACC_LIST_FILE OUTDIR [batch=40]
Skips accessions whose <acc>.fna.gz already exists (or is listed in OUTDIR/../screened_done.txt).
"""
import gzip
import io
import shutil
import sys
import time
import urllib.request
import zipfile
from pathlib import Path

URL = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{}/download?include_annotation_type=GENOME_FASTA"


def fetch(batch, out):
    for i in range(5):
        try:
            with urllib.request.urlopen(URL.format(",".join(batch)), timeout=900) as r:
                data = r.read()
            z = zipfile.ZipFile(io.BytesIO(data))
            got = set()
            for n in z.namelist():
                if n.endswith(".fna") and "/data/" in n:
                    acc = n.split("/")[-2]
                    if not acc.startswith("GC"):
                        continue
                    tmp = out / f".{acc}.fna.gz.part"
                    with z.open(n) as src, gzip.open(tmp, "wb", compresslevel=4) as dst:
                        shutil.copyfileobj(src, dst)
                    tmp.rename(out / f"{acc}.fna.gz")
                    got.add(acc)
            return got
        except Exception as e:
            print("retry", i, e, flush=True)
            time.sleep(5 + 10 * i)
    return set()


def main():
    accs = [l.split()[0] for l in open(sys.argv[1]) if l.strip() and not l.startswith("#")]
    out = Path(sys.argv[2])
    out.mkdir(parents=True, exist_ok=True)
    bs = int(sys.argv[3]) if len(sys.argv) > 3 else 40
    done_f = out.parent / "screened_done.txt"
    done = set(open(done_f).read().split()) if done_f.exists() else set()
    todo = [a for a in accs if a not in done and not (out / f"{a}.fna.gz").exists()]
    print("to fetch", len(todo), flush=True)
    for k in range(0, len(todo), bs):
        b = todo[k:k + bs]
        got = fetch(b, out)
        miss = [a for a in b if a not in got]
        print(time.strftime("%H:%M:%S"), k + len(b), "/", len(todo), "got", len(got), "missing", miss[:5], flush=True)
        time.sleep(0.5)


if __name__ == "__main__":
    main()
