"""Write ARAGORN tRNA calls in the genome-cache layout that src/analysis/cgg_test.py and src/recheck/ read:
<CODON_ATLAS>/products/<acc[:8]>/<acc>/trna.parquet, columns as in codon_atlas/trna.py (the settings of the GTDB-wide
tRNA screen). Genomes are fetched from NCBI into <CODON_ATLAS>/genomes if missing.
usage: uv run python -I src/trna_products.py <set.tsv with a gca column> | <acc> ..."""
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import duckdb  # noqa: E402

from codon_atlas import cache, trna  # noqa: E402
from paths import PRODUCTS  # noqa: E402


def write(acc):
    out = PRODUCTS / acc[:8] / acc / "trna.parquet"
    if out.exists():
        return out
    rows = trna.call(cache.read_fasta(cache.fetch(acc)))
    out.parent.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect()
    con.execute("create table t (" + ", ".join(f"{c} varchar" for c in trna.COLS) + ")")
    con.executemany(f"insert into t values ({', '.join('?' * len(trna.COLS))})", [[str(x) for x in r] for r in rows])
    con.execute(f"copy t to '{out}' (format parquet)")
    return out


if __name__ == "__main__":
    args = sys.argv[1:]
    if len(args) == 1 and args[0].endswith(".tsv"):
        args = [r["gca"] for r in csv.DictReader(open(args[0]), delimiter="\t")]
    for a in args:
        print(a, write(a), flush=True)
