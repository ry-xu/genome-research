"""Fetch NCBI assemblies into data/genomes/ (reuses tag-bacteria's fetcher: Datasets API, FTP fallback).
usage: python src/fetch.py <acc> [<acc> ...]"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT.parent / "tag-bacteria" / "src"))
import readthrough as rt

rt.GEN = ROOT / "data" / "genomes"


def fetch(acc):
    census = ROOT.parent / "tag-bacteria" / "data" / "genomes" / f"{acc}.fna.gz"
    if census.exists():
        return census
    return rt.fetch(acc)


if __name__ == "__main__":
    for a in sys.argv[1:]:
        try:
            print(a, fetch(a), flush=True)
        except Exception as e:
            print(a, "FAILED", e, flush=True)
