"""Fetch the genomes of a set into the shared atlas cache, 3 at a time (NCBI). Appends to <set>.fetch.log.
usage: uv run python -I src/fetch_set.py data/family_set.tsv"""
import csv
import sys
import time
from concurrent.futures import ThreadPoolExecutor

from codon_atlas.cache import fetch, genome_path

rows = [r for r in csv.DictReader(open(sys.argv[1]), delimiter="\t") if not r["local"]]
todo = [r["gca"] for r in rows if not genome_path(r["gca"]).exists()]
log = open(sys.argv[1] + ".fetch.log", "a")


def one(acc):
    try:
        fetch(acc)
        return acc, "ok"
    except Exception as e:
        return acc, f"fail {type(e).__name__}: {str(e)[:100]}"


with ThreadPoolExecutor(3) as ex:
    for i, (acc, st) in enumerate(ex.map(one, todo)):
        log.write(f"{time.strftime('%H:%M:%S')}\t{acc}\t{st}\n")
        log.flush()
print("done", len(todo))
