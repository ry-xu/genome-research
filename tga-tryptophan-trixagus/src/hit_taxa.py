"""Where do a genome's read-through homologs come from? Tally Swiss-Prot species mnemonics (and mitochondrial /
plastid gene names) of the hits for one stop. A flag carried by host mitochondria or a eukaryote shows up here.
usage: python3 src/hit_taxa.py <stop> <acc> [...]
"""
import csv
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
ORG = ("MOUSE HUMAN BOVIN PIG RAT CHICK XENLA DANRE YEAST SCHPO ARATH MAIZE ORYSJ TOBAC CAEEL DROME".split())


def main():
    stop = sys.argv[1]
    for acc in sys.argv[2:]:
        f = ROOT / "results" / "rt_full" / f"{acc}.tsv"
        if not f.exists():
            f = ROOT / "results" / "rt" / f"{acc}.tsv"
        rows = [r for r in csv.DictReader(open(f), delimiter="\t") if r["stop"] == stop]
        sp = Counter(r["sseqid"].split("_")[-1] for r in rows)
        mito = sum(1 for r in rows if r["sseqid"].split("|")[-1].split("_")[0][-1:] == "M" or "mitochondr" in
                   r["stitle"].lower() or r["sseqid"].split("_")[-1] in ORG)
        print(acc, stop, len(rows), "mito/euk-like", mito, sp.most_common(8))


if __name__ == "__main__":
    main()
