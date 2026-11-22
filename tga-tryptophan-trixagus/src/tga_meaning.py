"""What does each TGA-flagged genome put at TGA? The read-through tables (results/rt_full, else results/rt) record
the Swiss-Prot residue aligned to the read-through stop (sres). Code 4 predicts W, code 25 predicts G.
usage: python3 src/tga_meaning.py  -> results/tga_meaning.tsv (one row per TGA-flagged genome in results/labels.tsv)"""
import csv
import os
from collections import Counter
from pathlib import Path

R = Path(__file__).resolve().parent.parent / "results"
rows = [r for r in csv.DictReader(open(R / "labels.tsv"), delimiter="\t") if r["stop"] == "TGA"]
with open(R / "tga_meaning.tsv", "w") as out:
    out.write("acc\tknown\torganism\tn\tW\tG\ttop4\tcall\n")
    for r in rows:
        p = next((R / d / f"{r['acc']}.tsv" for d in ("rt_full", "rt") if os.path.exists(R / d / f"{r['acc']}.tsv")), None)
        c = Counter(x["sres"] for x in csv.DictReader(open(p), delimiter="\t") if x["stop"] == "TGA") if p else Counter()
        n = sum(c.values())
        w, g = (c["W"] / n, c["G"] / n) if n else (0, 0)
        call = "Trp" if w > 2 * g and c.most_common(1)[0][0] == "W" else "Gly" if g > 2 * w and c.most_common(1)[0][0] == "G" else "?"
        out.write(f"{r['acc']}\t{r['known']}\t{r['organism']}\t{n}\t{w:.2f}\t{g:.2f}\t"
                  f"{' '.join(f'{k}{v}' for k, v in c.most_common(4))}\t{call}\n")
