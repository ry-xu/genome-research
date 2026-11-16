"""One genome per GTDB R232 family: the best-assembled representative, preferring one already on disk.
-> data/family_set.tsv (reps.tsv columns + score). Score = CheckM2 completeness - 5 x contamination - 10 x log10(contigs)."""
import csv
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
rows = list(csv.DictReader(open(ROOT / "data/reps.tsv"), delimiter="\t"))


def score(r):
    try:
        return float(r["compl"]) - 5 * float(r["contam"]) - 10 * math.log10(max(1, int(r["contigs"])))
    except ValueError:
        return -1e9


fam = {}
for r in rows:
    fam.setdefault((r["domain"], r["phylum"], r["family"]), []).append(r)
out = []
for k, v in fam.items():
    loc = [r for r in v if r["local"]]
    best = max(v, key=score)
    # take a local genome unless it is much worse than the best one (score margin 20)
    if loc and score(max(loc, key=score)) >= score(best) - 20:
        best = max(loc, key=score)
    out.append(best | {"score": f"{score(best):.1f}"})
cols = list(rows[0].keys()) + ["score"]
with open(ROOT / "data/family_set.tsv", "w") as fh:
    w = csv.DictWriter(fh, cols, delimiter="\t")
    w.writeheader()
    w.writerows(sorted(out, key=lambda r: (r["domain"], r["phylum"], r["family"])))
print(len(out), "families;", sum(1 for r in out if r["local"]), "local;",
      sum(1 for r in out if not r["local"]), "to fetch")
