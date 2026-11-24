"""Ledger run order: control genomes first (the RNA survey's 81 units, L02's GCA_025839675.1), then Asgard, then the
rest. Writes data/run1.tsv (cached genomes) and data/run2.tsv (fetched Asgard reps). Paths made absolute."""
import csv
ATLAS = "/Users/ryan/lab/shared/codon-atlas/"
meta = list(csv.DictReader(open("data/genomes_meta.tsv"), delimiter="\t"))
path = {}
for l in open("data/genomes.tsv"):
    a, p, _ = l.rstrip("\n").split("\t")
    path[a] = p if p.startswith("/") else ATLAS + p
survey = {r["assembly"] for r in csv.DictReader(open("/Users/ryan/lab/work/archaeal-introns/survey/data/units.tsv"), delimiter="\t")}
ctrl = survey | {"GCA_025839675.1"}
def key(d):
    return (0 if d["acc"] in ctrl else 1 if d["phylum"] == "Asgardarchaeota" else 2, d["phylum"], d["acc"])
r1 = sorted((d for d in meta if d["tier"] == "cached" and d["acc"] in path), key=key)
r2 = sorted((d for d in meta if d["tier"] == "asgard_fetch" and d["acc"] in path), key=key)
for fn, rows in (("data/run1.tsv", r1), ("data/run2.tsv", r2)):
    with open(fn, "w") as f:
        for d in rows:
            f.write(f"{d['acc']}\t{path[d['acc']]}\t0\n")
    print(fn, len(rows), "controls:", sum(d["acc"] in ctrl for d in rows), "asgard:", sum(d["phylum"] == "Asgardarchaeota" for d in rows))
print("survey units in list:", len(survey & set(path)), "of", len(survey))
