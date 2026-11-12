"""Find the TTA clade's MRCA in the GTDB R232 bac120 tree and its nearest outgroup genera (for negative controls).
usage: uv run --no-project --with treeswift python -I src/catalogues/tree_neighbours.py data/catalogues/gtdb/bac120_r232.tree.gz
"""
import gzip, sys, csv
import treeswift
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
CLADE = {"g__JADLFU01", "g__JANWXX01", "g__DASTPP01", "g__JAXLYQ01", "g__CAYNEK01"}
tax = {}
with gzip.open(ROOT / "data/gtdb_meta/bac120_metadata.tsv.gz", "rt") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        if "f__Polyangiaceae;" in r["gtdb_taxonomy"] and r["gtdb_representative"] == "t":
            tax[r["accession"]] = r["gtdb_taxonomy"].split(";")[5]
t = treeswift.read_tree_newick(gzip.open(sys.argv[1], "rt").read())
lab = {n.label: n for n in t.traverse_leaves()}
clade = [lab[a] for a, g in tax.items() if g in CLADE and a in lab]
print("clade reps in tree:", len(clade))
m = t.mrca({n.label for n in clade})
def tips(n): return [x.label for x in n.traverse_leaves()]
inside = tips(m)
print("MRCA tips:", len(inside), "label:", m.label)
for a in inside: print("  in", a, tax.get(a, "?"))
n, k = m, 0
while n.parent is not None and k < 4:
    p = n.parent; k += 1
    sib = [c for c in p.children if c is not n]
    print(f"--- level {k} parent label={p.label}")
    for s in sib:
        ts = tips(s)
        gs = {}
        for a in ts: gs[tax.get(a, "non-Polyangiaceae")] = gs.get(tax.get(a, "non-Polyangiaceae"), 0) + 1
        print("  sister subtree", len(ts), "tips; genera:", dict(sorted(gs.items(), key=lambda x: -x[1])[:15]))
        if len(ts) <= 30: print("   ", ts)
    n = p
