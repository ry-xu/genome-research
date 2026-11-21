"""Summarise loopsA1 tables: controls, bacterial loop counts (organelles excluded), non-GAGA genera.
usage: srl_silva_summary.py loopsA1_Ref.tsv"""
import csv, sys, collections

rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
def org(t): return "mito" if "Mitochondria" in t else "chloro" if "Chloroplast" in t else ""
ec = collections.Counter(r["loop"] for r in rows if r["taxonomy"].endswith(";Escherichia coli") or ";Escherichia coli " in r["taxonomy"])
print("PC1 E. coli (species field 'Escherichia coli*'):", dict(ec))
zy = collections.Counter((r["loop"], r["taxonomy"].split(";")[1] if ";" in r["taxonomy"] else "") for r in rows if "Zymomonas" in r["taxonomy"])
print("PC2 Zymomonas (loop, phylum):", dict(zy))
bac = [r for r in rows if r["taxonomy"].startswith("Bacteria") and not org(r["taxonomy"])]
read = [r for r in bac if r["loop"] and r["loop"] != "multi"]
c = collections.Counter(r["loop"] for r in read)
print(f"Bacteria, no organelle: {len(bac)} entries, {len(read)} read, multi {sum(r['loop']=='multi' for r in bac)}")
print("  loops:", c.most_common())
g = collections.Counter((r["loop"], r["taxonomy"].split(";")[5] if len(r["taxonomy"].split(";")) > 5 else r["taxonomy"]) for r in read if r["loop"] != "GAGA")
for (lp, gen), n in sorted(g.items(), key=lambda x: (-x[1], x[0])):
    print(f"  {lp}\t{gen}\t{n}")
for o in ("mito", "chloro"):
    print(o, collections.Counter(r["loop"] for r in rows if org(r["taxonomy"]) == o and r["loop"]).most_common(6))
print("Archaea read:", sum(1 for r in rows if r["taxonomy"].startswith("Archaea") and r["loop"] not in ("", "multi")))
