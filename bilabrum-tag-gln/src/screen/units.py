"""Build data/units.tsv: one row per stored new-stops-euk nuclear scan (read-only inputs).

Transcriptomes: data/tsa and data/mmetsp, named from frameshift-codes-euk's final_calls.tsv (with its duplicate_of).
Genomes: data/genomes/scan2 (v2) where finished, else data/genomes/scan (v1); named from Clover's genus tables.
NCBI taxid and its nuclear genetic code id (taxdump nodes.dmp, Clover's copy of 2026-11-10) are the 'known' label.
usage: units.py > data/units.tsv
"""
import csv
import glob
import os

NS = "lab/work/new-stops-euk/data"
FC = "lab/work/frameshift-codes-euk/results/final_calls.tsv"

TAX = "lab/work/new-stops-euk/data/taxdump_nodes/nodes.dmp"
gc = {}
for line in open(TAX):
    x = line.split("\t|\t")
    gc[x[0]] = x[6]
names, taxid = {}, {}
for r in csv.DictReader(open(FC), delimiter="\t"):
    names[r["id"]] = (r["name"], r["group"], r["duplicate_of"])
    taxid[r["id"]] = r["taxid"]
for t in ("genus.tsv", "genus2.tsv", "overcap2.tsv"):
    p = os.path.join(NS, "genomes", t)
    if os.path.exists(p):
        for line in open(p):
            x = line.rstrip("\n").split("\t")
            if len(x) >= 5 and x[0].startswith("GC"):
                names.setdefault(x[0], (x[3], x[4], ""))
                taxid.setdefault(x[0], x[2])

print("unit\tstratum\tprefix\tname\tgroup\tduplicate_of\ttaxid\tgc")
for st in ("tsa", "mmetsp"):
    for f in sorted(glob.glob(os.path.join(NS, st, "*.codons.npz"))):
        u = os.path.basename(f)[:-len(".codons.npz")]
        if os.path.exists(os.path.join(NS, st, u + ".ends.tsv")):
            n, g, d = names.get(u, ("", "", ""))
            t = taxid.get(u, "")
            print(f"{u}\t{st}\t{f[:-len('.codons.npz')]}\t{n}\t{g}\t{d}\t{t}\t{gc.get(t, '')}")
seen = set()
for st in ("scan2", "scan"):
    for f in sorted(glob.glob(os.path.join(NS, "genomes", st, "*.codons.npz"))):
        u = os.path.basename(f)[:-len(".codons.npz")]
        if u in seen or not os.path.exists(os.path.join(NS, "genomes", st, u + ".ends.tsv")):
            continue
        seen.add(u)
        n, g, d = names.get(u, ("", "", ""))
        t = taxid.get(u, "")
        print(f"{u}\tgenome_{st}\t{f[:-len('.codons.npz')]}\t{n}\t{g}\t\t{t}\t{gc.get(t, '')}")
