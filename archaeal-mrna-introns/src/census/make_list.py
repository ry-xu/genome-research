"""Genome list for the census: every archaeal genome in the codon-atlas cache, plus GTDB Asgard species reps not yet
cached (fetched into the atlas cache). Writes data/genomes.tsv (acc, path, tier) and data/genomes_meta.tsv."""
import csv, sys
from codon_atlas import imports, cache

con = imports.connect()
idx = {r["acc"]: r["path"] for r in csv.DictReader(open("/Users/ryan/lab/shared/codon-atlas/index/genomes.tsv"), delimiter="\t")}
arc = con.sql("select acc, phylum, class, \"order\", family, genus, species, gtdb_representative, checkm2_completeness, "
              "checkm2_contamination, contig_count, genome_size, gc_percentage from gtdb where gtdb_taxonomy like 'd__Archaea%'").fetchall()
cols = ["acc", "phylum", "class", "order", "family", "genus", "species", "rep", "compl", "contam", "contigs", "size", "gc"]
rows, fetch = [], []
for r in arc:
    d = dict(zip(cols, r))
    if d["acc"] in idx:
        d["tier"] = "cached"; rows.append(d)
    elif d["phylum"] in ("Asgardarchaeota", "p__Asgardarchaeota") and d["rep"] == "t":
        d["tier"] = "asgard_fetch"; rows.append(d); fetch.append(d["acc"])
print(len(rows), "genomes;", len(fetch), "to fetch", file=sys.stderr)
if "--fetch" in sys.argv:
    bad = []
    for i, a in enumerate(fetch):
        try:
            p = cache.fetch(a)
            idx[a] = str(p)
        except Exception as e:
            bad.append((a, repr(e)[:200]))
        if i % 25 == 0:
            print(i, a, file=sys.stderr, flush=True)
    for a, e in bad:
        print("FAILED", a, e, file=sys.stderr)
with open("data/genomes.tsv", "w") as f, open("data/genomes_meta.tsv", "w") as g:
    g.write("\t".join(cols + ["tier"]) + "\n")
    for d in rows:
        if d["acc"] in idx:
            p = idx[d["acc"]]; p = p if p.startswith("/") else "/Users/ryan/lab/shared/codon-atlas/" + p
            f.write(f"{d['acc']}\t{p}\t0\n")
            g.write("\t".join(str(d[c]) for c in cols + ["tier"]) + "\n")
