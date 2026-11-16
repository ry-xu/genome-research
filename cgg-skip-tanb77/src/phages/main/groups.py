"""Build T (CRISPR host o__TANB77) and R (CRISPR host c__Clostridia, other orders) ID lists per RULES.md.
usage: python -I groups.py METADATA.tsv.zst OUTDIR"""
import csv, io, sys, subprocess
from collections import Counter

meta, outd = sys.argv[1], sys.argv[2]
p = subprocess.Popen(["zstd", "-dcq", meta], stdout=subprocess.PIPE, text=True)
rows = csv.DictReader(p.stdout, delimiter="\t")
OKQ = {"Complete", "High-quality", "Medium-quality"}
T, R, rec, stat = [], [], [], Counter()
for r in rows:
    h = r["crispr_host_lineage_gtdb_r207"]
    if "c__Clostridia;" not in h:
        continue
    order = h.split("o__")[1].split(";")[0] if "o__" in h else ""
    stat["clos", r["genetic_code"], r["checkv_quality"], order == "TANB77"] += 1
    if r["checkv_quality"] not in OKQ:
        continue
    if r["genetic_code"] != "11":
        rec.append((r["uhgv_genome"], r["genetic_code"], order))
        continue
    grp = "T" if order == "TANB77" else "R"
    (T if grp == "T" else R).append((r["uhgv_genome"], float(r["gc_percent"]), h))
for name, lst in (("T", T), ("R", R)):
    with open(f"{outd}/{name}.tsv", "w") as f:
        f.write("id\tgc\thost\n")
        for i, g, h in lst:
            f.write(f"{i}\t{g}\t{h}\n")
with open(f"{outd}/recoded_clos.tsv", "w") as f:
    f.write("id\tgenetic_code\torder\n")
    for x in rec:
        f.write("\t".join(x) + "\n")
print("T", len(T), "R", len(R), "recoded", len(rec))
print("T genera", len({h.split("g__")[1].split(";")[0] if "g__" in h else h for _, _, h in T}), Counter("g__" in h for _, _, h in T))
print("T gc range", min(g for _, g, _ in T), max(g for _, g, _ in T))
for k, v in sorted(stat.items(), key=lambda x: -x[1])[:15]:
    print(k, v)
