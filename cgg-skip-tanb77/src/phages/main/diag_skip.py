"""Post-freeze diagnostics (not a frozen test): skip/O by set, GC bin, CheckV host markers / trimmed / quality.
usage: python -I diag_skip.py"""
import csv, subprocess
from collections import defaultdict
sets = {r["id"]: (r["set"], int(r["gc_bin"])) for r in csv.DictReader(open("data/groups/run_ids.tsv"), delimiter="\t")}
st = {}
for f in ("results/main/target.tsv", "results/main/ctrl11_plant0.tsv"):
    for r in csv.DictReader(open(f), delimiter="\t"):
        st[r["id"]] = r
p = subprocess.run(["zstd", "-dcq", "data/uhgv/votus_metadata_extended.tsv.zst"],
                   capture_output=True, text=True).stdout.splitlines()
ix = {k: i for i, k in enumerate(p[0].split("\t"))}
meta = {}
for line in p[1:]:
    f = line.split("\t")
    if f[0] in st:
        meta[f[0]] = {k: f[ix[k]] for k in ("checkv_host_markers", "checkv_trimmed", "checkv_quality", "lifestyle",
                                             "genomad_plasmid_hallmarks", "cds_density")}
def agg(ids):
    O = sum(int(st[i]["O"]) for i in ids); S = sum(int(st[i]["skip_sites"]) for i in ids)
    n1 = sum(int(st[i]["skip_sites"]) > 0 for i in ids)
    return f"n={len(ids):5d} O={O:6d} skip={S:5d} skip/O={S / O if O else float('nan'):.4f} phages_with_skip={n1}"
print("by set, bins 24-30 only")
for s in ("T1", "T0", "Tf", "R"):
    ids = [i for i in st if sets[i][0] == s and 24 <= sets[i][1] <= 30]
    print(" ", s, agg(ids))
print("by set x bin")
for s in ("T1", "T0", "Tf", "R"):
    for b in range(24, 54, 2):
        ids = [i for i in st if sets[i] == (s, b)]
        if ids:
            print(" ", s, b, agg(ids))
print("by CheckV host markers / trimmed / quality / lifestyle, bins 24-30")
for k, f in (("host_markers=0", lambda m: m["checkv_host_markers"] in ("0", "")),
             ("host_markers>0", lambda m: m["checkv_host_markers"] not in ("0", "")),
             ("trimmed=Yes", lambda m: m["checkv_trimmed"] in ("Yes", "True", "1")),
             ("trimmed=No", lambda m: m["checkv_trimmed"] not in ("Yes", "True", "1")),
             ("Complete", lambda m: m["checkv_quality"] == "Complete"),
             ("temperate", lambda m: m["lifestyle"] == "temperate"),
             ("virulent", lambda m: m["lifestyle"] == "virulent")):
    for s in ("T1", "T0", "Tf", "R"):
        ids = [i for i in st if sets[i][0] == s and 24 <= sets[i][1] <= 30 and f(meta[i])]
        print(" ", k, s, agg(ids))
print("values seen:", {k: sorted({m[k] for m in meta.values()})[:6] for k in ("checkv_trimmed", "lifestyle")})
