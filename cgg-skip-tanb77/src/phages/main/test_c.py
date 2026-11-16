"""Frozen test C (RULES.md; implemented after the freeze, spec unchanged): share of phages whose CRISPR host range
crosses orders (hrb_rank order/class/phylum/domain, crispr_spacer rows) in T1+Tf vs R; Fisher exact; per-GC-bin table;
skip signature of each crossing T phage. usage: python -I test_c.py"""
import csv, subprocess
from collections import defaultdict
from scipy.stats import fisher_exact
sets = {r["id"]: (r["set"], int(r["gc_bin"])) for r in csv.DictReader(open("data/groups/run_ids.tsv"), delimiter="\t")}
p = subprocess.run(["zstd", "-dcq", "data/uhgv/votus_metadata_extended.tsv.zst"],
                   capture_output=True, text=True).stdout.splitlines()
ix = {k: i for i, k in enumerate(p[0].split("\t"))}
votu = {f[0]: f[ix["uhgv_votu"]] for f in (l.split("\t") for l in p[1:]) if f[0] in sets}
hrb = {}
for r in csv.DictReader(subprocess.run(["zstd", "-dcq", "data/uhgv/host_range_breadth.tsv.zst"], capture_output=True,
                                       text=True).stdout.splitlines(), delimiter="\t"):
    if r["method"] == "crispr_spacer":
        hrb[r["uhgv_votu"]] = r["hrb_rank"]
st = {}
for f in ("results/main/target.tsv", "results/main/ctrl11_plant0.tsv"):
    for r in csv.DictReader(open(f), delimiter="\t"):
        st[r["id"]] = r
CROSS = {"order", "class", "phylum", "domain"}
cnt = defaultdict(lambda: [0, 0, 0])  # crossing, total, missing hrb
byb = defaultdict(lambda: defaultdict(lambda: [0, 0]))
crossers = []
for i, (s, b) in sets.items():
    g = "T" if s in ("T1", "Tf") else "R" if s == "R" else None
    if g is None:
        continue
    rk = hrb.get(votu.get(i))
    if rk is None:
        cnt[g][2] += 1; continue
    c = rk in CROSS
    cnt[g][0] += c; cnt[g][1] += 1; byb[b][g][0] += c; byb[b][g][1] += 1
    if c and g == "T":
        crossers.append((i, s, b, rk, st[i]["O"], st[i]["skip_sites"]))
for g in ("T", "R"):
    print(g, "crossing", cnt[g][0], "of", cnt[g][1], f"= {cnt[g][0] / max(1, cnt[g][1]):.4f}", "(no hrb row:", cnt[g][2], ")")
t, r = cnt["T"], cnt["R"]
print("Fisher (T vs R):", fisher_exact([[t[0], t[1] - t[0]], [r[0], r[1] - r[0]]]))
for b in sorted(byb):
    print(" bin", b, {g: f"{v[0]}/{v[1]}" for g, v in byb[b].items()})
print("crossing T phages (id, set, bin, hrb_rank, O, skip):")
for c in crossers:
    print(" ", *c)
