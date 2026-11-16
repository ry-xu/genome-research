# Post hoc (descriptive, after the 18:27 post): skip/O by phage family; T1/Tf vs R without vFAM-00025;
# high-GC T1 phages (GC >= 32%) by vFAM and host.
import csv, subprocess, io, collections
META = "data/uhgv/votus_metadata_extended.tsv.zst"
sets = {r["id"]: r for r in csv.DictReader(open("data/groups/run_ids.tsv"), delimiter="\t")}
stats = {}
for f in ("results/main/target.tsv", "results/main/ctrl11_plant0.tsv"):
    for r in csv.DictReader(open(f), delimiter="\t"):
        if r["id"] in sets: stats[r["id"]] = r
txt = subprocess.run(["zstd", "-dcq", META], capture_output=True, text=True).stdout
rd = csv.DictReader(io.StringIO(txt), delimiter="\t")
fam = {}; host = {}
for r in rd:
    i = r[rd.fieldnames[0]]
    if i in stats:
        t = r["uhgv_taxonomy"].split(";")
        fam[i] = next((x for x in t if x.startswith("vFAM-")), "none")
        host[i] = r["crispr_host_lineage_gtdb_r207"]
print("ids with meta:", len(fam), "of", len(stats))
def pool(ids, lo=24, hi=30):
    O = sum(int(stats[i]["O"]) for i in ids if lo <= float(stats[i]["gc"]) * 100 < hi + 2)
    S = sum(int(stats[i]["skip_sites"]) for i in ids if lo <= float(stats[i]["gc"]) * 100 < hi + 2)
    n = sum(1 for i in ids if lo <= float(stats[i]["gc"]) * 100 < hi + 2)
    return n, S, O, (S / O if O else float("nan"))
by = collections.defaultdict(list)
for i, r in sets.items():
    if i in stats: by[r["set"]].append(i)
print("\n== skip/O, GC 24-31.99 (bins 24-30), all vs without vFAM-00025")
for s in ("T1", "T0", "Tf", "R"):
    a = pool(by[s]); b = pool([i for i in by[s] if fam.get(i) != "vFAM-00025"])
    print(f"{s}: all n={a[0]} {a[1]}/{a[2]}={a[3]:.3f} | no vFAM-00025 n={b[0]} {b[1]}/{b[2]}={b[3]:.3f}")
print("\n== T1+Tf skip/O by vFAM (bins 24-30; vFAMs with O >= 20); R skip/O for the same vFAM")
c = collections.Counter(fam[i] for i in by["T1"] + by["Tf"])
for f_, n in c.most_common():
    ids = [i for i in by["T1"] + by["Tf"] if fam[i] == f_]; a = pool(ids)
    r = pool([i for i in by["R"] if fam[i] == f_])
    if a[2] >= 20: print(f"{f_}: T n={a[0]} {a[1]}/{a[2]}={a[3]:.3f} | R n={r[0]} {r[1]}/{r[2]}")
print("\n== T1 phages GC >= 32%:", sum(1 for i in by["T1"] if float(stats[i]["gc"]) >= 0.32))
hg = [i for i in by["T1"] if float(stats[i]["gc"]) >= 0.32]
O = sum(int(stats[i]["O"]) for i in hg); S = sum(int(stats[i]["skip_sites"]) for i in hg)
print(f"skip/O {S}/{O}; vFAMs:", collections.Counter(fam[i] for i in hg).most_common(6))
print("hosts:", collections.Counter(host[i].split(";")[-2] if ";" in host[i] else host[i] for i in hg).most_common(8))
lo = [i for i in by["T1"] if float(stats[i]["gc"]) < 0.32]
print("T1 <32% vFAMs:", collections.Counter(fam[i] for i in lo).most_common(6))
print("R phages in the high-GC T1 vFAMs:", {f_: sum(1 for i in by["R"] if fam[i] == f_) for f_, _ in collections.Counter(fam[i] for i in hg).most_common(4)})
