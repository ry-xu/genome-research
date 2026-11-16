# Post hoc (descriptive): skip/O by host-call evidence strength (CRISPR connections) and by the independent
# k-mer host call (TANB77 or not), in T1+Tf and R.
import csv, subprocess, io, collections
META = "data/uhgv/votus_metadata_extended.tsv.zst"
sets = {r["id"]: r["set"] for r in csv.DictReader(open("data/groups/run_ids.tsv"), delimiter="\t")}
stats = {}
for f in ("results/main/target.tsv", "results/main/ctrl11_plant0.tsv"):
    for r in csv.DictReader(open(f), delimiter="\t"):
        if r["id"] in sets: stats[r["id"]] = r
txt = subprocess.run(["zstd", "-dcq", META], capture_output=True, text=True).stdout
m = {r["uhgv_genome"]: r for r in csv.DictReader(io.StringIO(txt), delimiter="\t") if r["uhgv_genome"] in stats}
def km(i):
    l = m[i]["kmers_host_lineage_gtdb_r207"]
    return "none" if l in ("", "NULL") else ("TANB77" if "o__TANB77" in l else "other")
def cbin(i):
    c = int(m[i]["crispr_total_connections"]); return "1" if c <= 1 else "2-4" if c <= 4 else ">=5"
acc = collections.defaultdict(lambda: [0, 0, 0])
for i, s in sets.items():
    if i not in stats: continue
    g = "T" if s in ("T1", "Tf") else s
    if g not in ("T", "R"): continue
    gc = "lowGC" if float(stats[i]["gc"]) < 0.32 else "GC>=32"
    for key in ((g, "conn", cbin(i)), (g, "kmer", km(i)), (g, "kmer x GC", km(i) + "/" + gc), (g, "conn x GC", cbin(i) + "/" + gc)):
        a = acc[key]; a[0] += 1; a[1] += int(stats[i]["skip_sites"]); a[2] += int(stats[i]["O"])
for k in sorted(acc):
    n, S, O = acc[k]; print(*k, f"n={n} skip/O={S}/{O}={S/O if O else float('nan'):.3f}", sep="\t")
