"""Assign run sets (RULES draft): T1 skip-genus hosts, T0 non-skip TANB77 genera, Tf CAG-508 family-only, R sample
(<= 500 per 2-pt GC bin, bins holding T), P2 code-4 phages (all hosts), code-15 Clostridia phages.
usage: python -I pick.py OUT.tsv"""
import csv, random, subprocess, sys
from collections import defaultdict
call = {}
for r in csv.DictReader(open("scratch/cgg_test_T.tsv"), delimiter="\t"):
    call[r["genus"]] = r["call"]
T = list(csv.DictReader(open("data/groups/T.tsv"), delimiter="\t"))
R = list(csv.DictReader(open("data/groups/R.tsv"), delimiter="\t"))
b = lambda g: int(float(g) // 2) * 2
out = open(sys.argv[1], "w"); out.write("id\tset\tgc_bin\thost\n")
bins = set()
for r in T:
    h = r["host"].split(";")
    fam = h[4][3:] if len(h) > 4 else ""
    gen = h[5][3:] if len(h) > 5 else ""
    c = call.get(gen)
    s = "T1" if c == "1" else "T0" if c == "0" else ("Tf" if fam == "CAG-508" and not gen else "Tx")
    out.write(f"{r['id']}\t{s}\t{b(r['gc'])}\t{r['host']}\n"); bins.add(b(r["gc"]))
byb = defaultdict(list)
for r in R:
    byb[b(r["gc"])].append(r)
rng = random.Random(20261115)
for k in sorted(bins):
    for r in rng.sample(byb[k], min(500, len(byb[k]))):
        out.write(f"{r['id']}\tR\t{k}\t{r['host']}\n")
p = subprocess.run(["zstd", "-dcq", "data/uhgv/votus_metadata_extended.tsv.zst"],
                   capture_output=True, text=True).stdout.splitlines()
hd = p[0].split("\t"); ix = {k: i for i, k in enumerate(hd)}
for line in p[1:]:
    f = line.split("\t")
    if f[ix["checkv_quality"]] not in ("Complete", "High-quality", "Medium-quality"):
        continue
    if f[ix["genetic_code"]] == "4":
        out.write(f"{f[0]}\tP2\t{b(f[ix['gc_percent']])}\t{f[ix['crispr_host_lineage_gtdb_r207']]}\n")
    elif f[ix["genetic_code"]] == "15" and "c__Clostridia;" in f[ix["crispr_host_lineage_gtdb_r207"]]:
        out.write(f"{f[0]}\tC15\t{b(f[ix['gc_percent']])}\t{f[ix['crispr_host_lineage_gtdb_r207']]}\n")
