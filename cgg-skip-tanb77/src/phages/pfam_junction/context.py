"""Sequence context of SKIP/CTRL sites listed in out/term_pos.tsv (sid, hmm_pos). Rebuilds set order as build.py."""
import csv, random, sys
CW = "phages_work"
rows = lambda p: list(csv.DictReader(open(p), delimiter="\t"))
cgg = rows(f"{CW}/kill1/out/meta.cgg_sites.tsv")
skip = [r for r in cgg if r["set"] in ("T1", "Tf") and r["skip"] == "1"]
meta = {r["id"]: r for r in rows(f"{CW}/kill1/out/meta_sub.tsv")}
want = {}
for l in open(sys.argv[1]):
    sid, pos = l.split()
    if sid.startswith("SKIP"):
        want[int(sid.split("_")[1])] = (sid, pos)
need = {skip[i]["id"] for i in want}
seqs, name, buf = {}, None, []
for l in open(f"{CW}/data/phages.fna"):
    if l[0] == ">":
        if name in need: seqs[name] = "".join(buf).upper()
        name, buf = l[1:].split()[0], []
    else: buf.append(l.strip())
if name in need: seqs[name] = "".join(buf).upper()
COMP = str.maketrans("ACGTN", "TGCAN")
for i in sorted(want, key=lambda i: int(want[i][1])):
    r = skip[i]; g = seqs[r["id"]]; L = len(g)
    b, e, st, k = int(r["begin"]), int(r["end"]), int(r["strand"]), int(r["k"])
    s = g if st == 1 else g.translate(COMP)[::-1]
    p = (b - 1 if st == 1 else L - e) + 3 * k
    m = meta[r["id"]]
    print(want[i][0], want[i][1], r["id"], m["vgenus"], m["vfam"], m["quality"], s[p-15:p].lower() + " " + s[p:p+3] + " " + s[p+3:p+18].lower(), f"k={k} n={r['n']}", sep="\t")
