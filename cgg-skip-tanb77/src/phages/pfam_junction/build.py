"""Build fused (skip-read) proteins for SKIP, CTRL and RCGG sets (RULES.md, frozen 18:43).
usage: python -I src/build.py OUTDIR"""
import csv, random, sys
from collections import Counter, defaultdict

CW = "phages_work"
out = sys.argv[1]
grp = {r["id"]: r for r in csv.DictReader(open(f"{CW}/kill1/out/my_groups.tsv"), delimiter="\t")}


def rows(p):
    return list(csv.DictReader(open(p), delimiter="\t"))


cgg = rows(f"{CW}/kill1/out/meta.cgg_sites.tsv")
ctl = rows(f"{CW}/kill1/out/meta.ctrl_sites.tsv")
skip = [r for r in cgg if r["set"] in ("T1", "Tf") and r["skip"] == "1"]
rcgg = [r for r in cgg if r["set"] == "R" and r["skip"] == "1"]
gkey = lambda r: (r["id"], r["begin"], r["end"], r["strand"])
skip_genes = {gkey(r) for r in skip}
tphage = {r["id"] for r in cgg + ctl if r["set"] in ("T1", "Tf")}
bygene = defaultdict(list)
for r in ctl:
    if r["id"] in tphage and r["set"] in ("T1", "Tf") and r["skip"] == "1" and gkey(r) not in skip_genes:
        bygene[gkey(r)].append(r)
rng = random.Random(1)
dists = [int(r["n"]) - 1 - int(r["k"]) for r in skip]
ctrl = []
for g in sorted(bygene):
    d = rng.choice(dists)
    ctrl.append(min(bygene[g], key=lambda r: (abs(int(r["n"]) - 1 - int(r["k"]) - d), int(r["k"]))))
sets = [("SKIP", skip), ("CTRL", ctrl), ("RCGG", rcgg)]
need = {r["id"] for _, rs in sets for r in rs}

COMP = str.maketrans("ACGTN", "TGCAN")
CODE = dict(zip([a + b + c for a in "TCAG" for b in "TCAG" for c in "TCAG"],
                "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"))
seqs, name, buf = {}, None, []
for l in open(f"{CW}/data/phages.fna"):
    if l[0] == ">":
        if name in need:
            seqs[name] = "".join(buf).upper()
        name, buf = l[1:].split()[0], []
    else:
        buf.append(l.strip())
if name in need:
    seqs[name] = "".join(buf).upper()


def translate(s):
    aa = []
    for i in range(0, len(s) - 2, 3):
        a = CODE.get(s[i:i + 3], "X")
        if a == "*":
            break
        aa.append(a)
    return "".join(aa)


stat = Counter()
with open(f"{out}/fused.faa", "w") as fa, open(f"{out}/sites.tsv", "w") as ft:
    ft.write("sid\tset\tid\tquality\tcodon\tk\tn\tj\tlen_aa\n")
    for S, rs in sets:
        for i, r in enumerate(rs):
            g = seqs[r["id"]]
            L = len(g)
            b, e, st, k = int(r["begin"]), int(r["end"]), int(r["strand"]), int(r["k"])
            s = g if st == 1 else g.translate(COMP)[::-1]
            s0 = b - 1 if st == 1 else L - e
            if s[s0 + 3 * k:s0 + 3 * k + 3] != r["codon"]:
                stat[(S, "codon_mismatch")] += 1
                continue
            prot = translate(s[s0:s0 + 3 * k]) + translate(s[s0 + 3 * k + 1:])
            if len(prot) < k + 30:
                stat[(S, "short_plus1")] += 1
                continue
            sid = f"{S}_{i}"
            fa.write(f">{sid}\n{prot}\n")
            ft.write(f"{sid}\t{S}\t{r['id']}\t{grp[r['id']]['quality']}\t{r['codon']}\t{k}\t{r['n']}\t{k}\t{len(prot)}\n")
            stat[(S, "written")] += 1
for S, rs in sets:
    print(S, "rows", len(rs), {k[1]: v for k, v in stat.items() if k[0] == S})
