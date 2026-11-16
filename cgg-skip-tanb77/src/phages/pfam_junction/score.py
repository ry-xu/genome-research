"""Score junction-spanning Pfam domains per set (RULES.md, frozen 18:43; Amendment 1).
usage: python -I src/score.py SITES.tsv SAMPLE.faa DOMTBL OUT.tsv"""
import csv, re, sys
from collections import Counter, defaultdict

sites_p, faa_p, dom_p, out_p = sys.argv[1:5]
PH = "data/pfam/phage.hmm"
HALL = re.compile(r"capsid|portal|terminase|tail|baseplate|head|prohead|scaffold", re.I)
phage, name = {}, None
for l in open(PH):
    if l.startswith("NAME "):
        name = l.split()[1]
    elif l.startswith("DESC "):
        phage[name] = l[6:].strip()
site = {r["sid"]: r for r in csv.DictReader(open(sites_p), delimiter="\t")}
sample = [l[1:].strip() for l in open(faa_p) if l[0] == ">"]
span = defaultdict(list)
for l in open(dom_p):
    if l[0] == "#":
        continue
    f = l.split()
    sid, pf, iev, a0, a1 = f[0], f[3], float(f[12]), int(f[17]), int(f[18])
    j = int(site[sid]["j"]) + 1  # 1-based junction residue (first +1-frame residue)
    if iev <= 1e-5 and a0 <= j - 10 and a1 >= j + 10:
        span[sid].append(pf)
n, s, h, sq, nq = Counter(), Counter(), Counter(), Counter(), Counter()
with open(out_p, "w") as fo:
    fo.write("sid\tset\tid\tquality\tcodon\tspanned\thallmark\tdomains\n")
    for sid in sample:
        r = site[sid]
        S, q = r["set"], r["quality"]
        doms = sorted(set(span[sid]))
        hall = [d for d in doms if d in phage and (HALL.search(d) or HALL.search(phage[d]))]
        n[S] += 1; nq[(S, q)] += 1
        s[S] += bool(doms); sq[(S, q)] += bool(doms); h[S] += bool(hall)
        fo.write(f"{sid}\t{S}\t{r['id']}\t{q}\t{r['codon']}\t{int(bool(doms))}\t{','.join(hall)}\t{','.join(doms)}\n")
for S in ("SKIP", "CTRL", "RCGG"):
    if n[S]:
        print(f"{S}\tsites {n[S]}\tSPANNED {s[S]} ({s[S] / n[S]:.3f})\tHALLMARK {h[S]}")
        for q in ("Complete", "High-quality", "Medium-quality"):
            if nq[(S, q)]:
                print(f"  {q}\t{nq[(S, q)]}\tSPANNED {sq[(S, q)]} ({sq[(S, q)] / nq[(S, q)]:.3f})")
