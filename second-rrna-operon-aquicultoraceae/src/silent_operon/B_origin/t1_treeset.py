"""Assemble the T1 tree sets (23S and 16S): family + order copies (>= 60% of model), top outside hits of each query,
and controls/outgroups. Writes out/t1_<kind>_in.fna with short names  <id>__<genus-or-taxon>__<label>.
Usage: python -I src/t1_treeset.py
"""
import csv, subprocess
from collections import defaultdict
from Bio import SeqIO

tax = {}
with open("lab/work/frozen-core/results/lsu/scan23/gtdb_tax.tsv") as fh:
    for line in fh:
        p = line.rstrip("\n").split("\t"); tax[p[0]] = p[3]
loci = {r["id"]: r for r in csv.DictReader(open("out/loci.tsv"), delimiter="\t")}


def tops(fn, ssu, n, minL):
    hs = defaultdict(list)
    for line in open(fn):
        p = line.rstrip("\n").split("\t")
        q, s = p[0], p[1]
        hs[(q, s)].append((float(p[9]), int(p[3]), float(p[2]), int(p[4]), int(p[5]), p[10] if ssu else tax.get(s, "?")))
    best = defaultdict(list)
    for (q, s), lst in hs.items():
        lst.sort(reverse=True); took = []
        for h in lst:
            if all(h[4] < t[3] or h[3] > t[4] for t in took):
                took.append(h)
        L = sum(t[1] for t in took); pid = sum(t[1] * t[2] for t in took) / L
        t = took[0][5]
        if "Aquicultoraceae" in t or L < minL:
            continue
        best[q].append((pid, s, t))
    want = {}
    for q, lst in best.items():
        lst.sort(reverse=True)
        k = 0; seen_g = set()
        for pid, s, t in lst:
            g = t.split("g__")[1].split(";")[0] if "g__" in t else s
            if g in seen_g:
                continue
            seen_g.add(g); want[s] = t; k += 1
            if k >= n:
                break
    return want


def tag(t):
    p = t.split(" ")[1] if " " in t and t.startswith("RS_") is False and "d__" in t.split(" ")[-1] is False else t
    ph = t.split("p__")[1].split(";")[0] if "p__" in t else "?"
    cl = t.split("c__")[1].split(";")[0] if "c__" in t else "?"
    g = t.split("g__")[1].split(";")[0] if "g__" in t else "?"
    return f"{ph}-{cl}-{g}".replace(" ", "")


for kind, fn, ssu, minL, db in [("23S", "out/q23_vs_census.tsv", False, 1200, "lab/work/frozen-core/results/lsu/scan23/census_lsu_all.fna"),
                                ("16S", "out/q16_vs_ssuall.tsv", True, 800, "data/ssu_all_r232.fna")]:
    want = tops(fn, ssu, 4, minL)
    out = open(f"out/t1_{kind}_in.fna", "w")
    nm = 0
    allf = "out/lsu_all.fna" if kind == "23S" else "out/ssu_all.fna"
    for r in SeqIO.parse(allf, "fasta"):
        gid = r.id.split("|")[0]; L = loci[gid]
        if float(L["frac"]) < 0.6:
            continue
        fam = "AQ" if L["family"] == "Aquicultoraceae" else L["family"]
        out.write(f">{gid}__{fam}-{L['genus']}__{L['label'] if fam == 'AQ' else 'order'}\n{r.seq}\n"); nm += 1
    ids = set(want)
    for r in SeqIO.parse(db, "fasta"):
        key = r.id
        if key in ids:
            out.write(f">{key.split('~')[0]}__{tag(want[key])}__out\n{r.seq}\n"); nm += 1
    ctlf = "out/ctl_lsu.fna" if kind == "23S" else "out/ctl_ssu.fna"
    keep = {"23S": ["ec_23S_NC_000913.3_225761+", "vv_23S_NC_014965.1_50140+", "tb_23S_NC_014165.1_445262+",
                    "tt_23S_NZ_CP041345.1_533890+", "tt_23S_NZ_CP041345.1_1989407-"],
            "16S": ["ec_16S_NC_000913.3_223774+", "vv_16S_NC_014965.1_48092+", "tb_16S_NC_014165.1_443553+",
                    "tb_16S_NC_014165.1_2893291-", "tt_16S_NZ_CP041345.1_531966+", "tt_16S_NZ_CP041345.1_1992611-"]}[kind]
    for r in SeqIO.parse(ctlf, "fasta"):
        if r.id in keep:
            out.write(f">{r.id.replace('+','p').replace('-','m')}__ctl__ctl\n{r.seq}\n"); nm += 1
    out.close()
    print(kind, nm, "sequences;", len(want), "outside hits")
