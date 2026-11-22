"""T1 vertical inheritance: topology among genera of the divergent copies vs normal copies vs the GTDB species tree.
One leaf per genus; for every pair of trees on their shared genera, lists each split with support and whether the
other tree has it (conflict = the other tree has an incompatible split).
Usage: python -I src/t1_vertical.py
"""
import itertools
from Bio import Phylo

G = {"Aquicultor_A": "AqA", "CADDYS01": "CAD", "DASXVE01": "DXV", "Kmv33": "KMV", "DASQQM01": "DQQ", "JASEIR01": "JAS",
     "JAMCWL01": "JAM", "Aquicultor": "Aqu"}
tax = {}
for line in open("lab/work/frozen-core/results/lsu/scan23/gtdb_tax.tsv"):
    p = line.rstrip("\n").split("\t"); tax[p[0]] = p[3]


def lab_div(n):
    if n.endswith("__div"):
        return G.get(n.split("__")[1].split("-", 1)[1])
    return "JAM" if n.startswith("16S_024__") else None


def lab_norm(n):
    if n.endswith("__norm") and "__AQ-" in n and not n.startswith("23S_054"):
        return G.get(n.split("__")[1].split("-", 1)[1])
    return "JAM" if ("JAMCWL01" in n and n.endswith("__order") and not n.startswith("16S_024")) else None


AQG = {"Aquicultor_A", "CADDYS01", "DASXVE01", "Kmv33", "DASQQM01", "JASEIR01", "JAMCWL01", "Aquicultor"}
def lab_gtdb(n):
    # sp_labels tree leaves are species names, e.g. 's__DASXVE01 sp036504445'
    n = (n or "").strip("'")
    if not n.startswith("s__"):
        return None
    g = n[3:].split(" ")[0]
    return G.get(g) if g in AQG else None


def gsplits(t, lab, support_from_name=False):
    leaves = {}
    for l in t.get_terminals():
        g = lab(l.name)
        if g and g not in leaves:
            leaves[g] = l
    under = {}
    def rec(c):
        s = set()
        if c.is_terminal():
            for g, l in leaves.items():
                if l is c:
                    s.add(g)
        for ch in c.clades:
            s |= rec(ch)
        under[id(c)] = (c, frozenset(s))
        return s
    rec(t.root)
    out = {}
    for c, s in under.values():
        if c.is_terminal():
            continue
        sup = c.confidence
        if support_from_name and c.name:
            try:
                sup = float(str(c.name).split(":")[0])
            except ValueError:
                sup = None
        out.setdefault(s, sup)
    return set(leaves), out


def restrict(sp, shared):
    shared = frozenset(shared); res = {}
    for s, sup in sp.items():
        x = s & shared
        if 1 < len(x) < len(shared) - 1:
            key = min(x, shared - x, key=lambda z: sorted(z))
            if key not in res or (sup or 0) > (res[key] or 0):
                res[key] = sup
    return res


def compatible(a, b, shared):
    A, B = a, b; Ac, Bc = shared - a, shared - b
    return not (A & B and A & Bc and Ac & B and Ac & Bc)


T = {}
for k in ("23S", "16S"):
    t = Phylo.read(f"out/t1_{k}.contree", "newick")
    T[f"div{k}"] = gsplits(t, lab_div)
    T[f"norm{k}"] = gsplits(t, lab_norm)
gt = Phylo.read("lab/work/frozen-core/data/gtdb/bac120_r232.sp_labels.tree", "newick")
T["GTDB"] = gsplits(gt, lab_gtdb, support_from_name=True)
for nm, (tx, sp) in T.items():
    print(f"{nm:8s} genera {sorted(tx)}")
for a, b in [("div23S", "norm23S"), ("div23S", "GTDB"), ("div16S", "norm16S"), ("div16S", "GTDB"), ("norm23S", "GTDB"),
             ("norm16S", "GTDB"), ("div23S", "div16S")]:
    shared = frozenset(T[a][0] & T[b][0])
    sa, sb = restrict(T[a][1], shared), restrict(T[b][1], shared)
    print(f"\n{a} vs {b} on {sorted(shared)}")
    for s, sup in sa.items():
        if s in sb:
            print(f"   {sorted(s)} ({sup}) shared ({sb[s]})")
        else:
            bad = [(sorted(x), sb[x]) for x in sb if not compatible(s, x, shared)]
            print(f"   {sorted(s)} ({sup}) {'CONFLICT with ' + str(bad) if bad else 'not resolved in ' + b}")
