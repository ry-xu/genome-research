"""T3: where along the molecule do divergent copies change, compared with functional copies at matched divergence.
Usage: python -I src/t3_profile.py 23S|16S   -> out/t3_<kind>_pairs.tsv, out/t3_<kind>_regions.tsv, out/t3_<kind>_sites.tsv
"""
import sys, itertools, statistics as st
sys.path.insert(0, "src")
from alnlib import read_sto, read_cols, short

kind = sys.argv[1]
aln = read_sto("out/aln23.sto" if kind == "23S" else "out/aln16.sto")
cols = read_cols(kind)
names = {short(n) if "|" in n else n: n for n in aln}
A = {k: aln[v] for k, v in names.items()}

if kind == "23S":
    DIV = {"AqA": "23S_004", "CADa": "23S_042", "CADb": "23S_062", "DQQ": "23S_052", "DXV": "23S_055", "KMV": "23S_016"}
    OWN = {"AqA": "23S_003", "CADa": "23S_043", "DQQ": "23S_051", "DXV": "23S_053"}
    FAMN = ["23S_003", "23S_023", "23S_043", "23S_051", "23S_053", "23S_058", "23S_021"]
    REG = {"I": (1, 579), "II": (580, 1261), "III": (1262, 1646), "IV": (1647, 2010), "V": (2011, 2625), "VI": (2626, 2904)}
    SITES = {"PTC": [(2057, 2063), (2250, 2254), (2439, 2456), (2496, 2507), (2550, 2557), (2571, 2589), (2599, 2611)],
             "P-loop H80": [(2248, 2256)], "A-loop H92": [(2548, 2558)], "SRL H95": [(2653, 2667)],
             "GAC H43-44": [(1051, 1108)], "H69": [(1906, 1924)], "A-site finger H38": [(865, 905)],
             "L1 stalk H76-78": [(2093, 2195)], "H89": [(2457, 2495)], "H1 ends": [(1, 10), (2895, 2904)],
             "tunnel 2058-2062/750s": [(745, 765), (2057, 2062), (2609, 2611)]}
else:
    DIV = {"AqA": "16S_004", "CADa": "16S_030", "CADb": "16S_051", "DQQ": "16S_043", "DXV": "16S_045", "KMV": "16S_015"}
    OWN = {"AqA": "16S_003", "DXV": "16S_044"}
    FAMN = ["16S_003", "16S_020", "16S_026", "16S_027", "16S_044", "16S_049", "16S_000"]
    REG = {"5'": (1, 556), "central": (557, 912), "3'major": (913, 1396), "3'minor": (1397, 1542)}
    SITES = {"decoding h44": [(1400, 1410), (1490, 1500)], "530 loop h18": [(515, 535)], "790 loop h24": [(785, 795)],
             "P-site 1338/h23 690s": [(690, 700), (1336, 1342)], "h34 1050/1200": [(1046, 1067), (1189, 1211)],
             "h45 1516-1519": [(1506, 1529)], "anti-SD tail": [(1530, 1542)], "h28 neck": [(921, 933), (1384, 1396)],
             "h1 5' end": [(1, 20)]}

epos = [r["ecoli_pos"] for r in cols]
fcons = [r["f_cons"] for r in cols]
cons = [r["cons"].replace("T", "U") for r in cols]
ok_col = [r["f_gap"] < 0.2 for r in cols]


def in_ranges(p, rr):
    return any(a <= p <= b for a, b in rr)


def compare(x, y, sel=None):
    n = d = 0
    for i, (a, b) in enumerate(zip(x, y)):
        if not ok_col[i] or a in "-." or b in "-." or a == "N" or b == "N":
            continue
        if sel is not None and not sel(i):
            continue
        n += 1; d += a != b
    return n, d


BINS = [("inv>=0.99", 0.99, 1.01), ("0.95-0.99", 0.95, 0.99), ("0.80-0.95", 0.80, 0.95), ("var<0.80", 0.0, 0.80)]


def pair_stats(x, y):
    n, d = compare(x, y)
    out = {"n": n, "pid": 1 - d / n if n else float("nan")}
    for lab, lo, hi in BINS:
        bn, bd = compare(x, y, lambda i: lo <= fcons[i] < hi)
        out[lab] = (bn, bd)
    ni, di = out["inv>=0.99"]
    # invariant-site depletion: share of differences at invariant columns / share of compared columns that are invariant
    out["inv_depl"] = (di / d) / (ni / n) if d and n and ni else float("nan")
    out["reg"] = {}
    for r, (a, b) in REG.items():
        out["reg"][r] = compare(x, y, lambda i: a <= epos[i] <= b)
    out["site"] = {}
    for s, rr in SITES.items():
        out["site"][s] = compare(x, y, lambda i: in_ranges(epos[i], rr))
    return out


pairs = []
dv = list(DIV.items())
for (g1, a), (g2, b) in itertools.combinations(dv, 2):
    if a in A and b in A and not (g1.startswith("CAD") and g2.startswith("CAD")):
        pairs.append(("div-div", f"{g1}:{g2}", A[a], A[b]))
for g, a in DIV.items():
    if g in OWN and a in A:
        pairs.append(("div-ownnorm", g, A[a], A[OWN[g]]))
fam = [n for n in FAMN if n in A]
for a, b in itertools.combinations(fam, 2):
    pairs.append(("norm-norm fam", f"{a}:{b}", A[a], A[b]))
# order-level normals (other families) and functional controls
other = [k for k, v in names.items() if "|" in v and v.endswith("|norm") and k not in FAMN
         and "Aquicultoraceae" not in v and sum(c not in "-." for c in A[k]) > 0.8 * len(A[k])]
fun = other + [k for k in A if k.startswith(("ec_", "vv_")) ] [:0]
reps_ctl = {}
for k in A:
    if k[:3] in ("ec_", "vv_", "tb_", "tt_"):
        reps_ctl.setdefault(k[:2], k)
famrep = ["23S_053" if kind == "23S" else "16S_044", "23S_003" if kind == "23S" else "16S_003"]
for a, b in itertools.combinations(other + famrep + list(reps_ctl.values()), 2):
    pairs.append(("functional", f"{a}:{b}", A[a], A[b]))
# divergent-operon controls inside one genome
ctlpairs = {"23S": [("tt div-norm", "tt_23S_NZ_CP041345.1_1989407-", "tt_23S_NZ_CP041345.1_533890+"),
                    ("Arsenophonus decaying-norm", "ar_23S_JBQMRP010000030.1_14720+", "ar_23S_JBQMRP010000004.1_315056+"),
                    ("Arsenophonus2 decaying-norm", "as_23S_JBHEQW010000092.1_13747+", "as_23S_JBHEQW010000001.1_72701-")],
            "16S": [("tt div-norm", "tt_16S_NZ_CP041345.1_1992611-", "tt_16S_NZ_CP041345.1_531966+"),
                    ("tb typeB-typeA", "tb_16S_NC_014165.1_2893291-", "tb_16S_NC_014165.1_443553+"),
                    ("Arsenophonus decaying-norm", "ar_16S_JBQMRP010000030.1_13534+", "ar_16S_JBQMRP010000004.1_313238+")]}
for lab, a, b in ctlpairs[kind]:
    if a in A and b in A:
        pairs.append((lab, f"{a}:{b}", A[a], A[b]))

res = []
for cat, lab, x, y in pairs:
    s = pair_stats(x, y)
    if s["n"] < 0.5 * len(cols):
        continue
    res.append((cat, lab, s))

with open(f"out/t3_{kind}_pairs.tsv", "w") as fh:
    fh.write("cat\tpair\tn\tpid\t" + "\t".join(f"rate_{b[0]}" for b in BINS) + "\tinv_depl\t" +
             "\t".join(f"rel_{r}" for r in REG) + "\n")
    for cat, lab, s in res:
        tot = 1 - s["pid"]
        rr = []
        for r in REG:
            n, d = s["reg"][r]
            rr.append(f"{(d / n) / tot:.2f}" if n and tot else "nan")
        fh.write(f"{cat}\t{lab}\t{s['n']}\t{s['pid']:.3f}\t" +
                 "\t".join(f"{s[b[0]][1] / s[b[0]][0]:.4f}" if s[b[0]][0] else "nan" for b in BINS) +
                 f"\t{s['inv_depl']:.3f}\t" + "\t".join(rr) + "\n")

# summary by category
print(f"== {kind}: per category (median over pairs)")
cats = []
for c, _, _ in res:
    if c not in cats:
        cats.append(c)
for c in cats:
    ss = [s for cc, _, s in res if cc == c]
    pid = [s["pid"] for s in ss]
    dep = [s["inv_depl"] for s in ss if s["inv_depl"] == s["inv_depl"]]
    rinv = [s["inv>=0.99"][1] / s["inv>=0.99"][0] for s in ss]
    rvar = [s["var<0.80"][1] / s["var<0.80"][0] for s in ss]
    print(f"{c:28s} n={len(ss):3d} pid {min(pid):.3f}-{max(pid):.3f}  rate_inv {st.median(rinv):.4f}  rate_var {st.median(rvar):.3f}"
          f"  inv_depl median {st.median(dep) if dep else float('nan'):.3f}")
