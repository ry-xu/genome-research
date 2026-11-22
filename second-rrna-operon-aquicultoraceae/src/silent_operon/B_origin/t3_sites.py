"""T3 sites: per named functional site, census-conserved columns (f_cons >= 0.95) that deviate from the census consensus
in each divergent copy, against expectation from that copy's molecule-wide deviation rate at the same columns class.
Also base-pair maintenance: among census-paired columns, fraction canonical (WC/GU) per sequence, and for div-div pairs,
the share of paired-column changes that keep both pairs canonical.
Usage: python -I src/t3_sites.py 23S|16S
"""
import sys, math, itertools
sys.path.insert(0, "src")
from alnlib import read_sto, read_cols, short

kind = sys.argv[1]
exec(open("src/t3_profile.py").read().split("epos = ")[0].split("aln = read_sto")[0])  # imports only
aln = read_sto("out/aln23.sto" if kind == "23S" else "out/aln16.sto")
cols = read_cols(kind)
A = {(short(n) if "|" in n else n): s for n, s in aln.items()}
src = open("src/t3_profile.py").read()
blk = src[src.index('if kind == "23S":'):src.index("epos = ")]
exec(blk)

epos = [r["ecoli_pos"] for r in cols]; f = [r["f_cons"] for r in cols]
cons = [r["cons"].replace("T", "U") for r in cols]; part = [r["partner"] for r in cols]
ctl = {"23S": {"tt_div": "tt_23S_NZ_CP041345.1_1989407-", "tt_norm": "tt_23S_NZ_CP041345.1_533890+",
               "Ars_decay": "ar_23S_JBQMRP010000030.1_14720+", "Ars_norm": "ar_23S_JBQMRP010000004.1_315056+",
               "Ecoli": "ec_23S_NC_000913.3_225761+"},
       "16S": {"tt_div": "tt_16S_NZ_CP041345.1_1992611-", "tt_norm": "tt_16S_NZ_CP041345.1_531966+",
               "tb_typeB": "tb_16S_NC_014165.1_2893291-", "Ars_decay": "ar_16S_JBQMRP010000030.1_13534+",
               "Ecoli": "ec_16S_NC_000913.3_223774+"}}[kind]
seqs = {**{f"div_{g}": A[k] for g, k in DIV.items() if k in A}, **{f"norm_{g}": A[k] for g, k in OWN.items() if k in A},
        **{k: A[v] for k, v in ctl.items() if v in A}}

CONS_T = 0.95


def poisson_cdf(k, lam):
    return sum(math.exp(-lam) * lam ** i / math.factorial(i) for i in range(k + 1))


def poisson_sf(k, lam):
    return 1 - poisson_cdf(k - 1, lam) if k > 0 else 1.0


print(f"== {kind} sites: deviations from census consensus at columns with f_cons >= {CONS_T}")
rates = {}
for nm, s in seqs.items():
    n = d = 0
    for i in range(len(cols)):
        if f[i] >= CONS_T and s[i] not in "-.N":
            n += 1; d += s[i] != cons[i]
    rates[nm] = d / n
print("molecule-wide deviation rate at conserved columns: " + ", ".join(f"{k} {v:.3f}" for k, v in rates.items()))
hdr = "site".ljust(24) + "ncol " + " ".join(f"{k[:9]:>9s}" for k in seqs) + "   div_sum obs/exp  P(kept) P(drift)"
print(hdr)
out = open(f"out/t3_{kind}_sites.tsv", "w")
out.write("site\tncol\t" + "\t".join(seqs) + "\tdiv_obs\tdiv_exp\tp_kept\tp_drift\n")
for site, rr in SITES.items():
    idx = [i for i in range(len(cols)) if any(a <= epos[i] <= b for a, b in rr) and f[i] >= CONS_T]
    cells, dobs, dexp = [], 0, 0.0
    for nm, s in seqs.items():
        cov = [i for i in idx if s[i] not in "-.N"]
        dev = sum(s[i] != cons[i] for i in cov)
        gaps = len(idx) - len(cov)
        cells.append(f"{dev}/{len(cov)}" + ("" if not gaps else ""))
        if nm.startswith("div_"):
            dobs += dev; dexp += rates[nm] * len(cov)
    pk, pd = poisson_cdf(dobs, dexp), poisson_sf(dobs, dexp)
    print(f"{site[:24]:24s}{len(idx):4d} " + " ".join(f"{c:>9s}" for c in cells) + f"   {dobs:3d}/{dexp:5.1f}  {pk:.3g}  {pd:.3g}")
    out.write(f"{site}\t{len(idx)}\t" + "\t".join(cells) + f"\t{dobs}\t{dexp:.2f}\t{pk:.3g}\t{pd:.3g}\n")

# base pairing
CAN = {"AU", "UA", "GC", "CG", "GU", "UG"}
pairs = [(i, part[i]) for i in range(len(cols)) if part[i] > i and f[i] >= 0.0]
print("\nfraction of census base pairs canonical (both bases present):")
for nm, s in seqs.items():
    n = sum(1 for i, j in pairs if s[i] not in "-.N" and s[j] not in "-.N")
    c = sum(1 for i, j in pairs if s[i] not in "-.N" and s[j] not in "-.N" and s[i] + s[j] in CAN)
    print(f"  {nm:12s} {c}/{n} = {c/n:.3f}")
print("\ndiv-div: paired columns where the two copies differ; kept canonical in both / changed pairs; double changes")
tot = [0, 0, 0, 0]
for (g1, a), (g2, b) in itertools.combinations(DIV.items(), 2):
    if a not in A or b not in A or (g1.startswith("CAD") and g2.startswith("CAD")):
        continue
    x, y = A[a], A[b]
    ch = keep = dbl = dblkeep = 0
    for i, j in pairs:
        if any(c in "-.N" for c in (x[i], x[j], y[i], y[j])):
            continue
        if x[i] != y[i] or x[j] != y[j]:
            ch += 1
            k = (x[i] + x[j] in CAN) and (y[i] + y[j] in CAN)
            keep += k
            if x[i] != y[i] and x[j] != y[j]:
                dbl += 1; dblkeep += k
    tot = [tot[0] + ch, tot[1] + keep, tot[2] + dbl, tot[3] + dblkeep]
print(f"  all div-div pairs: changed pairs {tot[0]}, both canonical {tot[1]} ({tot[1]/tot[0]:.2f}); double changes {tot[2]}, of which canonical {tot[3]} ({tot[3]/max(tot[2],1):.2f})")
