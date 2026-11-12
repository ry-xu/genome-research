"""Coriobacteriales RF2 lead, clade contrast (PREREG.md, "Pre-registered lead: Coriobacteriales RF2", unit as
clarified 19:28): TGG lr (Arm A) and TGG ratio (Arm B) in Coriobacteriales vs the other Actinomycetota orders, one
genome per GTDB family on both sides (stratified-sample rule: most complete genes, ties by accession), one-sided
Mann-Whitney (Coriobacteriales greater) at p < 0.01, effect = difference of medians (Hodges-Lehmann shift also given).
Families whose genome has < 20 Arm B events are left out of the Arm B comparison and counted; genomes with < 300
complete genes have no Arm A score and are left out of the Arm A comparison and counted.
Reported, not tested: the same with every Coriobacteriales genome; a split by RF2 motif (SPT vs SP[FY] carriers)
across Actinomycetota families.

Arm B rows: results/corio/armB_{aa,ab,rerun,subst}.tsv (later files replace) and results/gtdb/armB_*.tsv via
data/gtdb_armB/labels.tsv; a RefSeq-twin label (GCF_x) stands for GTDB's GB_GCA_x; data/corio/substitutes.tsv maps
the six Coriobacteriales genomes NCBI serves no sequence for to their stand-ins (same length and scaffolds).

usage: uv run --no-project --with numpy,scipy python -I src/corio_contrast.py [A]   (A: Arm A side only)
Refuses the Arm B side while any family's genome has no Arm B rows (ALLOW_PARTIAL=1 overrides). Fixes after the
independent audit (scratch/audit_eval/): GTDB files are armB_<i>.tsv only; the refusal above.
"""
import collections
import csv
import gzip
import os
import sys
from pathlib import Path

import numpy as np
from scipy.stats import mannwhitneyu

ROOT = Path(__file__).resolve().parent.parent
G = ROOT / "results" / "gtdb"
MIN_EVENTS, ALPHA = 20, 0.01


def rank(tax, r):
    return next((x for x in tax.split(";") if x.startswith(r + "__")), "")


def genomes():
    out = {}
    for r in csv.DictReader(open(G / "armA_genomes.tsv"), delimiter="\t"):
        out[r["acc"]] = dict(n=int(r["ncomplete"]), tax=r["taxonomy"])
    return out


def family_picks(gen):
    pick = {}
    for a, v in gen.items():
        fam = rank(v["tax"], "f")
        if fam in ("", "f__"):
            continue
        b = pick.get(fam)
        if b is None or (-v["n"], a) < (-gen[b]["n"], b):
            pick[fam] = a
    return pick


def arm_a_tgg():
    lr = {}
    with gzip.open(G / "armA_rows.tsv.gz", "rt") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["codon"] == "TGG":
                lr[r["acc"]] = float(r["lr"])
    return lr


def arm_b_tgg():
    """GTDB accession -> (events, TGG ratio)."""
    lab = {}
    for lt in (ROOT / "data" / "gtdb_armB" / "labels.tsv", ROOT / "data" / "corio" / "substitutes.tsv"):
        if lt.exists():
            lab.update({r["label"]: r["gtdb_acc"] for r in csv.DictReader(open(lt), delimiter="\t")})
    out = {}

    def to_gtdb(label):
        return ("GB_GCA_" + label[4:]) if label.startswith("GCF_") else None

    files = [ROOT / "results" / "corio" / f"armB_{p}.tsv" for p in ("aa", "ab", "rerun", "subst")]
    files += sorted(p for p in G.glob("armB_*.tsv") if p.name[5:-4].isdigit())  # armB_<i>.tsv only
    for fn in files:
        if not fn.exists():
            continue
        for r in csv.DictReader(open(fn), delimiter="\t"):
            if r["codon"] != "TGG":
                continue
            l = r["label"]
            cands = [x for x in (lab.get(l), "RS_" + l, "GB_" + l, to_gtdb(l)) if x]  # a GTDB label wins
            out[l] = (cands, int(r["events"]), float(r["ratio"]))
    return out


def resolve(b, gen):
    res = {}
    for l, (cands, ev, ratio) in b.items():
        a = next((x for x in cands if x in gen), None)
        if a:
            res[a] = (ev, ratio)
    return res


def test(name, x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    p = mannwhitneyu(x, y, alternative="greater").pvalue
    hl = np.median(np.subtract.outer(x, y))
    print(f"{name}: Coriobacteriales n={len(x)} median {np.median(x):.3f}; other Actinomycetota n={len(y)} median "
          f"{np.median(y):.3f}; difference of medians {np.median(x) - np.median(y):+.3f} (HL {hl:+.3f}); "
          f"one-sided MW p = {p:.3g} -> {'PASS' if p < ALPHA else 'not significant'} at {ALPHA}")


def rf2_state(accs):
    st = collections.defaultdict(set)
    for r in csv.DictReader(open(G / "rf_motifs.tsv"), delimiter="\t"):
        if r["class"] != "RF2" or r["acc"] not in accs or float(r["ref_cover"]) < 0.6 or "-" in r["residues"]:
            continue
        st[r["acc"]].add("SPT" if r["residues"] == "SPT" else "SP[FY]" if r["canonical"] == "1" else "other")
    return {a: ",".join(sorted(s)) for a, s in st.items()}


def main(only_a=False):
    gen = genomes()
    pick = family_picks(gen)
    act = {f: a for f, a in pick.items() if rank(gen[a]["tax"], "p") == "p__Actinomycetota"}
    is_c = lambda a: rank(gen[a]["tax"], "o") == "o__Coriobacteriales"
    cfam = sorted(f for f, a in act.items() if is_c(a))
    ofam = sorted(f for f, a in act.items() if not is_c(a))
    print(f"Actinomycetota families {len(act)}: Coriobacteriales {len(cfam)}, other {len(ofam)}")
    lr = arm_a_tgg()
    xa = [lr[act[f]] for f in cfam if act[f] in lr]
    ya = [lr[act[f]] for f in ofam if act[f] in lr]
    print(f"Arm A: left out (< 300 complete genes) Coriobacteriales {len(cfam) - len(xa)}, other {len(ofam) - len(ya)}")
    test("Arm A TGG lr, one genome per family", xa, ya)
    call = [lr[a] for a in gen if is_c(a) and a in lr]
    test("Arm A TGG lr, every Coriobacteriales genome (reported, not tested)", call, ya)
    if only_a:
        return
    b = resolve(arm_b_tgg(), gen)
    nob = [f for f in cfam + ofam if act[f] not in b]
    print(f"Arm B: families with no Arm B run yet {len(nob)} {nob[:10]}")
    if nob and os.environ.get("ALLOW_PARTIAL") != "1":
        sys.exit("Arm B rows missing for some families (ALLOW_PARTIAL=1 to run the Arm B side anyway)")
    low = [f for f in cfam + ofam if act[f] in b and b[act[f]][0] < MIN_EVENTS]
    print(f"Arm B: left out (< {MIN_EVENTS} events) Coriobacteriales {sum(f in cfam for f in low)}, "
          f"other {sum(f in ofam for f in low)}")
    xb = [b[act[f]][1] for f in cfam if act[f] in b and b[act[f]][0] >= MIN_EVENTS]
    yb = [b[act[f]][1] for f in ofam if act[f] in b and b[act[f]][0] >= MIN_EVENTS]
    test("Arm B TGG ratio, one genome per family", xb, yb)
    callb = [v[1] for a, v in b.items() if is_c(a) and v[0] >= MIN_EVENTS]
    test("Arm B TGG ratio, every Coriobacteriales genome (reported, not tested)", callb, yb)
    st = rf2_state(set(act.values()))
    for arm, val in (("Arm A TGG lr", lambda a: lr.get(a)),
                     ("Arm B TGG ratio", lambda a: b[a][1] if a in b and b[a][0] >= MIN_EVENTS else None)):
        groups = collections.defaultdict(list)
        for f, a in act.items():
            v = val(a)
            if v is not None:
                groups[st.get(a, "none/missing")].append(v)
        print(f"{arm} by RF2 motif across Actinomycetota families (reported, not tested): " +
              "; ".join(f"{k} n={len(v)} median {np.median(v):.3f}" for k, v in sorted(groups.items())))


if __name__ == "__main__":
    main(len(sys.argv) > 1 and sys.argv[1] == "A")
