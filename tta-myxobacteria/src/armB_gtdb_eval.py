"""GTDB Arm B under the frozen pass rule (PREREG.md: O_fam >= 10, ratio >= 3, log10p <= -10; the codon need not rank
first). Rows: results/gtdb/armB_<i>.tsv (labels via data/gtdb_armB/labels.tsv) and, for genomes reused from the
Coriobacteriales runs (data/gtdb_armB/reuse.tsv), results/corio/armB_{aa,ab,rerun,subst}.tsv.

Checks every selected genome has all its codon rows (61, or 62 under table 4) and was run with its selected table
(column `used`), and refuses otherwise (ALLOW_PARTIAL=1 evaluates the complete genomes anyway). Reports passes by
selection reason, whether each Arm A flag is confirmed (its codon passes in the same genome; flags with no Arm B
row are counted apart), genomes with < 20 events, passes at the rounding edge of the rule, and the most extreme
non-passing log10p. Fixes after the independent audit (scratch/audit_eval/): file order, completeness check,
reason dedup, untestable flags.

usage: uv run --no-project python -I src/armB_gtdb_eval.py  -> results/gtdb/armB_passes.tsv
"""
import collections
import csv
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
G = ROOT / "results" / "gtdb"
D = ROOT / "data" / "gtdb_armB"
MIN_EVENTS = 20


def ok(r):
    return int(r["O_fam"]) >= 10 and float(r["ratio"]) >= 3 and float(r["log10p"]) <= -10


def tsv(fn):
    return list(csv.DictReader(open(fn), delimiter="\t"))


def main():
    sel = {r["gtdb_acc"]: r for r in tsv(G / "armB_select.tsv")}
    label_of = {r["gtdb_acc"]: r["label"] for r in tsv(D / "labels.tsv")}
    for r in tsv(D / "reuse.tsv"):
        label_of["RS_" + r["acc"] if "RS_" + r["acc"] in sel else "GB_" + r["acc"]] = r["from"].split(":", 1)[1]
    rows = collections.defaultdict(dict)
    # Corio files first: a Corio genome rerun here with another table (armB_plan) must take the GTDB rows
    files = [ROOT / "results" / "corio" / f"armB_{p}.tsv" for p in ("aa", "ab", "rerun", "subst")]
    files += sorted(p for p in G.glob("armB_*.tsv") if p.name[5:-4].isdigit())
    for fn in files:
        by = collections.defaultdict(dict)
        for r in tsv(fn):
            by[r["label"]][r["codon"]] = r
        rows.update(by)  # later files replace a label's rows
    flags = collections.defaultdict(set)
    for r in tsv(G / "armA_flags.tsv"):
        flags[r["acc"]].add(r["codon"])

    norun = [a for a in sel if label_of.get(a) not in rows]
    used = {a: next(iter(rows[label_of[a]].values()))["used"] for a in sel if a not in norun}
    wrong = [a for a in used if used[a] != sel[a]["table"]]
    part = [a for a in used if len(rows[label_of[a]]) != 61 + (used[a] == "4")]  # every sense codon has a row
    print(f"{len(sel)} selected; {len(norun)} with no Arm B rows {norun[:8]}; {len(wrong)} run with another table "
          f"{wrong[:8]}; {len(part)} with a partial set of codon rows {part[:8]}")
    if (norun or wrong or part) and os.environ.get("ALLOW_PARTIAL") != "1":
        sys.exit("Arm B incomplete (ALLOW_PARTIAL=1 to evaluate anyway)")
    passes, low, border, worst = [], [], [], (0.0, None)
    by_reason = collections.Counter()
    conf = collections.Counter()
    for a, s in sel.items():
        if a in norun or a in part or a in wrong:
            continue
        cod = rows[label_of[a]]
        ev = int(next(iter(cod.values()))["events"])
        if ev < MIN_EVENTS:
            low.append(a)
        pc = sorted(c for c, r in cod.items() if ok(r))
        for c in pc:
            r = cod[c]
            passes.append((label_of[a], a, s["reasons"], s["table"], c, r["aa"], r["O"], r["O_fam"], r["E"], r["ratio"],
                           r["log10p"], r["rank"], str(ev), ",".join(sorted(flags[a])), s["taxonomy"]))
        for x in {y.split(":")[0] for y in s["reasons"].split(",")}:
            by_reason[(x, bool(pc))] += 1
        for c in flags[a]:  # a flag with no Arm B row (TGA in a suspect table-4 genome run with 11) is not testable
            conf["untestable" if c not in cod else c in pc] += 1
        for c in pc:
            r = cod[c]
            if float(r["log10p"]) > -10.01 or float(r["ratio"]) < 3.001:
                border.append(f"{a} {c} ratio {r['ratio']} log10p {r['log10p']}")
        for c, r in cod.items():
            if c not in pc and float(r["log10p"]) < worst[0]:
                worst = (float(r["log10p"]), f"{a} {c} O_fam {r['O_fam']} ratio {r['ratio']}")
    with open(G / "armB_passes.tsv", "w") as f:
        f.write("label\tgtdb_acc\treasons\ttable\tcodon\taa\tO\tO_fam\tE\tratio\tlog10p\trank\tevents\tarmA_flags\t"
                "taxonomy\n")
        f.writelines("\t".join(p) + "\n" for p in passes)
    print(f"genomes with < {MIN_EVENTS} events: {len(low)}")
    print(f"passing (genome, codon) pairs: {len(passes)} in {len({p[1] for p in passes})} genomes")
    for x in ("flag", "strat", "actino", "armc"):
        print(f"  {x}: {by_reason[(x, True)]} of {by_reason[(x, True)] + by_reason[(x, False)]} genomes pass")
    print(f"Arm A flags confirmed by Arm B: {conf[True]} of {conf[True] + conf[False]} testable "
          f"({conf['untestable']} flags have no Arm B row under the Arm B table)")
    print(f"passes at the rounding edge of the rule (output is rounded): {len(border)} {border[:10]}")
    print(f"most extreme non-passing log10p: {worst[0]:.2f} ({worst[1]})")
    for p in passes[:40]:
        print("  " + " ".join(p[1:13]) + " " + p[14].split(";")[1] + ";" + ";".join(p[14].split(";")[4:6]))


if __name__ == "__main__":
    main()
