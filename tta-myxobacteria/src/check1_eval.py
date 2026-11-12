"""Check 1 (clade) verdicts from the added Arm B runs (data/check1/plan.tsv -> results/check1/armB.tsv), frozen
rule as in src/armB_gtdb_eval.py. Per candidate: added genomes run / passing X, same-species ones apart, and every
added genome that passes any codon.

usage: uv run --no-project python -I src/check1_eval.py  -> results/check1/verdict.txt
"""
import collections
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "results" / "check1"


def ok(r):
    return int(r["O_fam"]) >= 10 and float(r["ratio"]) >= 3 and float(r["log10p"]) <= -10


def main():
    plan = list(csv.DictReader(open(ROOT / "data" / "check1" / "plan.tsv"), delimiter="\t"))
    rows = collections.defaultdict(dict)
    for r in csv.DictReader(open(OUT / "armB.tsv"), delimiter="\t"):
        rows[r["label"]][r["codon"]] = r
    lines = []
    by = collections.defaultdict(list)
    for p in plan:
        by[(p["candidate"], p["genus"], p["codon"])].append(p)
    for (cand, g, x), ps in sorted(by.items(), key=lambda kv: kv[0][1]):
        run = [p for p in ps if x in rows.get(p["acc"], {})]
        hit = [p for p in run if ok(rows[p["acc"]][x])]
        same = [p for p in run if p["same_species"] == "1"]
        lines.append(f"{g} {x} (candidate {cand}): {len(hit)} of {len(run)} added genomes pass {x} "
                     f"({sum(p in hit for p in same)} of {len(same)} same species; {len(ps) - len(run)} not run) "
                     f"-> check 1 {'MET' if hit else 'not met'}")
        for p in run:
            r = rows[p["acc"]][x]
            lines.append(f"    {p['gtdb_acc']} {p['species']} cm2 {p['completeness']}/{p['contamination']} "
                         f"events {r['events']} {x} O {r['O']} O_fam {r['O_fam']} ratio {r['ratio']} "
                         f"log10p {r['log10p']}{' PASS' if ok(r) else ''}")
    other = [(lab, c, r) for lab, d in rows.items() for c, r in d.items() if ok(r)]
    lines.append(f"all passes in added genomes (any codon): {len(other)}")
    lines += [f"    {lab} {c} O_fam {r['O_fam']} ratio {r['ratio']} log10p {r['log10p']}" for lab, c, r in other]
    (OUT / "verdict.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
