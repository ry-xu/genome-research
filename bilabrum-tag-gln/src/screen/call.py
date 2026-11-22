"""Frozen stop->sense calls (README "Rules", frozen 2026-11-11 21:17 from the controls only).

usage: call.py SCORE.tsv > calls.tsv      (SCORE.tsv from score.py)
Per unit and stop codon:
  call   n >= 20 inner aligned columns, in >= 10 distinct reference proteins; rate_al >= 0.20 per 1,000 inner aligned
         columns; the profile meaning aligns there in >= 5 distinct reference proteins; rel >= 0.10.
  watch  not a call, but n >= 10, >= 5 proteins, rate_al >= 0.05, meaning in >= 3 proteins, rel >= 0.10.
  -      anything else.
Per unit: powered = inner >= 100,000 (a codon used at 0.2 per 1,000 then gives >= 20 aligned stops).
TGA aligned to Sec (U) is excluded from n and only reported (sec).
"""
import csv
import sys

CALL = dict(n=20, prot=10, rate=0.20, pm=5, rel=0.10)
WATCH = dict(n=10, prot=5, rate=0.05, pm=3, rel=0.10)
POWER = 100_000


def tier(r):
    if not r["rel"]:
        return "-"
    n, prot, rate, pm, rel = int(r["n"]), int(r["prot"]), float(r["rate_al"]), int(r["prot_meaning"]), float(r["rel"])
    for name, t in (("call", CALL), ("watch", WATCH)):
        if n >= t["n"] and prot >= t["prot"] and rate >= t["rate"] and pm >= t["pm"] and rel >= t["rel"]:
            return name
    return "-"


def main(path):
    rows = list(csv.DictReader(open(path), delimiter="\t"))
    cols = list(rows[0].keys()) + ["powered", "tier"]
    print("\t".join(cols))
    for r in rows:
        r["powered"] = "yes" if int(r["inner"]) >= POWER else "no"
        r["tier"] = tier(r)
        print("\t".join(str(r[c]) for c in cols))


if __name__ == "__main__":
    main(sys.argv[1])
