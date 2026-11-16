"""Per-genome readout for the frozen CGG skip test (README Gate, "CGG skip test, order TANB77").
usage: python -I src/analysis/cgg_test.py <set.tsv> <outdir> [filter]  -> TSV on stdout
filter: 'tanb77' (order TANB77 only), 'control' (class Clostridia, order not TANB77, GC <= 35%), or none.
cgg: clean +1 breaks whose old-frame codon is CGG; p1: all clean +1 breaks; ccg: an atlas trna row with anticodon CCG
(empty if the trna product is missing); call: cgg >= 3 and cgg / p1 >= 0.25."""
import csv
import gzip
import sys
from pathlib import Path

import duckdb

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from summarize import MIN_ALN, flag  # noqa: E402

from paths import PRODUCTS as ATLAS  # noqa: E402


def trna(gca):
    p = ATLAS / gca[:8] / gca / "trna.parquet"
    if not p.exists():
        return None
    return [a.upper() for (a,) in duckdb.sql(f"select anticodon from '{p}'").fetchall() if a]


def main():
    setf, outdir = sys.argv[1], Path(sys.argv[2])
    filt = sys.argv[3] if len(sys.argv) > 3 else None
    w = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    w.writerow(["gca", "family", "genus", "gc", "loci", "p1", "cgg", "share", "call", "ccg", "n_trna"])
    for r in csv.DictReader(open(setf), delimiter="\t"):
        gc = float(r["gc"])
        gc = gc / 100 if gc > 1 else gc
        if filt == "tanb77" and r["order"] != "TANB77":
            continue
        if filt == "control" and not (r["class"] == "Clostridia" and r["order"] != "TANB77" and gc <= 0.35):
            continue
        fsf = outdir / f"{r['gca']}.fs.tsv.gz"
        if not fsf.exists():
            continue
        fs = list(csv.DictReader(gzip.open(fsf, "rt"), delimiter="\t"))
        flag(fs)
        loci = sum(1 for h in csv.DictReader(gzip.open(outdir / f"{r['gca']}.hits.tsv.gz", "rt"), delimiter="\t")
                   if int(h["alen"]) >= MIN_ALN)
        p1 = [x for x in fs if x["clean"] and x["type"] == "+1"]
        cgg = sum(1 for x in p1 if x["old_codon"] == "CGG")
        call = cgg >= 3 and cgg / max(len(p1), 1) >= 0.25
        t = trna(r["gca"])
        w.writerow([r["gca"], r["family"], r["genus"], f"{gc:.3f}", loci, len(p1), cgg,
                    f"{cgg / len(p1):.2f}" if p1 else "", int(call),
                    "" if t is None else int("CCG" in t), "" if t is None else len(t)])


if __name__ == "__main__":
    main()
