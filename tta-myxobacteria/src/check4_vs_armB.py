"""WITHDRAWN 2026-11-11 (independent audit, scratch/audit_check4/REPORT.md item 3): not like-for-like. The
observed O includes TTA after a canonical stop at e*+1, which the prediction does not, and the denominator
(all hits without a stop at e*+1) is about 5x the hits that are Arm B events. Do not cite its output. Kept for the
record.

Post-unblinding comparison (not a frozen check): does the check-4 TTA count in each TTA candidate match what its
Arm B rate predicts? Arm B rate = share of the genome's Arm B events with TTA at e*+1..e*+3. Predicted check-4 O =
that rate x the genome's tblastn hits without a canonical stop at e*+1, i.e. the hits that run past the homolog end,
as every Arm B event does. Totals per clade are compared by Poisson. Also lists the families that carry TTA at
e*+1..e*+3 in the most genomes.

usage: uv run --no-project --with scipy python -I src/check4_vs_armB.py  -> results/check4_vs_armB.txt
"""
import collections
import csv
from pathlib import Path

from scipy.stats import poisson

ROOT = Path(__file__).resolve().parent.parent
POLY = ("DASTPP01", "JADLFU01", "JANWXX01", "JAXLYQ01")
STOPS = {"TAA", "TAG", "TGA"}


def tsv(fn):
    return list(csv.DictReader(open(fn), delimiter="\t"))


def main():
    key = {r["gtdb_acc"]: r["id"] for r in tsv(ROOT / "results" / "check4_key.tsv")}
    cand = {r["gtdb_acc"]: r for r in tsv(ROOT / "results" / "gtdb" / "candidates.tsv") if r["codon"] == "TTA"}
    clade = {a: "Polyangiaceae" if c["genus"][3:] in POLY else "Streptosporangiaceae" for a, c in cand.items()}
    ev = collections.defaultdict(list)
    for fn in ("armB_0_events.tsv", "armB_1_events.tsv"):
        for r in tsv(ROOT / "results" / "gtdb" / fn):
            ev[r["label"]].append(r)
    lines, tot = [], collections.defaultdict(lambda: [0, 0.0])
    fam = collections.defaultdict(lambda: collections.defaultdict(set))
    for acc, c in sorted(cand.items(), key=lambda kv: (clade[kv[0]], kv[1]["genus"], kv[0])):
        es = ev[acc[3:]]
        tta = [e for e in es if "TTA" in e["window"].split(",")[2:]]
        for e in tta:
            fam[clade[acc]][e["family"]].add(acc)
        rate = len(tta) / len(es)
        hits = tsv(ROOT / "check4" / "work" / "hits" / f"{key[acc]}.hits.tsv")
        over = sum(h["w1"] not in STOPS for h in hits)
        o = sum("TTA" in (h["w1"], h["w2"], h["w3"]) for h in hits)
        tot[clade[acc]][0] += o
        tot[clade[acc]][1] += rate * over
        lines.append(f"{clade[acc]:20s} {c['genus'][3:]:10s} {acc} Arm B {len(tta)}/{len(es)} = {rate:.3f}; "
                     f"check 4 hits past the end {over}/{len(hits)}; predicted O {rate * over:.1f}, observed {o}")
    for cl, (o, e) in sorted(tot.items()):
        lines.append(f"{cl}: observed {o}, predicted {e:.1f}, P(X <= observed) {poisson.cdf(o, e):.2g}")
    for cl, d in sorted(fam.items()):
        top = sorted(d.items(), key=lambda kv: (-len(kv[1]), kv[0]))[:10]
        lines.append(f"{cl}: families with TTA at e*+1..e*+3 in most genomes: "
                     + ", ".join(f"{f} {len(g)}" for f, g in top))
    (ROOT / "results" / "check4_vs_armB.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
