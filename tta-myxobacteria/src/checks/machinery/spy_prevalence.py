"""How common is RF2 SPY across GTDB R232, and with which residue at E. coli position 201? (post hoc, fact-check)

Reads results/gtdb/rf_motifs.tsv (RF2 rows; the window holds 10 residues either side of the motif, so position 201
is 4 residues before the motif's S205) and the Arm A TTA flags.

usage: uv run --no-project python -I src/checks/machinery/spy_prevalence.py
"""
import collections
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
motifs, p201 = collections.defaultdict(set), collections.defaultdict(set)
for r in csv.DictReader(open(ROOT / "results/gtdb/rf_motifs.tsv"), delimiter="\t"):
    if r["class"] != "RF2":
        continue
    w = r["window"]
    i = w.find(r["residues"])
    motifs[r["acc"]].add(r["residues"])
    if i >= 4:
        p201[r["acc"]].add((w[i - 4], r["residues"]))
flags = {r["acc"] for r in csv.DictReader(open(ROOT / "results/gtdb/armA_flags.tsv"), delimiter="\t")
         if r["codon"] == "TTA"}

n = len(motifs)
print("genomes with an RF2 motif call:", n)
print("motifs:", collections.Counter(m for s in motifs.values() for m in s).most_common(6))
for m in ("SPF", "SPY"):
    gs = [a for a, s in motifs.items() if m in s]
    k = sum(a in flags for a in gs)
    print(f"{m}: {len(gs)} genomes ({len(gs) / n:.1%}); TTA-flagged in Arm A {k} ({k / len(gs):.3%})")
print("position 201 + motif:", collections.Counter(x for s in p201.values() for x in s
                                                    if x[1] in ("SPF", "SPY")).most_common(10))
for combo in (("F", "SPY"), ("F", "SPF")):
    gs = sorted(a for a, s in p201.items() if combo in s)
    print(f"F201 + {combo[1]}: {len(gs)} genomes; TTA-flagged {[a for a in gs if a in flags]}")
    if len(gs) < 20:
        print("  ", gs)
