"""Recount the README's numbers from the tables in data/.
usage: uv run python -I src/check_numbers.py   (from the folder root)
"""
import ast
import csv
import json
import re
from collections import Counter
from pathlib import Path

D = Path(__file__).resolve().parent.parent / "data"
READERS = {"GCA_040210545.2": "Tser", "GCA_040211365.2": "Tder", "GCA_040409785.2": "Ttur"}
OTHERS = {"GCA_040211755.2": "Tcar", "GCA_040210955.2": "Tduv", "GCA_037094695.1": "Tmey", "GCA_040212085.2": "AUS03"}
TCAR = "GCA_040211755.2"


def tsv(p):
    return list(csv.DictReader(open(p), delimiter="\t"))


print("== Codon usage (data/symbionts/symbiont_profile.tsv; pyrodigal genes under table 4)")
prof = {r["acc"]: r for r in tsv(D / "symbionts" / "symbiont_profile.tsv")}
for acc, n in {**READERS, **OTHERS, "GCA_029851105.1": "S_bost"}.items():
    r = prof.get(acc)
    if r:
        tga, tgg = float(r["TGA_per_k"]), float(r["TGG_per_k"])
        print(f"  {n:6s} TGA {tga:.2f} per 1,000 codons; TGA share of Trp codons {100 * tga / (tga + tgg):.1f}%")

print("== Orthologs: codon in the T. carinifrons symbiont's ortholog (data/symbionts/ortholog_tga/)")
for acc, n in READERS.items():
    rows = tsv(D / "symbionts" / "ortholog_tga" / f"{acc}_vs_{TCAR}.tsv")
    for own in ("TGA", "TGG"):
        sub = [r for r in rows if r["codonA"] == own]
        trp = sum(r["codonB"] in ("TGG", "TGA") for r in sub)
        print(f"  {n} {own} sites: {trp}/{len(sub)} at Trp in Tcar ({100 * trp / len(sub):.0f}%)")

print("== Raw reads (data/symbionts/read_support/)")
for f in sorted((D / "symbionts" / "read_support").glob("*.txt")):
    t = f.read_text()
    cov = int(re.search(r"in-frame TGA sites: \d+; covered by >= 3 reads: (\d+)", t).group(1))
    fixed = int(re.search(r"TGA in >= 90% of reads: (\d+)", t).group(1))
    print(f"  {f.stem}: TGA fixed at {fixed}/{cov} covered sites ({100 * fixed / cov:.1f}%)")

print("== Codetta (data/symbionts/codetta/*.json)")
for f in sorted((D / "symbionts" / "codetta").glob("*.json")):
    d = json.load(open(f))
    ns = d.get("nonstandard", [])
    if isinstance(ns, str):
        ns = ast.literal_eval(ns)  # the file stores a Python repr
    calls = [f"{x['codon']}={x['call']} (n {x['n']}, margin {x['margin_vs_std']:.1f})" for x in ns] or ["none"]
    print(f"  {d['acc']}: {', '.join(calls)}")

print("== TGA sites inside conserved alignment (data/throscidae/tga_sites.tsv)")
c = Counter((r["genome"], r["class"]) for r in tsv(D / "throscidae" / "tga_sites.tsv"))
for g in sorted({k[0] for k in c}):
    tot = sum(v for k, v in c.items() if k[0] == g)
    print(f"  {g}: {c[(g, 'inside alignment')]} of {tot} inside alignment")

print("== GTDB R232 QC (data/screen/qc_failed_throscid.tsv)")
for r in tsv(D / "screen" / "qc_failed_throscid.tsv"):
    acc = r["Accession"][3:]
    print(f"  {READERS.get(acc, OTHERS.get(acc, acc))}: markers {r['Markers (%)']}%, CM2 completeness {r['CM2 completeness (%)']}%")

print("== RF2 read scan: see data/rf2_readscan/summary.txt")
