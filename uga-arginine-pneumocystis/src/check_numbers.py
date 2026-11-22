"""Recount the README's numbers from the tables in data/.
usage: uv run python -I src/check_numbers.py   (from the folder root)
"""
import ast
import csv
from collections import Counter
from pathlib import Path

D = Path(__file__).resolve().parent.parent / "data"
ARG = {"AGA", "AGG", "CGA", "CGC", "CGG", "CGT"}


def tsv(p):
    return list(csv.DictReader(open(p), delimiter="\t"))


# Positions: UGA at strict Arg / strict Trp columns, UGG at strict Trp, UGA at Arg-majority columns
cols = {r["acc"]: r for r in tsv(D / "sites" / "pneumo_cols.tsv")}
j = cols["NC_020331.1"]
print(f"P. jirovecii: UGA at strict Arg {j['TGA_strictR']}, at strict Trp {j['TGA_strictW']}; "
      f"UGG at strict Trp {j['TGG_strictW']}; UGA placed {j['TGA_all']}, at Arg-majority columns {j['TGA_allR']}")
for acc, r in cols.items():
    if acc != "NC_020331.1":
        print(f"  {acc}: UGA at strict Arg {r['TGA_strictR']}, strict Trp {r['TGA_strictW']}")

# The other 8 genomes at the 23 columns
same = tsv(D / "sites" / "pneumo_same_cols.tsv")
c = Counter(v for r in same for k, v in r.items() if k != "col")
n = sum(c.values())
print(f"other 8 genomes at the {len(same)} columns: {sum(c[x] for x in ARG)} of {n} Arg codons, "
      f"UGA {c['TGA']}, UGG {c['TGG']}, other {dict((k, v) for k, v in c.items() if k not in ARG | {'TGA'})}")
tga_by = Counter(k for r in same for k, v in r.items() if v == "TGA")
print(f"  UGA by genome: {dict(tga_by)}  (NC_013660.2 P. carinii, NC_020332.1 P. murina)")

# Isolates (JX499143.1 is the reference itself and is not counted). The flank method misses one site that is
# UGA by whole-genome alignment (JX855938 atp8:41, '-' here); the README counts it.
iso = [r for r in tsv(D / "sites" / "isolates.tsv") if r["acc"] != "JX499143.1"]
ic = Counter(r["codon"] for r in iso)
print(f"isolates: {len({r['acc'] for r in iso})} genomes x 32 sites = {len(iso)}; flank method {dict(ic)}; "
      f"UGA with the missed site: {ic['TGA'] + ic['-']} of {len(iso)}")
print("  AGA sites:", Counter((r["gene"], r["codon_index"]) for r in iso if r["codon"] == "AGA"))

# RNA: 10 P. jirovecii runs, UGA sites at Arg-consensus columns
tot, same_reads, cga, sites = Counter(), 0, 0, set()
reads = 0
for f in sorted((D / "rna").glob("pjir_*.tsv")):
    for r in tsv(f):
        if r["codon"] == "TGA" and r["cons"] == "R":
            sites.add((r["gene"], r["codon_index"]))
            reads += int(r["reads"]); same_reads += int(r["same"])
            cga += ast.literal_eval(r["other"]).get("CGA", 0)
print(f"RNA, {len(sites)} sites: {same_reads} of {reads} reads UGA ({same_reads / reads:.1%}), CGA {cga}")
carinii = [r for r in tsv(D / "rna" / "pcar_SRR2156857.tsv") if r["codon"] == "TGA" and r["strict"] == "1" and r["cons"] == "R"]
print(f"P. carinii RNA, {len(carinii)} strict-Arg UGA sites: {sum(int(r['same']) for r in carinii)} of "
      f"{sum(int(r['reads']) for r in carinii)} reads UGA")
