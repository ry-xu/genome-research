"""Recount the README's numbers from the tables in ../data/ (and the two GenBank records).
usage (from the folder root): uv run python -I src/check_numbers.py

Census numbers come from tables written by the census code (not included); the isolate pairing is recomputed here
from the two GenBank records with the same method as src/newleads/codon_pair.py.
"""
import ast
import csv
import os
from collections import Counter

from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.Seq import Seq
from scipy.stats import fisher_exact

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
DATA = os.path.join(ROOT, "data")
ARG = {"AGA", "CGT", "CGC", "CGA", "CGG"}  # AGG left out: it is the codon in question


def rows(path):
    with open(os.path.join(DATA, path)) as f:
        return list(csv.DictReader(f, delimiter="\t"))


# 1. Red algal census (data/census/rhodo_agg.tsv, from src/newleads/codon_by_clade.py)
red = rows("census/rhodo_agg.tsv")
met = {r["acc"]: int(r["AGG@M"]) for r in red if int(r["AGG@M"])}
print(f"red algal genomes: {len(red)}; AGG at strict Met: {sum(met.values())} {met}; "
      f"AGG at strict Arg: {sum(int(r['AGG@R']) for r in red)}")
sag = next(r for r in red if r["acc"] == "NC_062302.1")
print(f"SAG 1380-2 strict Met: ATG {sag['ATG@M']}, AGG {sag['AGG@M']}; strict Arg: AGG {sag['AGG@R']}, "
      f"AGA {sag['AGA@R']}")

# 2. Relatives at SAG 1380-2's AGG columns (data/census/porph_agg_columns.tsv, from src/newleads/same_column.py)
cols = rows("census/porph_agg_columns.tsv")
rel = [k for k in cols[0] if k not in ("model:col", "cons", "strict")]
cod = [c for r in cols for k in rel for c in r[k].split(",") if c not in ("-", "")]
print(f"relatives ({len(rel)}) at SAG's {len(cols)} AGG columns: ATG {sum(c == 'ATG' for c in cod)} of {len(cod)}, "
      f"Arg codon {sum(c in ARG for c in cod)}")

# 3. The two isolates (recomputed from data/genomes/, as codon_pair.py does)
def cds(path):
    rec = SeqIO.read(path, "genbank")
    return {(f.qualifiers.get("gene") or f.qualifiers.get("product") or ["?"])[0].lower(): str(f.extract(rec.seq)).upper()
            for f in rec.features if f.type == "CDS"}


a = cds(os.path.join(DATA, "genomes/NC_062302.1.gb"))
b = cds(os.path.join(DATA, "genomes/MT483995.1.gb"))
al = PairwiseAligner(); al.substitution_matrix = substitution_matrices.load("BLOSUM62")
al.open_gap_score, al.extend_gap_score, al.mode = -10, -0.5, "global"
tr = lambda nt: str(Seq(nt[:len(nt) // 3 * 3]).translate(table=4)).replace("*", "X")
pairs, same = [], 0
for g in sorted(set(a) & set(b)):
    aln = al.align(tr(a[g]), tr(b[g]))[0]
    for (s1, e1), (s2, _) in zip(*aln.aligned):
        for k in range(e1 - s1):
            i, j = s1 + k, s2 + k
            pairs.append((a[g][3 * i:3 * i + 3], b[g][3 * j:3 * j + 3]))
print(f"shared genes: {len(set(a) & set(b))}; same codon at {sum(x == y for x, y in pairs)} of {len(pairs)} aligned "
      f"positions ({100 * sum(x == y for x, y in pairs) / len(pairs):.0f}%)")
agg = Counter("AGG" if y == "AGG" else "ATG" if y == "ATG" else "Arg" if y in ARG else "other"
              for x, y in pairs if x == "AGG")
argc = Counter("same" if y == x else "other Arg" if y in ARG else "ATG/AGG" if y in ("ATG", "AGG") else "other"
               for x, y in pairs if x in ARG)
print(f"at SAG's AGG ({sum(agg.values())} sites), CCMP1948 has: {dict(agg)}")
print(f"at SAG's AGA/CGN ({sum(argc.values())} sites), CCMP1948 has: {dict(argc)}")
p = fisher_exact([[agg["ATG"], agg["Arg"]], [argc["ATG/AGG"], argc["other Arg"]]], alternative="greater")[1]
print(f"one-sided Fisher, changed AGG sites (ATG vs Arg) against changed Arg sites: p {p:.1e}")

# 4. Reads (data/reads/, from src/codes_census/reads_codons.py)
def reads(path):
    out = []
    for r in rows(path):
        r["reads"], r["same"] = int(r["reads"]), int(r["same"])
        r["other"] = ast.literal_eval(r["other"])
        out.append(r)
    return out


dna = [r for r in reads("reads/paer_SRR15114507.tsv") if r["codon"] == "AGG" and r["cons"] == "M"]
st = [r for r in dna if r["strict"] == "1"]
lo = min(st, key=lambda r: r["same"] / r["reads"])
print(f"DNA: {len(st)} strict-Met AGG sites, {min(r['reads'] for r in st)}-{max(r['reads'] for r in st)} reads each, "
      f"lowest share {lo['same']} of {lo['reads']} ({100 * lo['same'] / lo['reads']:.2f}%)")
print(f"DNA: all AGG at Met-consensus columns ({len(dna)} sites): {sum(r['same'] for r in dna)} of "
      f"{sum(r['reads'] for r in dna)} reads AGG; ATG itemized {sum(r['other'].get('ATG', 0) for r in dna)}")
rna = [r for r in reads("reads/paer_SRR1296983.tsv") if r["codon"] == "AGG" and r["cons"] == "M"]
rs = [r for r in rna if r["strict"] == "1" and r["reads"]]
print(f"RNA: {sum(r['same'] for r in rs)} of {sum(r['reads'] for r in rs)} reads AGG at {len(rs)} covered strict sites; "
      f"{sum(r['same'] for r in rna)} of {sum(r['reads'] for r in rna)} at Met-consensus sites")
