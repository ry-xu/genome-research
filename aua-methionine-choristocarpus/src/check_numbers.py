"""Recount the README's numbers from the tables in ../data/ (no downloads, no BLAST).
usage (from the folder root): uv run python -I src/check_numbers.py

- data/codons/sites_all.parquet: one row per target codon x reference protein (best tblastn HSP), written by
  src/kill/k1_audit.py, the kill agent's re-implementation of src/read/b_read.py. The rules are those of
  data/checks/RULES.md: skip reference residues 1-10; a position is Met (or Ile) when >= 2 references align there and
  >= 2/3 of them agree. The strict rule: >= 5 references, all agreeing.
- data/composition/: Altair's alignment-free codon usage (comp.out), the per-gene-group excess (groups.out) and the
  residues facing each codon in the relatives (cons.out).
- data/rna/k12_rnaseq.out, data/screen/d_test.out: printed as written.
"""
import re
from collections import Counter, defaultdict
from math import comb

import polars as pl

D = "data/"
NAMES = {"cten_WGS": "Choristocarpus (Genoscope assembly)", "POS_NC_085326": "Choristocarpus (GenBank NC_085326)",
         "dmes_WGS": "Discosporangium", "NEG_Phaeostrophion": "Phaeostrophion", "NEG_Protohalopteris": "Protohalopteris",
         "NEG_Sargassum": "Sargassum", "NEG_Saccharina": "Saccharina", "NEG_Petalonia": "Petalonia"}

s = pl.read_parquet(D + "codons/sites_all.parquet").filter(pl.col("qpos") > 10)
calls = {}  # target -> list of (codon, residue, n_refs, n_agree, gene)
for (t, contig, strand, pos), g in s.group_by(["target", "contig", "strand", "pos"]):
    res = Counter(g["res"].to_list())
    top, c = res.most_common(1)[0]
    n = g.height
    if n >= 2 and c / n >= 2 / 3:
        gene = Counter(q.split("_", 1)[1].lower() for q in g["qgene"].to_list()).most_common(1)[0][0]
        calls.setdefault(t, []).append((g["codon"][0], top, n, c, gene))

print("== AUN codons at conserved Met / Ile positions (8 brown-algal references; rule of RULES.md)")
for t in NAMES:
    C = Counter((cod, r) for cod, r, *_ in calls[t] if r in "MI")
    strict = Counter((cod, r) for cod, r, n, c, _ in calls[t] if r in "MI" and n >= 5 and c == n)
    print(f"{NAMES[t]:38s}", "  ".join(f"{x}: M {C[(x, 'M')]} / I {C[(x, 'I')]}" for x in ("ATA", "ATG", "ATT", "ATC")),
          f"| strict ATA: M {strict[('ATA', 'M')]} / I {strict[('ATA', 'I')]}")

print("\n== Choristocarpus (Genoscope)")
cm = [x for x in calls["cten_WGS"] if x[1] == "M"]
ci = [x for x in calls["cten_WGS"] if x[1] == "I"]
ata_m = sum(x[0] == "ATA" for x in cm); ata_i = sum(x[0] == "ATA" for x in ci)
print(f"AUA: {ata_m + ata_i} at Met or Ile positions, {ata_m} at Met")
genes = defaultdict(Counter)
for cod, r, n, c, gene in cm:
    genes[gene][cod] += 1
with_m = [g for g, v in genes.items() if v["ATA"] + v["ATG"] > 0]
print(f"genes with a conserved Met coded AUA or AUG: {len(with_m)}; with AUA at a Met position: "
      f"{sum(genes[g]['ATA'] > 0 for g in with_m)}")
for t in ("cten_WGS", "dmes_WGS", "NEG_Phaeostrophion", "NEG_Saccharina"):
    m = [x[0] for x in calls[t] if x[1] == "M"]
    print(f"Ile-type codons (AUA, AUU, AUC) at conserved Met positions, {NAMES[t]}: "
          f"{sum(c in ('ATA', 'ATT', 'ATC') for c in m)} of {len(m)} = {sum(c in ('ATA', 'ATT', 'ATC') for c in m) / len(m):.1%}")
# substitution null of the kill check: AUA's share among Ile codons at Ile positions vs at Met positions
ii = Counter(x[0] for x in ci); mm = Counter(x[0] for x in cm)
p0 = ii["ATA"] / (ii["ATA"] + ii["ATT"] + ii["ATC"]); nm = mm["ATA"] + mm["ATT"] + mm["ATC"]
pv = sum(comb(nm, k) * p0 ** k * (1 - p0) ** (nm - k) for k in range(mm["ATA"], nm + 1))
print(f"AUA share of Ile codons at Ile positions {p0:.3f}; at Met positions {mm['ATA']} of {nm}; binomial P {pv:.1e}")

print("\n== Codon usage without alignment (data/composition/comp.out: annotated CDS, 212 stramenopile mito genomes)")
comp = pl.read_csv(D + "composition/comp.out", separator="\t")
ch = comp.filter(pl.col("name") == "Choristocarpus tenellus")
br = comp.filter((pl.col("group") == "Phaeophyceae") & (pl.col("name") != "Choristocarpus tenellus"))
print(f"genomes: {comp.height}; lowest AUG per 1,000 codons: {comp.sort('ATG/1k')['name'][0]} {comp['ATG/1k'].min()}")
print(f"other brown algae: {br.height}, AUG per 1,000 {br['ATG/1k'].min()} to {br['ATG/1k'].max()}")
print(f"Choristocarpus AUG + AUA per 1,000 (Met if AUA = Met): {ch['ATG/1k'][0] + ch['ATA/1k'][0]:.1f}")

print("\n== Where the extra methionine sits (data/composition/groups.out)")
exc = {}
for line in open(D + "composition/groups.out"):
    m = re.match(r"(\S+).*ATG\s+(\d+) ATA\s+(\d+).*relatives Met\s+([\d.]+)", line)
    if m:
        exc[m.group(1)] = int(m.group(2)) + int(m.group(3)) - float(m.group(4))
tot = sum(exc.values())
for k, v in exc.items():
    print(f"{k}: {v:.1f} Met above the relatives' median ({v / tot:.0%} of the excess)")
for line in open(D + "composition/cons.out"):
    if line.startswith("POS_NC_085326") and " ATA " in line:
        print(line.rstrip())

print("\n== RNA (data/rna/k12_rnaseq.out)")
for line in open(D + "rna/k12_rnaseq.out"):
    if line.startswith(("reads total", "ATA@M")):
        print(line.rstrip())
print("\n== Screen (data/screen/d_test.out)")
for line in open(D + "screen/d_test.out"):
    if line.startswith(("population", "recall")):
        print(line.rstrip())
