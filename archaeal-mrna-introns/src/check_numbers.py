"""Recount the README's headline numbers directly from the tables in data/.

Usage (from the folder root):
  uv run python -I src/check_numbers.py
"""
import csv
from collections import Counter
from pathlib import Path

D = Path(__file__).resolve().parent.parent / "data"


def tsv(p):
    with open(p) as f:
        return list(csv.DictReader(f, delimiter="\t"))


# 1. 3,666-genome archaeal DNA census
genomes = tsv(D / "census" / "genomes_main.tsv")
assert len(genomes) == 3666, len(genomes)
ASGARD = {"Asgardarchaeota"}
TACK = {"Thermoproteota", "Korarchaeota"}
DPANN = {
    "Nanobdellota", "Micrarchaeota", "Aenigmatarchaeota", "Altiarchaeota",
    "Undinarchaeota", "Huberarchaeota", "Nanohalarchaeota", "Iainarchaeota",
}


def supergroup(phylum):
    if phylum in ASGARD:
        return "Asgard"
    if phylum in TACK:
        return "TACK"
    if phylum in DPANN:
        return "DPANN"
    return "Euryarchaeota side"


grp_counts = Counter()
grp_placed = Counter()
grp_bhb = Counter()
grp_exp = Counter()
grp_with1 = Counter()
for r in genomes:
    g = supergroup(r["phylum"])
    grp_counts[g] += 1
    grp_placed[g] += int(r["n_placed"])
    obs = int(r["obs"])
    grp_bhb[g] += obs
    grp_exp[g] += float(r["exp_raw"])
    if obs >= 1:
        grp_with1[g] += 1

print(f"Census genomes: {len(genomes)}")
for g in ("Asgard", "TACK", "DPANN", "Euryarchaeota side"):
    ratio = grp_bhb[g] / grp_exp[g]
    pct = 100.0 * grp_with1[g] / grp_counts[g]
    print(f"  {g:18s}: {grp_counts[g]:4d} genomes, {grp_bhb[g]:3d} / {grp_placed[g]:3d} placed BHB, "
          f"exp {grp_exp[g]:.2f}, obs/exp {ratio:.2f}x, >=1 hit {pct:.1f}%")
assert grp_counts["Asgard"] == 400 and grp_bhb["Asgard"] == 373 and grp_placed["Asgard"] == 641
assert grp_counts["TACK"] == 896 and grp_bhb["TACK"] == 517
assert grp_counts["DPANN"] == 1068 and grp_bhb["DPANN"] == 74
assert grp_counts["Euryarchaeota side"] == 1302 and grp_bhb["Euryarchaeota side"] == 32

# 2. Asgard B35 L02 RNA & DNA molecules
rna_rows = [r for r in tsv(D / "asgard_b35" / "rna_mol_templates.tsv") if r["site"] == "L02" and r["place"] == "0"]
rna_cls = Counter(r["class"] for r in rna_rows)
dna_rows = tsv(D / "asgard_b35" / "fig_l02_dna_templates.tsv")
dna_cls = Counter(r["class"] for r in dna_rows)
anyread = next(
    int(r.split("\t")[4])
    for r in open(D / "asgard_b35" / "rna_mol_gaps.tsv")
    if r.split("\t")[:3] == ["CP104013.1", "4487605", "31"]
)
print(f"\nAsgard B35 L02: RNA {rna_cls['joined']} lacking / {rna_cls['through']} through "
      f"(total {sum(rna_cls.values())}, {anyread} exact 31 bases); "
      f"DNA {dna_cls['joined']} lacking / {dna_cls['through']} through")
assert (rna_cls["joined"], rna_cls["through"], anyread) == (385, 26, 383)
assert (dna_cls["joined"], dna_cls["through"]) == (0, 46)

# 3. B35 PCR check (c1.out) and DNA hairpin check (c2.out)
c1_text = (D / "asgard_b35" / "c1.out").read_text().strip()
c2_lines = (D / "asgard_b35" / "c2.out").read_text().splitlines()
print(f"  B35 C1 (PCR): {c1_text.splitlines()[2]}")
print(f"  B35 C2 (DNA folds): {c2_lines[1]}")
assert "expected mixed 13.47; observed 0" in c1_text
assert "loci DNA dG <= -3.6: 3735; with share >= 0.10: 0" in c2_lines[1]

# 4. TACK Methanosuratincola verstraetei LCB70 c0 RNA, DNA, and 22-nt RNA hairpin background
f1_events = tsv(D / "tack_c0" / "f1.events.tsv")
c0_row = next(r for r in f1_events if r["study"] == "S2" and r["id"] == "GCA_049570475.1:c0")
metag_rows = tsv(D / "tack_c0" / "metag_words.tsv")
c0_dna = next(r for r in metag_rows if r["event"] == "GCA_049570475.1:c0")
hp_loci = [
    r for r in tsv(D / "tack_c0" / "hairpin_s2_W22_loci.tsv")
    if r["event"] == "" and float(r["dG"]) <= -5.6
]
hp_rem = sum(int(r["removal"]) for r in hp_loci)
hp_span = sum(int(r["spanning"]) for r in hp_loci)
print(f"\nTACK LCB70 c0: RNA {c0_row['joined']} lacking / {c0_row['through']} through "
      f"(total {int(c0_row['joined']) + int(c0_row['through'])}); "
      f"DNA {c0_dna['joined']} lacking / {c0_dna['through']} through")
print(f"  LCB70 22-nt RNA folds <= -5.6 kcal/mol: {len(hp_loci)} loci, {hp_rem} removal reads / {hp_span} spanning")
assert (int(c0_row["joined"]), int(c0_row["through"])) == (45, 6)
assert (int(c0_dna["joined"]), int(c0_dna["through"])) == (0, 58)
assert len(hp_loci) == 84 and hp_rem == 0

# 5. TACK c0 25-genome Methanomethylicaceae BHB fold check
hom_lines = (D / "tack_c0" / "c0_homologs2.out").read_text().splitlines()
print("\nTACK c0 25-genome fold summary:")
for l in hom_lines:
    if any(k in l for k in ("central", "3' bulge", "25/25", "copies")):
        print(" ", l)

print("\nAll assertions passed.")
