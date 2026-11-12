"""Blinded genome panel for an outside re-check by Clover (post hoc; not one of the frozen checks).

Panel: the 8-tip Polyangiaceae clade (the seven TTA-flagged genomes plus CAYNEK01) and the 8 Polyangiaceae tips
outside it that are nearest to the clade's MRCA by branch length on the GTDB R232 bac120 tree. The 16 genomes
are shuffled with a fixed seed and renamed p01..p16. Contigs are renamed pXX_c1.. so headers do not carry the
accession. Output goes to blind_panel/ (FASTA only) and the key to results/blind_panel_key.tsv (not for Clover).

usage: uv run --no-project --with biopython,requests,pyrodigal,numpy python -I src/blind_panel.py
"""
import csv
import gzip
import random
import sys
from io import StringIO
from pathlib import Path

from Bio import Phylo

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
import fetch  # noqa: E402

TREE = ROOT.parent / "tag-bacteria" / "data" / "bac120_r232.tree.gz"
N_CTRL = 8
SEED = 20261111


def family_subtree():
    s = gzip.open(TREE, "rt").read()
    i = s.find("f__Polyangiaceae")
    j = i - 1
    while s[j] != ")":
        j -= 1
    depth, k = 0, j
    while True:
        if s[k] == ")":
            depth += 1
        elif s[k] == "(":
            depth -= 1
            if depth == 0:
                break
        k -= 1
    return Phylo.read(StringIO(s[k:j + 1] + ";"), "newick")


def rf2(r):
    for x in r["rf"].split():
        if x.startswith("RF2:"):
            return x[4:]
    return "-"


def main():
    rows = {r["acc"]: r for r in csv.DictReader(open(ROOT / "results/gtdb/clade_TTA_Streptosporangiaceae.tsv"),
                                                delimiter="\t") if r["family"] == "f__Polyangiaceae"}
    t = family_subtree()
    tips = t.get_terminals()
    core = [tp for tp in tips if tp.name in rows and (rf2(rows[tp.name]) == "SPY" or "TTA" in rows[tp.name]["armA_flags"])
            and rows[tp.name]["genus"] != "g__JAFGBI01"]
    m = t.common_ancestor(core)
    clade = m.get_terminals()
    assert len(clade) == 8, len(clade)
    inside = {tp.name for tp in clade}
    others = sorted((t.distance(m, tp), tp.name) for tp in tips if tp.name not in inside)
    ctrl = [n for _, n in others[:N_CTRL]]
    dist = {n: d for d, n in others}
    panel = [(n, "clade") for n in sorted(inside)] + [(n, "control") for n in ctrl]
    random.Random(SEED).shuffle(panel)
    out = ROOT / "blind_panel"
    out.mkdir(exist_ok=True)
    key = []
    for i, (acc, role) in enumerate(panel, 1):
        pid = f"p{i:02d}"
        src = fetch.fetch(acc[3:])
        n = 0
        with gzip.open(src, "rt") as f, open(out / f"{pid}.fna", "w") as g:
            for line in f:
                if line.startswith(">"):
                    n += 1
                    g.write(f">{pid}_c{n}\n")
                else:
                    g.write(line)
        r = rows.get(acc, {})
        d = dist.get(acc) if role == "control" else None
        key.append({"id": pid, "acc": acc, "role": role, "genus": r.get("genus", "?"), "rf2": rf2(r) if r else "?",
                    "armA_TTA": int("TTA" in r.get("armA_flags", "")), "dist_to_mrca": f"{d:.4f}" if d is not None else "",
                    "contigs": n})
        print(pid, acc, role, r.get("genus", "?"), n, "contigs", flush=True)
    with open(ROOT / "results" / "blind_panel_key.tsv", "w") as f:
        w = csv.DictWriter(f, fieldnames=list(key[0]), delimiter="\t")
        w.writeheader()
        w.writerows(key)


if __name__ == "__main__":
    main()
