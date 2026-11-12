"""Are TTA ends inherited? Codon at the same end in the other clade genomes (post hoc, Basil 2026-11-11).

usage: shared_ends.py <redteam_dir> <outdir>

Ends come from the red team's close-relative test (<redteam_dir>/ortho/runon_events.tsv, R = sister): a gene of
genome X ends where its best sister ortholog ends at a canonical stop. Kinds:
  TTA end        class runon, codon at e*+1 = TTA
  canonical end  class same, X's own stop codon (TAA / TAG / TGA)
For each end and each other genome Y of the same role (clade with clade, sister with sister) where X's gene has an
RBH ortholog (<redteam_dir>/ortho/rbh.tsv), the two proteins (own Pyrodigal calls, TTA = Leu) are aligned locally
(BLOSUM62, open -11, extend -1); X's residue e* must be aligned to a residue j of Y. Y's codon after j is
  Y's stop codon if j is Y's last residue, else Y's codon j+1 (sense; TTA or other).
Ends where X's residue e* falls in a gap or outside the alignment are counted as unresolved.

Writes <outdir>/pairs.tsv (one row per end x Y) and prints the tables by X kind, Y codon and genus relation.
"""
import collections
import csv
import sys

from Bio import Align
from Bio.Align import substitution_matrices

STOPS = {"TAA", "TAG", "TGA"}


def main(rd, outdir):
    panel = {r[0]: r for r in csv.reader(open(f"{rd}/panel.tsv"), delimiter="\t")}
    role = {k: ("clade" if r[2] == "clade" else "sister") for k, r in panel.items()}
    genus = {k: r[3] for k, r in panel.items()}
    nt, prot = {}, {}
    for lab in panel:
        for r in csv.DictReader(open(f"{rd}/genes/{lab}.genes.tsv"), delimiter="\t"):
            nt[f"{lab}|{r['gid']}"] = r["nt"]
        name = None
        for line in open(f"{rd}/genes/{lab}.faa"):
            if line.startswith(">"):
                name = line[1:].split()[0]
                prot[name] = []
            else:
                prot[name].append(line.strip())
    prot = {k: "".join(v).rstrip("*") for k, v in prot.items()}
    rbh = collections.defaultdict(dict)  # gene -> {genome: partner}
    for r in csv.DictReader(open(f"{rd}/ortho/rbh.tsv"), delimiter="\t"):
        rbh[r["q"]][r["s"].split("|")[0]] = r["s"]
    al = Align.PairwiseAligner(mode="local", open_gap_score=-11, extend_gap_score=-1)
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    out = []
    for r in csv.DictReader(open(f"{rd}/ortho/runon_events.tsv"), delimiter="\t"):
        x = r["B"]
        if r["R"] != "sister":
            continue
        if r["class"] == "runon" and r["w_p1"] == "TTA":
            kind = "TTA"
        elif r["class"] == "same":
            kind = r["w_p1"]
        else:
            continue
        g1 = f"{x}|{r['gene']}"
        es = int(r["e_star"])  # 1-based residue of X aligned to the sister's last residue
        for y, g2 in rbh[g1].items():
            if y == x or role[y] != role[x]:
                continue
            p1, p2 = prot[g1], prot[g2]
            a = al.align(p1, p2)[0]
            j = None
            for (s1, e1), (s2, e2) in zip(*a.aligned):
                if s1 < es <= e1:
                    j = s2 + (es - 1 - s1) + 1  # 1-based residue of Y
                    break
            n2 = nt[g2]
            ql2 = len(n2) // 3 - 1  # residues; codon index ql2 is the stop
            if j is None:
                yc = "unresolved"
            elif j == ql2:
                yc = n2[3 * ql2:3 * ql2 + 3]
            elif j < ql2:
                yc = n2[3 * j:3 * j + 3]
            else:
                yc = "unresolved"
            rel = "same genus" if genus[x] == genus[y] else "other genus"
            out.append(dict(role=role[x], X=x, Y=y, relation=rel, kind=kind, gene=g1, partner=g2, e_star=es,
                            Y_codon=yc, Y_class=("stop" if yc in STOPS else yc if yc in ("TTA", "unresolved")
                                                 else "sense")))
    with open(f"{outdir}/pairs.tsv", "w") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0]), delimiter="\t")
        w.writeheader()
        w.writerows(out)
    for rl in ("clade", "sister"):
        for rel in ("same genus", "other genus", "any"):
            print(f"\n## {rl} X vs {rl} Y, {rel}")
            print("X end\tn\tY TTA\tY TAA\tY TAG\tY TGA\tY sense\tunresolved\tY same codon %")
            for kind in ("TTA", "TAA", "TAG", "TGA"):
                rows = [o for o in out if o["role"] == rl and o["kind"] == kind
                        and (rel == "any" or o["relation"] == rel)]
                c = collections.Counter(o["Y_codon"] if o["Y_codon"] in STOPS | {"TTA", "unresolved"} else "sense"
                                        for o in rows)
                res = len(rows) - c["unresolved"]
                same = 100 * c[kind] / res if res else float("nan")
                print(kind, len(rows), c["TTA"], c["TAA"], c["TAG"], c["TGA"], c["sense"], c["unresolved"],
                      f"{same:.1f}", sep="\t")
    print("\n## clade TTA ends: per X genome, Y TTA / resolved, other-genus Y only")
    per = collections.defaultdict(lambda: [0, 0])
    for o in out:
        if o["role"] == "clade" and o["kind"] == "TTA" and o["relation"] == "other genus" and o["Y_class"] != "unresolved":
            per[(o["X"], o["Y"])][1] += 1
            per[(o["X"], o["Y"])][0] += o["Y_codon"] == "TTA"
    for (x, y), (k, n) in sorted(per.items()):
        print(x, genus[x], y, genus[y], f"{k}/{n}", sep="\t")


if __name__ == "__main__":
    main(*sys.argv[1:])
