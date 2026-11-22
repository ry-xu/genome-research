"""What does a relative have where genome A has an in-frame TGA? Table-4 genes of A and B (pyrodigal), reciprocal best
DIAMOND hits, each pair aligned globally (Biopython, BLOSUM62), and the codon B has at the residue aligned to each of A's
internal TGA codons. TGA = Trp in A predicts B's codon is TGG/TGA (a Trp position). Assembly errors or frameshift junk predict
no consistent residue. The same tally at A's TGG sites is the control.
usage: uv run --no-project --with pyrodigal,biopython,requests python -I src/ortholog_tga.py <accA> <accB>
genomes from data/sym_genomes (fetched if missing); prints the tally and writes results/ortholog_tga/<A>_vs_<B>.tsv"""
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path

import pyrodigal
from Bio import Align
from Bio.Align import substitution_matrices

HERE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(HERE.parent / "stop-codon-census" / "src"))
import readthrough as rt  # noqa: E402

rt.GEN = HERE / "data" / "sym_genomes"
OUT = HERE / "results" / "ortholog_tga"


def genes(acc):
    contigs = [s.upper().encode() for _, s in rt.read_fasta(rt.fetch(acc))]
    f = pyrodigal.GeneFinder()
    f.train(*contigs, translation_table=4)
    out = {}
    for i, c in enumerate(contigs):
        for j, g in enumerate(f.find_genes(c)):
            nt = g.sequence()
            out[f"{i}_{j}"] = (nt, g.translate(translation_table=4).rstrip("*"))
    return out


def rbh(ga, gb):
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        for name, g in (("a", ga), ("b", gb)):
            (td / f"{name}.faa").write_text("".join(f">{k}\n{p}\n" for k, (_, p) in g.items()))
            subprocess.run(["diamond", "makedb", "--in", str(td / f"{name}.faa"), "-d", str(td / name), "--quiet"], check=True)
        best = {}
        for q, d in (("a", "b"), ("b", "a")):
            subprocess.run(["diamond", "blastp", "-q", str(td / f"{q}.faa"), "-d", str(td / d), "-o", str(td / f"{q}{d}.m8"),
                            "--quiet", "-k", "1", "-e", "1e-5", "--sensitive", "--outfmt", "6", "qseqid", "sseqid"], check=True)
            best[q] = dict(l.split()[:2] for l in open(td / f"{q}{d}.m8"))
    return [(a, b) for a, b in best["a"].items() if best["b"].get(b) == a]


def main(acc_a, acc_b):
    ga, gb = genes(acc_a), genes(acc_b)
    pairs = rbh(ga, gb)
    al = Align.PairwiseAligner(mode="global", substitution_matrix=substitution_matrices.load("BLOSUM62"),
                               open_gap_score=-11, extend_gap_score=-1)
    tally = {"TGA": Counter(), "TGG": Counter()}
    rows = []
    for a, b in pairs:
        nta, pa = ga[a]
        ntb, pb = gb[b]
        aln = al.align(pa, pb)[0]
        amap = {}
        for (a0, a1), (b0, b1) in zip(*aln.aligned):
            for k in range(a1 - a0):
                amap[a0 + k] = b0 + k
        for i in range(1, len(pa)):  # skip the start codon
            cod = nta[3 * i:3 * i + 3]
            if cod not in tally:
                continue
            j = amap.get(i)
            other = ntb[3 * j:3 * j + 3] if j is not None else "gap"
            tally[cod][other] += 1
            rows.append(f"{a}\t{i}\t{cod}\t{b}\t{j if j is not None else ''}\t{other}")
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"{acc_a}_vs_{acc_b}.tsv").write_text("geneA\tcodon_idx\tcodonA\tgeneB\tcodon_idxB\tcodonB\n" + "\n".join(rows) + "\n")
    print(f"{acc_a} vs {acc_b}: {len(pairs)} reciprocal best pairs")
    for cod, t in tally.items():
        n = sum(t.values())
        trp = t["TGG"] + t["TGA"]
        print(f"  A has {cod} ({n} sites): B has TGG {t['TGG']}, TGA {t['TGA']} (Trp positions {trp / max(1, n):.0%});"
              f" top others {t.most_common(6)}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
