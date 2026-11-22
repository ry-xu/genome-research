"""Do the TGA-reading Throscidae Shikimatogenerans form one clade? Table-4 genes (pyrodigal) of each genome; orthologs are
reciprocal best DIAMOND hits to the first genome, kept when every genome has one; each family aligned with mafft, the
alignments concatenated, and a tree built with iqtree3 (-T 1, LG+F+G4, 1000 ultrafast bootstraps).
usage: uv run --no-project --with pyrodigal,requests python -I src/throscidae_tree.py <acc[:label]>...  (first = reference)
writes results/throscidae_tree/{concat.faa,concat.faa.treefile,...}"""
import subprocess
import sys
from pathlib import Path

import pyrodigal

HERE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(HERE.parent / "stop-codon-census" / "src"))
import readthrough as rt  # noqa: E402

rt.GEN = HERE / "data" / "sym_genomes"
OUT = HERE / "results" / "throscidae_tree"


def proteins(acc):
    contigs = [s.upper().encode() for _, s in rt.read_fasta(rt.fetch(acc))]
    f = pyrodigal.GeneFinder()
    f.train(*contigs, translation_table=4)
    return {f"{i}_{j}": g.translate(translation_table=4).rstrip("*")
            for i, c in enumerate(contigs) for j, g in enumerate(f.find_genes(c))}


def best_hits(q, d, out):
    subprocess.run(["diamond", "blastp", "-q", str(q), "-d", str(d), "-o", str(out), "--quiet", "-k", "1", "-e", "1e-5",
                    "--sensitive", "--outfmt", "6", "qseqid", "sseqid"], check=True)
    return dict(l.split()[:2] for l in open(out))


def main(args):
    OUT.mkdir(parents=True, exist_ok=True)
    genomes = [(a.partition(":")[0], a.partition(":")[2] or a.partition(":")[0]) for a in args]
    prot = {}
    for acc, lab in genomes:
        prot[lab] = proteins(acc)
        (OUT / f"{lab}.faa").write_text("".join(f">{k}\n{p}\n" for k, p in prot[lab].items()))
        subprocess.run(["diamond", "makedb", "--in", str(OUT / f"{lab}.faa"), "-d", str(OUT / lab), "--quiet"], check=True)
    ref = genomes[0][1]
    orth = {k: {ref: k} for k in prot[ref]}
    for _, lab in genomes[1:]:
        fwd = best_hits(OUT / f"{ref}.faa", OUT / lab, OUT / f"{ref}_{lab}.m8")
        rev = best_hits(OUT / f"{lab}.faa", OUT / ref, OUT / f"{lab}_{ref}.m8")
        for k in list(orth):
            b = fwd.get(k)
            if b and rev.get(b) == k:
                orth[k][lab] = b
            else:
                del orth[k]
    print(f"{len(orth)} ortholog families in all {len(genomes)} genomes", file=sys.stderr)
    concat = {lab: [] for _, lab in genomes}
    for k, members in orth.items():
        fam = OUT / f"fam_{k}.faa"
        fam.write_text("".join(f">{lab}\n{prot[lab][g]}\n" for lab, g in members.items()))
        aln = subprocess.run(["mafft", "--auto", "--quiet", str(fam)], check=True, capture_output=True, text=True).stdout
        seqs, name = {}, None
        for line in aln.splitlines():
            if line.startswith(">"):
                name = line[1:].strip()
                seqs[name] = []
            else:
                seqs[name].append(line.strip())
        for lab in concat:
            concat[lab].append("".join(seqs[lab]))
        fam.unlink()
    (OUT / "concat.faa").write_text("".join(f">{lab}\n{''.join(s)}\n" for lab, s in concat.items()))
    print(f"{len(''.join(concat[ref]))} alignment columns", file=sys.stderr)
    subprocess.run(["iqtree3", "-s", str(OUT / "concat.faa"), "-m", "LG+F+G4", "-B", "1000", "-T", "1", "--redo",
                    "-o", genomes[-1][1], "--quiet"], check=True)
    print((OUT / "concat.faa.treefile").read_text())


if __name__ == "__main__":
    main(sys.argv[1:])
