"""Diagnostic (not part of the rule): what do a run's reads share with the 24 TANB77 genomes?
  .venv/bin/python -I kmer_census.py <run> <max_spots>
Patterns (exact 31-mers, both orientations):
  rRNA  every 31-mer of the annotated rRNA features (5S, SSU, LSU) of the 24 genomes;
  CDS   31-mers every 128 nt along every catalogue CDS of the 22 primary genomes (a read lying >= 158 nt inside a
        CDS holds at least one);
  other every 31-mer of the other annotated ncRNA features (tRNA, tmRNA, RNase P, ...).
Each read is classed by which pattern kinds it holds; counts per class are written, nothing else.
Output: diag/<run>.census
"""
import csv
import glob
import subprocess
import sys
from collections import Counter

import ahocorasick

R = "results/rna_cgg/recheck/"
MDIR = "data/proteomics/mgnify_pig/"
COMP = str.maketrans("ACGTN", "TGCAN")
K, STRIDE = 31, 128


def rc(s):
    return s.translate(COMP)[::-1]


def fasta(path):
    seqs, name, buf = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf).upper()
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf).upper()
    return seqs


def main():
    run, xmax = sys.argv[1], int(sys.argv[2])
    primary = {r["genome"]: r["primary"] == "1" for r in csv.DictReader(
        open("results/proteomics_pig/genomes.tsv"), delimiter="\t")}
    kinds = {}
    for f in sorted(glob.glob(MDIR + "*.gff")):
        g = f.rsplit("/", 1)[1][:-4]
        seqs = fasta(MDIR + f"{g}.fna")
        for line in open(f):
            if line.startswith("##FASTA"):
                break
            if line[0] == "#":
                continue
            x = line.rstrip("\n").split("\t")
            if len(x) < 9:
                continue
            a, b = int(x[3]), int(x[4])
            s = seqs[x[0]][a - 1:b]
            if x[2] == "CDS":
                if not primary[g]:
                    continue
                kind, step = "CDS", STRIDE
            else:
                prod = dict(y.split("=", 1) for y in x[8].split(";") if "=" in y).get("product", "")
                kind, step = ("rRNA" if "ribosomal RNA" in prod else "other"), 1
            for i in range(0, len(s) - K + 1, step):
                w = s[i:i + K]
                for q in (w, rc(w)):
                    kinds.setdefault(q, set()).add(kind)
    A = ahocorasick.Automaton()
    for q, v in kinds.items():
        A.add_word(q, "+".join(sorted(v)))
    A.make_automaton()
    print("patterns", len(kinds), Counter("+".join(sorted(v)) for v in kinds.values()), flush=True)
    err = open(R + f"diag/{run}.census.err", "w")
    proc = subprocess.Popen(["fastq-dump", "--stdout", "--split-spot", "--skip-technical", "-X", str(xmax), run],
                            stdout=subprocess.PIPE, stderr=err, text=True, bufsize=1 << 20)
    cls, reads = Counter(), 0
    it = iter(proc.stdout)
    for h in it:
        seq = next(it)
        next(it)
        next(it)
        reads += 1
        hit = set()
        for _, v in A.iter(seq.rstrip()):
            hit.update(v.split("+"))
        cls["+".join(sorted(hit)) or "none"] += 1
    rcode = proc.wait()
    with open(R + f"diag/{run}.census", "w") as o:
        o.write(f"run\t{run}\nspots_requested\t{xmax}\nreads\t{reads}\nexit\t{rcode}\n")
        for k, v in cls.most_common():
            o.write(f"{k}\t{v}\n")
    print(run, reads, dict(cls), "exit", rcode, flush=True)


if __name__ == "__main__":
    main()
