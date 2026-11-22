"""Synthetic recoding control: rewrite a fraction of one sense codon to a stop inside a real genome's genes, so the
genome 'reads that stop as sense'. Genes are pyrodigal (table 11) calls on the original genome; only internal codons
(not the first or last 10 codons) are rewritten. The result is cached as data/genomes/<label>.fna.gz, where
readthrough.fetch finds it without downloading.

usage: uv run --no-project --with pyrodigal,requests python src/synth.py <acc> <from_codon> <to_stop> <fraction> <label>
"""
import gzip
import random
import sys
from pathlib import Path

import pyrodigal

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent.parent / "stop-codon-census" / "src"))
import readthrough as rt  # noqa: E402

rt.GEN = HERE.parent / "data" / "genomes"

acc, frm, to, frac, label = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4]), sys.argv[5]
contigs = list(rt.read_fasta(rt.fetch(acc)))
finder = pyrodigal.GeneFinder()
finder.train(*[s.encode() for _, s in contigs], translation_table=11)
rng = random.Random(0)
n_all = n_done = 0
out = []
for name, seq in contigs:
    s = list(seq)
    for g in finder.find_genes(seq.encode()):
        if g.partial_begin or g.partial_end:
            continue
        lo, hi = g.begin - 1, g.end  # 0-based [lo, hi)
        L = hi - lo
        for k in range(10, L // 3 - 10):
            if g.strand == 1:
                p = lo + 3 * k
                cod = "".join(s[p:p + 3])
                if cod == frm:
                    n_all += 1
                    if rng.random() < frac:
                        s[p:p + 3] = list(to)
                        n_done += 1
            else:
                p = hi - 3 * (k + 1)  # codon k on the minus strand occupies [p, p+3) reverse-complemented
                cod = rt.revcomp("".join(s[p:p + 3]))
                if cod == frm:
                    n_all += 1
                    if rng.random() < frac:
                        s[p:p + 3] = list(rt.revcomp(to))
                        n_done += 1
    out.append((name, "".join(s)))
with gzip.open(rt.GEN / f"{label}.fna.gz", "wt") as f:
    for name, seq in out:
        f.write(f">{name}\n")
        for i in range(0, len(seq), 80):
            f.write(seq[i:i + 80] + "\n")
print(f"{label}: rewrote {n_done}/{n_all} internal {frm} -> {to} in {acc}")
