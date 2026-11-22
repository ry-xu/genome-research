"""pyrodigal (meta mode) proteins for each genome FASTA; headers ACC|CONTIG|start|end|strand"""
import sys, os, glob
from multiprocessing import Pool
import pyrodigal
def run(fn):
    acc = os.path.basename(fn)[:-4]
    out = f"faa/{acc}.faa"
    if os.path.exists(out): return
    gf = pyrodigal.GeneFinder(meta=True)
    seqs, name = {}, None
    for l in open(fn):
        if l.startswith(">"): name = l[1:].split()[0]; seqs[name] = []
        else: seqs[name].append(l.strip())
    with open(out + ".tmp", "w") as o:
        for c, s in seqs.items():
            s = "".join(s).upper()
            for g in gf.find_genes(s.encode()):
                o.write(f">{acc}|{c}|{g.begin}|{g.end}|{'+' if g.strand==1 else '-'}\n{g.translate().rstrip('*')}\n")
    os.rename(out + ".tmp", out)
if __name__ == "__main__":
    with Pool(3) as p: p.map(run, sorted(glob.glob("../aq/fna/*.fna")))
