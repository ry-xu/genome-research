"""Write pyrodigal (single mode) proteins of a genome: python -I write_faa.py GENOME.fna.gz OUT.faa"""
import gzip, sys
import pyrodigal
seqs = {}; name = None
for l in gzip.open(sys.argv[1], "rt"):
    if l.startswith(">"): name = l[1:].split()[0]; seqs[name] = []
    else: seqs[name].append(l.strip().upper())
seqs = {k: "".join(v) for k, v in seqs.items()}
gf = pyrodigal.GeneFinder(); gf.train(*[s.encode() for s in seqs.values() if len(s) > 20000][:200])
with open(sys.argv[2], "w") as o:
    for k, s in seqs.items():
        for g in gf.find_genes(s.encode()):
            o.write(f">{k}|{g.begin}|{g.end}|{g.strand}|{len(s)}\n{g.translate(include_stop=False)}\n")
