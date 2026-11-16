"""Call genes in standard-code reference genomes with pyrodigal (table 11) and write proteins.
Residues translated from CGG codons are written as X so the reference never 'votes' on CGG itself.
Also decompress the 9 focal genomes to plain FASTA for miniprot."""
import gzip, sys, os
import pyrodigal
from Bio import SeqIO

W = "recheck_work"
G = "data/atlas/genomes"
REFS = {  # acc: short name (all NCBI table 11 isolates, GTDB RS_)
    "GCF_000009045.1": "Bsub", "GCF_000009205.2": "Cdif", "GCF_000981035.1": "Chon",
    "GCF_023347275.1": "Fbut", "GCF_034334985.1": "Aint", "GCF_002222595.2": "Bhan",
    "GCF_018118585.1": "Alac", "GCF_000146185.1": "Leli",
}
FOCAL = {"G1": "GCA_946496695.1", "G2": "GCF_012317385.1", "G3": "GCA_900553755.1", "G4": "GCA_963931445.1",
         "G5": "GCF_005473905.2", "G6": "GCA_964631135.1", "G7": "GCA_009785865.1", "G8": "GCF_900187165.1",
         "G9": "GCA_949283575.1"}
BASES = "TCAG"; AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
STD = {a + b + c: x for (a, b, c), x in zip([(a, b, c) for a in BASES for b in BASES for c in BASES], AA)}

def path(acc): return f"{G}/{acc[:8]}/{acc}.fna.gz"

for lab, acc in FOCAL.items():
    out = f"{W}/genomes/{lab}.fna"
    if not os.path.exists(out):
        with gzip.open(path(acc), "rt") as f, open(out, "w") as o:
            for r in SeqIO.parse(f, "fasta"):
                o.write(f">{r.id}\n{str(r.seq).upper()}\n")

allf = open(f"{W}/refs/refs_all.faa", "w")
for acc, nm in REFS.items():
    with gzip.open(path(acc), "rt") as f:
        recs = [(r.id, str(r.seq).upper()) for r in SeqIO.parse(f, "fasta")]
    finder = pyrodigal.GeneFinder(meta=False)
    finder.train(*[s.encode() for _, s in recs], translation_table=11)
    n = 0
    for cid, s in recs:
        for i, g in enumerate(finder.find_genes(s.encode())):
            nt = g.sequence()
            prot = []
            for k in range(0, len(nt) - 3, 3):  # skip final stop
                c = nt[k:k + 3]
                a = STD.get(c, "X")
                if c == "CGG" or a == "*": a = "X"
                prot.append(a)
            if k == 0 or not g.partial_begin:
                prot[0] = "M" if prot else ""
            p = "".join(prot)
            if len(p) >= 60:
                allf.write(f">{nm}|{cid}|{g.begin}|{g.end}|{g.strand}\n{p}\n"); n += 1
    print(nm, acc, n, file=sys.stderr)
allf.close()
