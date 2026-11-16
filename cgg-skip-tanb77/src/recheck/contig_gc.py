import sys
sys.path.insert(0, "src/recheck")
from mp import read_fasta
W = "recheck_work"
lab = sys.argv[1]; g = read_fasta(f"{W}/genomes/{lab}.fna")
tot = sum(len(s) for s in g.values()); gc = sum(s.count("G") + s.count("C") for s in g.values()) / tot
print(lab, "genome GC %.3f, %d contigs" % (gc, len(g)))
for c in sys.argv[2:]:
    s = g[c]; print(" ", c, len(s), "GC %.3f" % ((s.count("G") + s.count("C")) / len(s)), "CGG count (both strands)", s.count("CGG") + s.count("CCG"))
