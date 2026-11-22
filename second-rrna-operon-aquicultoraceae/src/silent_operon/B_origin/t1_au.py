"""T1 topology tests: constrained ML trees + AU test (IQ-TREE 3) for where the divergent clade attaches.
H_order : divergent copies + normal copies of Aquicultorales (Aquicultoraceae + JAMCWL01) form one clade
H_class : divergent copies + all Aquicultoria normal copies form one clade (duplication inside the class)
H_phylum: divergent copies + all Actinomycetota copies form one clade
Usage: python -I src/t1_au.py 23S|16S
"""
import sys, subprocess, re
from Bio import SeqIO

k = sys.argv[1]
names = [r.id for r in SeqIO.parse(f"out/t1_{k}.trim.fa", "fasta")]
# taxonomy of order-level copies from loci table families -> class Aquicultoria for all 65 genomes
contam = {"23S_054"}  # DASXVE01 third 23S: groups with Clostridia (contaminant), left unconstrained
div = [n for n in names if n.endswith("__div") or n.startswith("16S_024__")]  # JAMCWL01 divergent 16S
aq_norm = [n for n in names if "__AQ-" in n and n.endswith("__norm") and n.split("__")[0] not in contam]
jam_norm = [n for n in names if "JAMCWL01" in n and n.endswith("__order") and not n.startswith("16S_024__")]
order_norm = [n for n in names if n.endswith("__order") and n not in jam_norm and not n.startswith("16S_024__")
              and not n.startswith("16S_046__")]
class_out = [n for n in names if "Actinomycetota-Aquicultoria" in n]
actino = [n for n in names if "Actinomycetota" in n or n.startswith("tb_")]
H = {"H_order": div + aq_norm + jam_norm,
     "H_class": div + aq_norm + jam_norm + order_norm + class_out}
for h, grp in H.items():
    rest = [n for n in names if n not in grp]
    with open(f"out/t1_{k}_{h}.constr", "w") as fh:
        fh.write("((" + ",".join(grp) + ")," + ",".join(rest) + ");\n")
    import os
    if not os.path.exists(f"out/t1_{k}_{h}.treefile"):
        subprocess.run(["iqtree3", "-s", f"out/t1_{k}.trim.fa", "-m", "GTR+F+R4", "-g", f"out/t1_{k}_{h}.constr", "-T", "1",
                        "--prefix", f"out/t1_{k}_{h}", "-quiet"], check=True)
with open(f"out/t1_{k}_cand.trees", "w") as fh:
    fh.write(open(f"out/t1_{k}.treefile").read())
    for h in H:
        fh.write(open(f"out/t1_{k}_{h}.treefile").read())
subprocess.run(["iqtree3", "-s", f"out/t1_{k}.trim.fa", "-m", "GTR+F+R4", "-z", f"out/t1_{k}_cand.trees", "-te",
                f"out/t1_{k}.treefile", "-n", "0", "-zb", "10000", "-au", "-T", "1", "--prefix", f"out/t1_{k}_au",
                "-redo", "-quiet"], check=True)
txt = open(f"out/t1_{k}_au.iqtree").read()
blk = txt[txt.index("USER TREES"):]
print(k, "trees: 1=ML, 2=H_order, 3=H_class")
print("\n".join(blk.splitlines()[:16]))
