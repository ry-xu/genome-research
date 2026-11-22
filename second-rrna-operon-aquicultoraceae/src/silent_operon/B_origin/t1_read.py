"""Read a T1 tree: root on Pseudomonadota+Bacteroidota leaves (Gracilicutes), report divergent-clade monophyly,
its support, its sister group, and the nested position; draw an ASCII tree.
Usage: python -I src/t1_read.py out/t1_23S.contree
"""
import sys
from Bio import Phylo

t = Phylo.read(sys.argv[1], "newick")
leaves = t.get_terminals()
name = lambda c: c.name or ""
og = [l for l in leaves if any(x in name(l) for x in ("Pseudomonadota", "Bacteroidota", "ec_", "vv_", "tt_"))]
try:
    t.root_with_outgroup(og)
except Exception as e:
    print("rooting on full outgroup failed:", e)
    t.root_with_outgroup([l for l in leaves if name(l).startswith("ec_")])
div = [l for l in leaves if name(l).endswith("__div")]
print("divergent leaves:", len(div))
ca = t.common_ancestor(div)
inside = ca.get_terminals()
print("divergent copies monophyletic:", len(inside) == len(div), "| clade support:", ca.confidence)
if len(inside) != len(div):
    print("  intruders:", [name(l) for l in inside if l not in div])
path = t.get_path(ca)
parent = path[-2] if len(path) >= 2 else t.root
sis = [c for c in parent.clades if c is not ca]
for s in sis:
    print(f"sister clade (parent support {parent.confidence}): {len(s.get_terminals())} leaves")
    for l in s.get_terminals()[:40]:
        print("   ", name(l))
# walk up: list successive nodes with support and composition
print("\nsuccessive ancestors of the divergent clade:")
for node in reversed(path[:-1]):
    ts = node.get_terminals()
    comp = {}
    for l in ts:
        k = name(l).split("__")[1] if "__" in name(l) else name(l)
        k = "AQ-norm" if (k.startswith("AQ-") and name(l).endswith("__norm")) else "AQ-div" if name(l).endswith("__div") else \
            ("order" if name(l).endswith("__order") else k.split("-")[0] + "-" + (k.split("-")[1] if "-" in k else ""))
        comp[k] = comp.get(k, 0) + 1
    print(f"  support {node.confidence}: n={len(ts)}  " + ", ".join(f"{k}:{v}" for k, v in sorted(comp.items(), key=lambda x: -x[1])[:8]))
    if len(ts) > 70:
        break
# divergent subtree
print("\ndivergent subtree:")
Phylo.draw_ascii(ca, column_width=90)
# tip-to-root distances (branch length from divergent MRCA vs from normal-family MRCA)
nf = [l for l in leaves if "__AQ-" in name(l) and name(l).endswith("__norm")]
if nf:
    cn = t.common_ancestor(nf)
    print("family normal copies monophyletic:", len(cn.get_terminals()) == len(nf), "support", cn.confidence)
    print("mean root-to-tip from divergent MRCA: %.3f; from family-normal MRCA: %.3f; stem branch of divergent clade: %.3f" % (
        sum(t.distance(ca, l) for l in div) / len(div), sum(t.distance(cn, l) for l in nf) / len(nf), ca.branch_length or 0))
