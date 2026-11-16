"""List the CGG loci in frameshift genomes that are NOT followed by a frameshift (same frame), with residues."""
import sys, ast
W = "recheck_work"
for lab in sys.argv[1:]:
    fs = set()
    for l in list(open(f"{W}/out/{lab}.fs_cgg_sites.tsv"))[1:]:
        f = l.split("\t"); fs.add((f[0], f[1], int(f[2])))
    for l in list(open(f"{W}/out/{lab}.cgg_loci.tsv"))[1:]:
        f = l.rstrip("\n").split("\t")
        k = (f[0], f[1], int(f[2]))
        d = ast.literal_eval(f[8]) if f[8] else {}
        if k in fs or d.get("+1") or d.get("-1"): continue
        print(lab, *k, "cons=" + f[5], f[7], sep="\t")
