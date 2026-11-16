import json, sys
W="recheck_work"
d=json.load(open(f"{W}/out/fs_detail_summary.json"))
for lab in sys.argv[1:]:
    s=json.load(open(f"{W}/out/{lab}.summary.json"))["codons"]
    g={c[2]: s["GG"+c[2]]["n_loci"] for c in ["GGA","GGT","GGC","GGG"]}; n=sum(g.values())
    nx=d.get(lab,{}).get("next_nt_p1",{}); m=sum(nx.values()) or 1
    print(lab, "Gly codon 3rd base usage:", {k: round(v/n,2) for k,v in g.items()}, "| nt after CGG at +1 sites:", {k: round(nx.get(k,0)/m,2) for k in "ATCG"}, "n=",m)
