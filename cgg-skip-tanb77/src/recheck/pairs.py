import json, sys
W="recheck_work"
for lab in sys.argv[1:]:
    s=json.load(open(f"{W}/out/{lab}.summary.json"))
    for c in ["CGG","CGA","CGC","AGG"]:
        p=s["codons"][c]["pairs"]; n=sum(p.values())
        print(lab,c,"pairs",n, "R%.2f"%(p.get("R",0)/max(n,1)), list(p.items())[:6])
