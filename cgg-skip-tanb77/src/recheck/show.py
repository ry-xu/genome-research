import json, sys
W="recheck_work"
COD=["CGG","CGA","CGC","AGG","AGA","CGT","GCG","CCG","ACG","TCG","CTC","TGG","TAA","TAG","TGA"]
for lab in sys.argv[1:]:
    s=json.load(open(f"{W}/out/{lab}.summary.json"))
    print(lab, "loci",s["n_loci"],"events",s["n_events"],"+1",s["n_events_p1"],"-1",s["n_events_m1"],">=2refs",s["n_events_ge2refs"])
    for c in COD:
        x=s["codons"][c]
        top=", ".join(f"{k}{v}" for k,v in list(x["cons"].items())[:4])
        print(f"  {c} n={x['n_loci']:6d} in/out={x['inframe']}/{x['outframe']:.0f} cons_n={x['cons_n']:5d} [{top}] fs={x['fs']} D={x['D']} junc_old={x['junction_old']} junc_new={x['junction_new']}")
