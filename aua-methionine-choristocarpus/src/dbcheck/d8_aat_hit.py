"""The weak tRNAscan -O 'Ile AAT' hit in Choristocarpus (NC_085326 18773-18841, score 20.2): compare with the adjacent
trnI(GAU) and with every tRNA in the six genomes; print its structure lines from the .ss file."""
B = "lab/work/organelle-reading/checks/chor_wgs/dbcheck/"
exec(open(B + "d6_helpers.py").read())
ss = parse_ss(B + "trna/NC_085326.tscanO.ss")
q = ss[15][2]
print("query NC_085326 trna15", ss[15][:2], "len", len(q), q)
txt = open(B + "trna/NC_085326.tscanO.ss").read()
blk = [b for b in txt.split("\n\n") if b.startswith("NC_085326.1.trna15 ")][0]
print(blk)
best = []
for g in G:
    s2 = parse_ss(B + f"trna/{g}.tscanO.ss"); rows = {int(p[1]): p for p in parse_tscan(B + f"trna/{g}.tscanO.txt")}
    for n, (b, e, s) in s2.items():
        if g == "NC_085326" and n == 15:
            continue
        aln = al.align(q, s)[0]
        m = sum(q[i] == s[j] for (a1, a2), (b1, b2) in zip(*aln.aligned) for i, j in zip(range(a1, a2), range(b1, b2)))
        best.append((m / max(len(q), len(s)), m, g, b, e, rows[n][4] + "-" + rows[n][5]))
best.sort(reverse=True)
print("top 8 matches (matches / longer length):")
for x in best[:8]:
    print(f"   {x[0]:.3f} ({x[1]})  {NAME[x[2]]} {x[3]}-{x[4]} {x[5]}")
