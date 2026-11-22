"""Is the divergent operon at a conserved genomic position across genera/studies?
For each 23S locus (loci_blast.tsv, span>=1000) take the up-to-6 nearest protein genes on each side (within 15 kb).
For the DASXVE01 divergent locus and DASXVE01 normal locus (contig 10), count how many of their flank genes have their best
DIAMOND hit (in that target genome) among the flank genes of each target locus. Div-flank->div-locus vs div-flank->norm-locus
is the comparison; norm-flank rows are the control."""
import collections
loci = []
for l in open("loci_blast.tsv"):
    f = l.rstrip("\n").split("\t")
    if f[0] == "family" or int(f[6]) < 1000: continue
    loci.append((f[2], f[3], int(f[4]), int(f[5]), f[7], f[1]))
genes = collections.defaultdict(list)
for l in open("faa_all.faa"):
    if l.startswith(">"):
        g, c, a, b, s = l[1:].strip().split("|"); genes[(g, c)].append((int(a), int(b), l[1:].strip()))
def flanks(g, c, a, b, k=6, w=15000):
    gs = sorted(genes[(g, c)])
    up = [x[2] for x in gs if x[1] < a and x[1] > a - w][-k:]
    dn = [x[2] for x in gs if x[0] > b and x[0] < b + w][:k]
    return up, dn
best = {}
for l in open("dx_vs_all.tsv"):
    q, s, pid, ln, ql, sl, ev, bs = l.split("\t")
    tg = s.split("|")[0]
    if (q, tg) not in best or float(bs) > best[(q, tg)][1]: best[(q, tg)] = (s, float(bs), float(pid))
src = {"DASXVE_div": ("GCA_036504445.1", "DASXVE010000023.1", 4035, 7241), "DASXVE_norm": ("GCA_036504445.1", "DASXVE010000010.1", 42611, 45607)}
for nm, (g, c, a, b) in src.items():
    up, dn = flanks(g, c, a - 2000, b)   # a-2000: exclude 16S region; flank genes outside operon
    fl = up + dn
    print(f"# {nm}: {len(up)} upstream + {len(dn)} downstream flank genes")
    print("target_genome\tgenus\ttarget_contig\tlocus\tlabel\tflank_genes_in_target\tsource_flanks_hitting_them")
    for tg, tc, ta, tb, lab, gen in loci:
        if tg == g: continue
        tu, td = flanks(tg, tc, ta - 2000 if lab else ta, tb)
        tset = set(tu + td)
        n = sum(1 for q in fl if (q, tg) in best and best[(q, tg)][0] in tset)
        print(f"{tg}\t{gen}\t{tc}\t{ta}-{tb}\t{lab}\t{len(tset)}\t{n}/{len(fl)}")
