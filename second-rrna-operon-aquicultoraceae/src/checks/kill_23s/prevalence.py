"""Count Aquicultoraceae genomes with any 23S locus (>=300 bp blast locus) and with a divergent-labelled one (div_pid>=80)."""
import collections
tax = {l.split("\t")[0][3:]: (l.split("\t")[1], l.split("\t")[2]) for l in open("../aq/genomes.tsv") if not l.startswith("acc")}
any23, div = set(), collections.defaultdict(list)
for l in open("loci_blast.tsv"):
    f = l.rstrip("\n").split("\t")
    if f[0] != "f__Aquicultoraceae" or int(f[6]) < 300: continue
    any23.add(f[2])
    if f[7] == "div" and float(f[9]) >= 80: div[f[2]].append(int(f[6]))
aq = [g for g in tax if tax[g][0] == "f__Aquicultoraceae"]
print(f"Aquicultoraceae genomes {len(aq)}; with any 23S locus {len(any23)}; with divergent 23S {len(div)}")
for gen in sorted(set(tax[g][1] for g in aq)):
    gs = [g for g in aq if tax[g][1] == gen]
    print(f"  {gen}: genomes {len(gs)}, any23S {sum(g in any23 for g in gs)}, div {sum(g in div for g in gs)}")
other = [l for l in open("loci_blast.tsv") if not l.startswith(("family", "f__Aquicultoraceae")) and l.split("\t")[7] == "div"]
print("non-Aquicultoraceae loci labelled div:", len(other), [ (x.split('\t')[2], x.split('\t')[6], x.split('\t')[9]) for x in other])
