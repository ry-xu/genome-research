"""Do proteins near the divergent operon look like the rest of the DASXVE01 genome?
For each DASXVE protein: best hit outside DASXVE01 (family of hit genome, identity). Compare genes within 30 kb of the
divergent operon (DASXVE010000023.1:1-37000) vs genes on the normal-operon contig and vs the whole genome. Also GC and depth."""
import collections, statistics as st
tax = {l.split("\t")[0][3:]: l.split("\t")[1] for l in open("../aq/genomes.tsv") if not l.startswith("acc")}
best = {}
for l in open("dx_vs_all.tsv"):
    q, s, pid, ln, ql, sl, ev, bs = l.split("\t")
    g = s.split("|")[0]
    if g == "GCA_036504445.1": continue
    if q not in best or float(bs) > best[q][2]: best[q] = (tax[g], float(pid), float(bs), g)
genes = [l[1:].strip() for l in open("faa/GCA_036504445.1.faa") if l.startswith(">")]
def summ(sel, label):
    n = len(sel); hit = [best[q] for q in sel if q in best]
    fam = collections.Counter(h[0] for h in hit)
    aq = [h[1] for h in hit if h[0] == "f__Aquicultoraceae"]
    print(f"{label}: genes {n}, with hit {len(hit)} ({len(hit)/n:.2f}); best-hit family Aquicultoraceae {fam['f__Aquicultoraceae']/max(1,len(hit)):.2f};"
          f" median %id to Aquicultoraceae best hits {st.median(aq) if aq else float('nan'):.1f}; top families {fam.most_common(3)}")
def on(c, a, b): return [q for q in genes if q.split("|")[1] == c and int(q.split("|")[2]) >= a and int(q.split("|")[3]) <= b]
summ(genes, "whole genome")
summ(on("DASXVE010000023.1", 1, 37000), "div-operon contig 23, 1-37 kb")
summ(on("DASXVE010000023.1", 1, 79122), "div-operon contig 23, whole")
summ(on("DASXVE010000010.1", 1, 70000), "normal-operon contig 10, whole")
print("genes flanking the divergent operon (contig 23, first 20 kb):")
for q in on("DASXVE010000023.1", 1, 20000):
    print(" ", q, best.get(q))
# GC
seq = {}; name = None
for l in open("../kd/dxG.ref.fa"):
    if l.startswith(">"): name = l[1:].split()[0]; seq[name] = []
    else: seq[name].append(l.strip().upper())
seq = {k: "".join(v) for k, v in seq.items()}
gc = lambda s: (s.count("G") + s.count("C")) / max(1, len(s) - s.count("N"))
allseq = "".join(seq.values())
print(f"GC whole genome {gc(allseq):.3f}; contig23 {gc(seq['DASXVE010000023.1']):.3f}; contig23 without 1-8000 {gc(seq['DASXVE010000023.1'][8000:]):.3f};"
      f" div 23S 4036-7241 {gc(seq['DASXVE010000023.1'][4035:7241]):.3f}; normal 23S contig10 42611-45607 {gc(seq['DASXVE010000010.1'][42610:45607]):.3f}")
