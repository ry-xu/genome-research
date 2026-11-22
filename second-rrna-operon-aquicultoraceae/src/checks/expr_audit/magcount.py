"""Count RNA spots on the CG MAG (GCA_009993245.1) from the mate-split remap.
usage: magcount.py MAG.paf GENES.faa [MINID] [MINLEN]"""
import sys, collections
paf, faa = sys.argv[1], sys.argv[2]
minid = float(sys.argv[3]) if len(sys.argv) > 3 else 0.98
minlen = int(sys.argv[4]) if len(sys.argv) > 4 else 60
R = {"div16S(069:23800-25600)": ("JAACTV010000069.1", 23800, 25600), "divITS(069:23507-23799)": ("JAACTV010000069.1", 23507, 23799),
     "div23S(069:20469-23506)": ("JAACTV010000069.1", 20469, 23506),
     "norm16S(080:10281-11821)": ("JAACTV010000080.1", 10281, 11821), "normITS(080:11822-12237)": ("JAACTV010000080.1", 11822, 12237),
     "norm23S_5p(080:12238-14589)": ("JAACTV010000080.1", 12238, 14589), "norm23S_insert(080:14590-15574)": ("JAACTV010000080.1", 14590, 15574),
     "norm23S_3p(080:15575-16204)": ("JAACTV010000080.1", 15575, 16204)}
genes = collections.defaultdict(list)
for l in open(faa):
    if l.startswith(">"):
        _, c, a, b, s = l[1:].strip().split("|"); genes[c].append((int(a), int(b), s))
def region(c, a, b):
    for k, (rc, ra, rb) in R.items():
        if c == rc and a >= ra - 20 and b <= rb + 20: return k
        if c == rc and a < rb and b > ra: return k + "[partial]"
    return None
reg = collections.defaultdict(set); gene_spots = collections.defaultdict(set); idhist = collections.Counter()
other = set(); ident_all = []
for l in open(paf):
    f = l.split("\t"); q, c, ts, te = f[0].rsplit("/", 1)[0], f[5], int(f[7]), int(f[8]); m, al = int(f[9]), int(f[10])
    if al < minlen: continue
    ident = m / al
    r = region(c, ts, te)
    if r is None:
        idhist[round(ident, 2) if ident >= 0.9 else "<0.90"] += 1
    if ident < minid: continue
    if r: reg[r].add(q); continue
    hit = [g for g in genes[c] if ts < g[1] and te > g[0]]
    if hit:
        for g in hit: gene_spots[(c, g[0], g[1], g[2])].add(q)
    else: other.add(q)
print(f"# id>={minid} alen>={minlen}")
for k in sorted(reg): print(f"region\t{k}\t{len(reg[k])} spots")
tot = sum(len(v) for v in gene_spots.values())
print(f"protein-gene spots genome-wide: {tot} over {len(gene_spots)} genes; intergenic/other non-rRNA spots: {len(other)}")
for c, lab in (("JAACTV010000069.1", "div-operon contig"), ("JAACTV010000080.1", "norm-operon contig")):
    print(f"## {lab} genes with spots:")
    for k in sorted(x for x in gene_spots if x[0] == c): print("  ", k, len(gene_spots[k]))
print("non-rRNA alignment identity histogram (all alignments >=minlen, any id):", dict(sorted(idhist.items(), key=lambda x: str(x[0]))))
top = sorted(gene_spots.items(), key=lambda x: -len(x[1]))[:10]
print("top genes:", [(k[0], k[1], k[2], len(v)) for k, v in top])
