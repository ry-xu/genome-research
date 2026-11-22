"""Library strandedness and genome-wide background density on the MAG. usage: strand_density.py MAG.paf GENES.faa"""
import sys, collections
paf, faa = sys.argv[1], sys.argv[2]
genes = collections.defaultdict(list)
for l in open(faa):
    if l.startswith(">"):
        _, c, a, b, s = l[1:].strip().split("|"); genes[c].append((int(a), int(b), s))
clen = {}
sense = collections.Counter(); spots = {}
for l in open(paf):
    f = l.split("\t"); q, c = f[0], f[5]; clen[c] = int(f[6]); m, al = int(f[9]), int(f[10])
    if al < 60 or m / al < 0.98: continue
    ts, te = int(f[7]), int(f[8]); spot, mate = q.rsplit("/", 1)
    # sense strand of the fragment: R1 antisense under dUTP -> fragment strand = opposite of R1 strand
    fs = ({"+": "-", "-": "+"}[f[4]]) if mate == "1" else f[4]
    spots.setdefault(spot, []).append((c, ts, te, fs))
    if mate == "1":
        g = [x for x in genes[c] if ts >= x[0] and te <= x[1]]
        if len(g) == 1: sense["R1 opposite gene (sense frag)" if fs == g[0][2] else "R1 same as gene (antisense frag)"] += 1
print("R1 orientation within single genes:", dict(sense))
W = 300; dens = collections.Counter(); fragstrand = collections.Counter()
for s, v in spots.items():
    c = v[0][0]; mid = (min(x[1] for x in v) + max(x[2] for x in v)) // 2
    dens[(c, mid // W)] += 1
allw = []
for c, L in clen.items():
    for w in range(L // W): allw.append(dens.get((c, w), 0))
allw.sort(); n = len(allw)
print(f"300-nt windows on contigs with any hit: {n}; spots/window quantiles: 50%={allw[n//2]} 75%={allw[3*n//4]} 90%={allw[9*n//10]} 95%={allw[int(.95*n)]}; fraction of windows with 0 spots={sum(x==0 for x in allw)/n:.3f}")
for c, a, b, lab in (("JAACTV010000069.1", 23650, 24050, "div 16S3'+ITS"), ("JAACTV010000069.1", 20469, 23506, "div 23S"),
                     ("JAACTV010000080.1", 11600, 12238, "norm 16S3'+ITS"), ("JAACTV010000080.1", 12238, 14589, "norm 23S 5' part")):
    v = [s for s, x in spots.items() if any(y[0] == c and y[1] < b and y[2] > a for y in x)]
    st = collections.Counter(next(y[3] for y in spots[s] if y[0] == c) for s in v)
    print(f"{lab} {c}:{a}-{b}: {len(v)} spots, fragment strand {dict(st)}")
