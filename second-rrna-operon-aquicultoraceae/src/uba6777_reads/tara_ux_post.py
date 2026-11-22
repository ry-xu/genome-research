"""Post-hoc diagnostics (written 21:45 14 Nov AFTER the full-run Rule 1 result: 10856 end 9 spanning reads, 0.49x the
control median, narrowly failing the >= 0.5x rule). Not pre-registered; reported as supplementary.

usage: python3 -I tara_ux_post.py HITS.sam
1. Depth (reads at >= 0.98 overlapping each 50-nt bin) across 6000-6700 and 10500-11200 of DKGF01000210.1.
2. Mate bridges: read names with one mate (>= 0.98) ending at least 20 nt inside the 23S and the other mate starting at
   least 20 nt outside it, on the same contig, for each end; same count for a virtual boundary at each control position.
   (Gave 0 everywhere, controls included: mates mostly overlap, so this is uninformative.)
3. Fragment spans: for pairs with both mates at >= 0.98 on the contig, the fragment = first mate start to last mate end;
   count fragments covering each end +-15 nt vs the control positions; also the fragment length distribution.
"""
import collections
import re
import sys

CTG, D5, D3 = "MAG|DKGF01000210.1", 6356, 10856
CTRL = (2000, 3000, 4000, 5000, 12500, 14000, 15500, 17000, 18500, 20000)


def ident_end(f):
    ops = re.findall(r"(\d+)([MIDNSHP=X])", f[5])
    m = sum(int(n) for n, o in ops if o in "M=X")
    i = sum(int(n) for n, o in ops if o == "I")
    d = sum(int(n) for n, o in ops if o == "D")
    nm = int(re.search(r"NM:i:(\d+)", "\t".join(f[11:])).group(1))
    return 1 - nm / max(1, m + i + d), int(f[3]) + m + d - 1


def main(sam):
    reads = collections.defaultdict(list)
    for line in open(sam):
        if line[0] == "@":
            continue
        f = line.split("\t")
        if f[2] != CTG or int(f[1]) & (4 | 256 | 2048):
            continue
        idn, end = ident_end(f)
        if idn >= 0.98:
            reads[f[0]].append((int(f[3]), end))
    depth = collections.Counter()
    for v in reads.values():
        for s, e in v:
            for b in range(s // 50, e // 50 + 1):
                depth[b * 50] += 1
    for lo, hi in ((6000, 6700), (10500, 11200)):
        print(f"depth {lo}-{hi} (50-nt bins):", " ".join(f"{b}:{depth[b]}" for b in range(lo, hi, 50)))

    def bridges(x):
        n = 0
        for v in reads.values():
            if len(v) == 2:
                (a, b) = sorted(v)
                if a[1] <= x - 20 and b[0] >= x + 20:
                    n += 1
        return n
    print("mate bridges: 23S 3' end (6356):", bridges(D5), " 23S 5' end (10856):", bridges(D3))
    c = {x: bridges(x) for x in CTRL}
    print("controls:", c, "median", sorted(c.values())[len(c) // 2 - 1:len(c) // 2 + 1])
    print("pairs with both mates on contig:", sum(len(v) == 2 for v in reads.values()))
    frags = [(min(s for s, _ in v), max(e for _, e in v)) for v in reads.values() if len(v) == 2]
    fl = sorted(e - s + 1 for s, e in frags)
    print("fragment length quartiles:", fl[len(fl) // 4], fl[len(fl) // 2], fl[3 * len(fl) // 4])
    cov = lambda x: sum(1 for s, e in frags if s <= x - 15 and e >= x + 15)
    c = sorted(cov(x) for x in CTRL)
    print(f"fragments spanning: 6356: {cov(D5)}, 10856: {cov(D3)}; controls {c}, median {(c[4] + c[5]) / 2}")


if __name__ == "__main__":
    main(sys.argv[1])
