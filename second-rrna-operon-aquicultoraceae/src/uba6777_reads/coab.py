"""Co-abundance test for NODE_6 vs the UBA6777 host (rules in ux/NOTES.md, "Co-abundance test", fixed before mapping).

usage: python3 -I coab.py RUN1.ref4.sam RUN2.ref4.sam ...
Per sample, at identity >= 0.98 (1 - NM/(M+I+D)), primary alignments only: H = median reads/kb over MAG contigs >= 5 kb,
D = reads/kb on the divergent 23S (MAG|DKGF01000210.1:6356-10856, by overlap), and reads/kb of every ASM|NODE contig.
Prints the table, fold ranges (max/min over samples) of D/H, N6/H and each X/H, and the pre-registered reading.
"""
import collections
import re
import statistics
import sys

DIV = ("MAG|DKGF01000210.1", 6356, 10856)
N6 = "ASM|NODE_6"


def sample(sam):
    ln, n, dv = {}, collections.Counter(), 0
    for line in open(sam):
        f = line.split("\t")
        if line[0] == "@":
            if line.startswith("@SQ"):
                ln[f[1][3:]] = int(f[2][3:])
            continue
        if int(f[1]) & (4 | 256 | 2048):
            continue
        ops = re.findall(r"(\d+)([MIDNSHP=X])", f[5])
        m = sum(int(a) for a, o in ops if o in "M=X")
        i = sum(int(a) for a, o in ops if o == "I")
        d = sum(int(a) for a, o in ops if o == "D")
        nm = int(next(t[5:] for t in f[11:] if t.startswith("NM:i:")))
        if 1 - nm / max(1, m + i + d) < 0.98:
            continue
        n[f[2]] += 1
        s = int(f[3])
        if f[2] == DIV[0] and s <= DIV[2] and s + m + d - 1 >= DIV[1]:
            dv += 1
    host = statistics.median(n[k] / ln[k] * 1000 for k in ln if k.startswith("MAG|") and ln[k] >= 5000)
    asm = {k: n[k] for k in ln if k.startswith("ASM|")}
    return host, dv / (DIV[2] - DIV[1] + 1) * 1000, {k: (c, c / ln[k] * 1000) for k, c in asm.items()}, dv


def fold(xs):
    return max(xs) / min(xs) if min(xs) > 0 else float("inf")


RUNS = {"ERR599038", "ERR599030", "ERR599136"}


def main(sams):
    names = [p.split("/")[-1].split(".")[0] for p in sams]
    assert sorted(names) == sorted(RUNS), f"need exactly {sorted(RUNS)} once each, got {names}"
    S = [(nm, *sample(p)) for nm, p in zip(names, sams)]
    assert all(h > 0 for _, h, _, _, _ in S), "host median 0 in a sample"
    for name, h, d, asm, dv in S:
        print(f"{name}\tH {h:.1f}/kb\tD {d:.1f}/kb ({dv} reads)\tD/H {d / h:.3f}\tN6 {asm[N6][1]:.1f}/kb ({asm[N6][0]} reads)\tN6/H {asm[N6][1] / h:.3f}")
    fd = fold([d / h for _, h, d, _, _ in S])
    f6 = fold([a[N6][1] / h for _, h, _, a, _ in S])
    dec = {k: fold([a[k][1] / h for _, h, _, a, _ in S]) for k in S[0][3]
           if k != N6 and all(a[k][0] >= 20 for _, _, _, a, _ in S)}
    print(f"fold range D/H {fd:.3f}; N6/H {f6:.3f}")
    print("decoys (>= 20 reads everywhere):", len(dec), "fold ranges:", " ".join(f"{k[4:]}:{v:.3f}" for k, v in sorted(dec.items(), key=lambda x: x[1])))
    nvar = sum(v > 1.5 for v in dec.values())
    better = sum(f6 < v for v in dec.values())
    print(f"decoys with fold > 1.5: {nvar}; N6 smaller than {better}/{len(dec)} decoys")
    absent = [nm for nm, h, d, a, dv in S if a[N6][0] >= 20 and dv < 5]
    # rule order as amended in ux/NOTES.md (written before ERR599030/ERR599136 were mapped)
    if absent:
        print("READING: NODE_6 is not the host's (N6 >= 20 reads where D < 5):", absent)
    elif fd > 1.5:
        print("POSITIVE CONTROL FAILED (D/H fold > 1.5): no reading")
    elif nvar < 10:
        print("UNINFORMATIVE: fewer than 10 decoys vary > 1.5-fold")
    elif f6 > 3:
        print("READING: NODE_6 is not the host's (N6/H fold > 3)")
    elif f6 <= 1.5 and better >= 0.9 * len(dec):
        print("READING: NODE_6 tracks the host (co-occurrence, not linkage)")
    else:
        print("READING: inconclusive")


if __name__ == "__main__":
    main(sys.argv[1:])
