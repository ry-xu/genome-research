"""Pile up minimap2 --cs=long short-read alignments over chosen contig regions; report depth and disagreements with the reference.

usage: python -I src/contig_pileup.py REF.fa PAF CONTIG:START-END[,CONTIG:START-END...] [--min-id 0.98] [--min-mapq 20]
  Prints, per region, mean/min depth of passing alignments (identity = matches / block length, mapq >= min) and every
  position whose reference base is < 90% of passing calls. Also prints mean depth per contig over all passing alignments.
  Coordinates are 1-based inclusive.
"""
import collections
import re
import sys


def read_fasta(path):
    out, name = {}, None
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].split()[0]
            out[name] = []
        else:
            out[name].append(line.strip())
    return {k: "".join(v).upper() for k, v in out.items()}


def main():
    ref, paf, regs = sys.argv[1], sys.argv[2], sys.argv[3]
    min_id = float(sys.argv[sys.argv.index("--min-id") + 1]) if "--min-id" in sys.argv else 0.98
    min_q = int(sys.argv[sys.argv.index("--min-mapq") + 1]) if "--min-mapq" in sys.argv else 20
    g = read_fasta(ref)
    regions = []
    for r in regs.split(","):
        c, se = r.split(":")
        s, e = map(int, se.split("-"))
        regions.append((c, s, min(e, len(g[c]))))
    calls = {i: collections.defaultdict(collections.Counter) for i in range(len(regions))}
    bases = collections.Counter()
    for line in open(paf):
        f = line.rstrip("\n").split("\t")
        tn, ts, te = f[5], int(f[7]), int(f[8])
        if int(f[11]) < min_q or int(f[9]) / int(f[10]) < min_id:
            continue
        bases[tn] += te - ts
        hit = [i for i, (c, s, e) in enumerate(regions) if c == tn and ts < e and te >= s]
        if not hit:
            continue
        cs = next(x[5:] for x in f[12:] if x.startswith("cs:Z:"))
        cur = ts
        for op, s in re.findall(r"([=*+\-~:])([A-Za-z0-9]+)", cs):
            if op == "=":
                for k, b in enumerate(s.upper()):
                    for i in hit:
                        calls[i][cur + k + 1][b] += 1
                cur += len(s)
            elif op == "*":
                for i in hit:
                    calls[i][cur + 1][s[1].upper()] += 1
                cur += 1
            elif op == "-":
                for k in range(len(s)):
                    for i in hit:
                        calls[i][cur + k + 1]["-"] += 1
                cur += len(s)
            elif op == "+":
                for i in hit:
                    calls[i][cur]["+"] += 1
    for c in g:
        if bases[c]:
            print(f"contig {c} len {len(g[c])} mean depth {bases[c] / len(g[c]):.1f}")
    for i, (c, s, e) in enumerate(regions):
        d = [sum(v for k, v in calls[i][p].items() if k != "+") for p in range(s, e + 1)]
        print(f"region {c}:{s}-{e} mean {sum(d) / len(d):.1f} min {min(d)} max {max(d)} zero-depth {sum(x == 0 for x in d)}")
        for p in range(s, e + 1):
            t = calls[i][p]
            n = sum(v for k, v in t.items() if k != "+")
            rb = g[c][p - 1]
            if n and (t[rb] / n < 0.9 or t["+"] / n >= 0.1):
                print(f"  {c}:{p} ref {rb} calls {dict(t)}")


main()
