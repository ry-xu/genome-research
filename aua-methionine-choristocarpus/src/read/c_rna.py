"""RNA-editing kill test at ATA/ATG sites (RULES.md Amendment C). Reads SAM on stdin (minimap2 -ax sr to the mito contig).
usage: minimap2 -ax sr -t 1 MITO.fa READS.fq.gz | uv run python -I checks/chor_wgs/c_rna.py SITES.tsv OUT.tsv"""
import re, statistics, sys
from collections import Counter, defaultdict
import polars as pl

sites = pl.read_csv(sys.argv[1], separator="\t").filter(
    pl.col("codon").is_in(["ATA", "ATG"]) & (pl.col("ref_res").is_in(["M", "I"])))
COMP = str.maketrans("ACGT", "TGCA")
want = {}                                                  # forward coord of codon base 3 -> (group, strand)
for r in sites.iter_rows(named=True):
    p3 = r["pos"] + 2 if r["strand"] == "+" else r["pos"] - 2
    want[p3] = (f"{r['codon']}@{r['ref_res']}", r["strand"])
cnt = defaultdict(Counter)
for l in sys.stdin:
    if l[0] == "@":
        continue
    f = l.split("\t", 11)
    flag, pos, mapq, cig, seq, qual = int(f[1]), int(f[3]), int(f[4]), f[5], f[9], f[10]
    if flag & 0x904 or mapq < 20 or cig == "*":
        continue
    rp, qp = pos, 0
    for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cig):
        n = int(n)
        if op in "M=X":
            for k in range(n):
                if rp + k in want and ord(qual[qp + k]) - 33 >= 20:
                    b = seq[qp + k]
                    if want[rp + k][1] == "-":
                        b = b.translate(COMP)
                    cnt[rp + k][b] += 1
            rp += n; qp += n
        elif op in "IS":
            qp += n
        elif op in "DN":
            rp += n
rows = [(p, g, s, sum(cnt[p].values()), cnt[p]["A"], cnt[p]["G"], cnt[p]["C"], cnt[p]["T"]) for p, (g, s) in want.items()]
t = pl.DataFrame(rows, schema=["pos3", "group", "strand", "depth", "A", "G", "C", "T"], orient="row").sort("pos3")
t.write_csv(sys.argv[2], separator="\t")
for g in ["ATA@M", "ATA@I", "ATG@M", "ATG@I"]:
    x = t.filter((pl.col("group") == g) & (pl.col("depth") >= 10))
    gs = [r["G"] / r["depth"] for r in x.iter_rows(named=True)]
    print(g, "sites", t.filter(pl.col("group") == g).height, "depth>=10", x.height,
          "median G share", round(statistics.median(gs), 3) if gs else None,
          "max G share", round(max(gs), 3) if gs else None, "median depth", int(x["depth"].median()) if x.height else None)
