"""Kill check 7: RNA editing instead of a code change? Choristocarpus tenellus KU-1152 RNA-seq (ERR12954301, read 1 only,
first ~150 MB of the ENA fastq.gz) mapped with minimap2 (-t 1, -ax sr) to NC_085326. At each consensus site (k3 consensus
table, target POS_NC_085326) count read bases (MAPQ >= 20, base Q >= 20) at the three codon positions. If ATA@M were A->G
edited to ATG in the transcript, the third codon base would read G (plus-strand sites) or C (minus-strand sites, on the
plus-strand coordinate). Controls: ATG@M, ATT@I, ATA@I sites and the genome-wide mismatch spectrum.
usage: uv run python -I checks/chor_wgs/kill/k12_rnaseq.py RNA_DIR REF_FA"""
import gzip, os, re, subprocess, sys, zlib
from collections import Counter, defaultdict
import polars as pl

RNA, REF = sys.argv[1:3]
ref = "".join(l.strip().upper() for l in open(REF) if not l.startswith(">"))
fq = os.path.join(RNA, "reads_head.fq")
n = 0
with open(fq, "w") as out:
    for gz in sorted(x for x in os.listdir(RNA) if x.endswith("_head.fastq.gz")):
        k = 0
        with gzip.open(os.path.join(RNA, gz), "rt") as fh:
            buf = []
            try:
                for line in fh:
                    buf.append(line)
                    if len(buf) == 4:
                        if buf[0].startswith("@") and buf[2].startswith("+") and len(buf[1]) == len(buf[3]): out.write("".join(buf)); k += 1
                        buf = []
            except (EOFError, zlib.error, OSError):
                pass
        print(f"{gz} ({os.path.getsize(os.path.join(RNA, gz))} bytes): complete reads from the truncated gzip: {k}")
        n += k
print(f"reads total: {n}")
sam = os.path.join(RNA, "reads_vs_NC085326.sam")
subprocess.run(["minimap2", "-t", "1", "-ax", "sr", "--secondary=no", "-o", sam, REF, fq], check=True, capture_output=True)
cnt = defaultdict(Counter)   # 0-based ref pos -> base counts (plus strand)
mapped = q20 = 0
for l in open(sam):
    if l.startswith("@"): continue
    f = l.split("\t"); flag, pos, mq, cig, seq, qual = int(f[1]), int(f[3]) - 1, int(f[4]), f[5], f[9], f[10]
    if flag & 4 or flag & 256 or flag & 2048: continue
    mapped += 1
    if mq < 20: continue
    q20 += 1
    rp, qp = pos, 0
    for L, op in re.findall(r"(\d+)([MIDNSHX=])", cig):
        L = int(L)
        if op in "M=X":
            for i in range(L):
                if ord(qual[qp + i]) - 33 >= 20 and 0 <= rp + i < len(ref): cnt[rp + i][seq[qp + i]] += 1
            rp += L; qp += L
        elif op in "IS": qp += L
        elif op in "DN": rp += L
print(f"primary mapped reads: {mapped}; MAPQ>=20: {q20} ({q20 / n:.4f} of reads)")
cov = [sum(cnt[i].values()) for i in range(len(ref))]
print(f"positions with depth >= 10: {sum(1 for c in cov if c >= 10)}, >= 3: {sum(1 for c in cov if c >= 3)} of {len(ref)}; median depth {sorted(cov)[len(cov) // 2]}, mean {sum(cov) / len(cov):.1f}")
spec = Counter(); npos = 0
for i, c in cnt.items():
    d = sum(c.values())
    if d < 5: continue
    npos += 1
    alt, ac = max(((b, k) for b, k in c.items() if b != ref[i]), key=lambda x: x[1], default=(None, 0))
    if alt and ac / d >= 0.1: spec[f"{ref[i]}>{alt}" + ("" if ac / d < 0.5 else " (major)")] += 1
print(f"genome-wide: positions (depth>=5) with a non-reference base at >=10% of reads: {sum(spec.values())} of {npos}: {dict(spec.most_common())}")

C = pl.read_parquet("checks/chor_wgs/kill/consensus.parquet").filter(pl.col("target") == "POS_NC_085326")
COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}
def codon_cols(p, strand):
    """0-based plus-strand coordinates of codon positions 1,2,3 (k9 slicing: + g[p-1:p+2]; - revcomp g[p-3:p])"""
    return [p - 1, p, p + 1] if strand == "+" else [p - 1, p - 2, p - 3]
print("\n== consensus sites: sense-strand read bases at codon position 3, pooled over sites (Q>=20, MAPQ>=20); per-site summary for depth >= 3")
for codon, res in [("ATA", "M"), ("ATG", "M"), ("ATT", "I"), ("ATA", "I"), ("ATC", "I")]:
    s = C.filter((pl.col("codon") == codon) & (pl.col("cres") == res))
    rows = []
    for r in s.iter_rows(named=True):
        cols = codon_cols(r["pos"], r["strand"]); per = []
        for j, x in enumerate(cols):
            c = cnt.get(x, Counter())
            if r["strand"] == "-": c = Counter({COMP.get(b, "N"): k for b, k in c.items()})
            per.append(c)
        g = "".join(ref[x] if r["strand"] == "+" else COMP[ref[x]] for x in cols)
        assert g == codon, (r, g)
        rows.append((r["pos"], r["strand"], sum(per[2].values()), per))
    pool = [Counter(), Counter(), Counter()]
    for _, _, _, per in rows:
        for j in range(3): pool[j].update(per[j])
    tot = sum(pool[2].values())
    cov_sites = [x for x in rows if x[2] >= 1]; d3 = [x for x in rows if x[2] >= 3]
    gsite = [x for x in d3 if x[3][2]["G"] / x[2] >= 0.2]
    print(f"{codon}@{res}: {s.height} sites, {len(cov_sites)} with any read, {len(d3)} with depth >= 3; pooled 3rd-base counts {dict(pool[2].most_common())}"
          f" (G fraction {pool[2]['G'] / tot if tot else float('nan'):.3f}, genome base '{codon[2]}' fraction {pool[2][codon[2]] / tot if tot else float('nan'):.3f});"
          f" pooled 1st {dict(pool[0].most_common(3))} 2nd {dict(pool[1].most_common(3))}; depth>=3 sites with G >= 20%: {len(gsite)}")
    if codon == "ATA" and res == "M":
        for pos, st, d, p in sorted(d3, key=lambda x: (-x[3][2]["G"] / x[2], -x[2]))[:5]:
            print(f"   {st} pos {pos} depth {d} third-base {dict(p[2])}")
