"""Kill check 5: where do NC_085326 (Choi, SPAdes) and CAXNAD020000071.1 (Genoscope) differ?
blastn (megablast, 1 thread) of NC_085326 vs the WGS contig; walk every HSP; list mismatches and indels with NC position,
feature (from NC_085326 GenBank), homopolymer context; and whether any lies in a codon read at an M/I site.
usage: uv run python -I checks/chor_wgs/kill/k5_diff.py GBDIR"""
import subprocess, sys
from collections import Counter
import polars as pl
from Bio import SeqIO

GB = sys.argv[1]
RD = "scratchpad/chor_work"
nc = SeqIO.read(f"{GB}/NC_085326.1.gb", "genbank")
q, s = f"{RD}/NC_085326.fa", f"{RD}/cten_mito.fa"
out = subprocess.run(["blastn", "-query", q, "-subject", s, "-evalue", "1e-50", "-dust", "no",
                      "-outfmt", "6 qstart qend sstart send pident length mismatch gaps qseq sseq"], check=True, capture_output=True, text=True).stdout
feat = []
for f in nc.features:
    if f.type in ("CDS", "tRNA", "rRNA"):
        feat.append((int(f.location.start) + 1, int(f.location.end), f.type, f.qualifiers.get("gene", ["?"])[0]))
def where(p):
    h = [f"{t}:{g}" for a, b, t, g in feat if a <= p <= b]
    return ",".join(h) if h else "intergenic"
ncs = str(nc.seq).upper()
def homo(p):  # length of homopolymer run touching NC position p (1-based)
    b = ncs[p - 1]; i = j = p - 1
    while i > 0 and ncs[i - 1] == b: i -= 1
    while j < len(ncs) - 1 and ncs[j + 1] == b: j += 1
    return j - i + 1
diffs = []
print("HSPs (qstart qend sstart send pident length mismatch gaps):")
covered = set()
for l in out.splitlines():
    qs, qe, ss, se, pid, L, mm, gp, qa, sa = l.split("\t")
    print(" ", qs, qe, ss, se, pid, L, mm, gp)
    qp = int(qs) - 1
    for a, b in zip(qa, sa):
        if a != "-": qp += 1
        if a != "-" and qp in covered: continue
        if a != b:
            diffs.append((qp, a, b, "ins_in_WGS" if a == "-" else "del_in_WGS" if b == "-" else "subst", where(qp), homo(qp)))
    covered.update(range(int(qs), int(qe) + 1))
print("NC positions covered by HSPs", len(covered), "of", len(ncs))
D = pl.DataFrame(diffs, schema=["nc_pos", "nc", "wgs", "kind", "feature", "homopolymer_len"], orient="row").unique().sort("nc_pos")
print("differences", D.height, dict(Counter(D["kind"].to_list())))
print("in CDS", D.filter(pl.col("feature").str.contains("CDS")).height, "| in homopolymer >= 4:", D.filter(pl.col("homopolymer_len") >= 4).height)
for r in D.iter_rows(): print("  ", *r)
# overlap with read M/I codons in the NC control
si = pl.read_parquet("checks/chor_wgs/kill/consensus.parquet").filter((pl.col("target") == "POS_NC_085326") & pl.col("cres").is_in(["M", "I"]))
cod = set()
for r in si.iter_rows(named=True):
    rng = range(r["pos"], r["pos"] + 3) if r["strand"] == "+" else range(r["pos"] - 2, r["pos"] + 1)
    cod.update(rng)
print("differences inside a consensus M/I codon of NC_085326:", D.filter(pl.col("nc_pos").is_in(list(cod))).height)
