import csv, gzip, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from candidate import read_fasta
COMP = str.maketrans("ACGTN", "TGCAN")
setf, outdir, gdir = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3])
tg = tc = 0
for r in csv.DictReader(open(setf), delimiter="\t"):
    a = r["gca"]
    if not (outdir / f"{a}.fs.tsv.gz").exists():
        continue
    seqs = read_fasta(str(gdir / a[:8] / f"{a}.fna.gz"))
    st = {"+": seqs, "-": {k: v.translate(COMP)[::-1] for k, v in seqs.items()}}
    g = c = 0
    for h in csv.DictReader(gzip.open(outdir / f"{a}.hits.tsv.gz", "rt"), delimiter="\t"):
        if h["n_fs"] != "0":
            continue
        L = len(seqs[h["contig"]]); g0, g1 = int(h["gstart"]), int(h["gend"])
        lo, hi = min(g0, g1) - 1, max(g0, g1)
        s = st[h["strand"]][h["contig"]]; p0 = lo if h["strand"] == "+" else L - hi
        for p in range(p0 + 60, p0 + (hi - lo) - 60, 3):
            if s[p:p + 2] == "GG":
                g += s[p - 1] != "C"; c += s[p - 1] == "C"
    tg += g; tc += c
print("in-frame GG[ACGT] codons (frame-continuous hits, >=20 codons from ends): not after C", tg, "after C", tc)
