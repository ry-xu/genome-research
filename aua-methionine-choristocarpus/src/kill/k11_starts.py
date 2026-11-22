"""Kill check extra: annotated start codons of NC_085326 CDS vs the Sphacelariales/Discosporangiales-adjacent controls.
usage: uv run python -I checks/chor_wgs/kill/k11_starts.py"""
import re
from collections import Counter
import polars as pl
from Bio.Seq import Seq
RD = "scratchpad/chor_work/read"
cds = pl.read_parquet("data/cds_mito.parquet", columns=["acc", "gene", "loc"])
for label, acc in [("POS_NC_085326", "NC_085326.1"), ("NEG_Saccharina", "NC_013473.1"), ("NEG_Petalonia", "NC_025227.1"), ("NEG_Sargassum", "NC_024613.1")]:
    g = "".join(l.strip().upper() for l in open(f"{RD}/{label}_target.fa") if not l.startswith(">"))
    st = Counter()
    for r in cds.filter(pl.col("acc") == acc).iter_rows(named=True):
        m = re.fullmatch(r"(complement\()?(\d+)\.\.(\d+)\)?", r["loc"] or "")
        if not m: continue
        nt = g[int(m.group(2)) - 1:int(m.group(3))]
        if m.group(1): nt = str(Seq(nt).reverse_complement())
        st[nt[:3]] += 1
    print(label, "annotated start codons", dict(st.most_common()))
