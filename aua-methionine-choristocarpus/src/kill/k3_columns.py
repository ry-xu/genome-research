"""Kill check 2 + 6: are Choristocarpus ATA@M sites biased in column, gene or position; and what do other genomes
have at the same columns? Column key = Desmarestia NC_007684.1 (never dropped) gene + residue position.
usage: uv run python -I checks/chor_wgs/kill/k3_columns.py"""
from collections import Counter
import polars as pl

df = pl.read_parquet("checks/chor_wgs/kill/sites_all.parquet").filter(pl.col("qpos") > 10)
# consensus per target codon
cons = []
for (t, sid, st, p), g in df.group_by(["target", "contig", "strand", "pos"]):
    res = Counter(g["res"].to_list()); n = g.height
    top, c = res.most_common(1)[0]
    if not (n >= 2 and c / n >= 2 / 3): continue
    d = g.filter(pl.col("ref") == "NC_007684.1")
    dk = (d["qgene"][0].split("_", 1)[1], d["qpos"][0], d["qlen"][0], d["res"][0]) if d.height else (None, None, None, None)
    cons.append((t, sid, st, p, g["codon"][0], top, n, *dk))
C = pl.DataFrame(cons, schema=["target", "contig", "strand", "pos", "codon", "cres", "n", "dgene", "dpos", "dlen", "dres"], orient="row")
C.write_parquet("checks/chor_wgs/kill/consensus.parquet")
T = ["cten_WGS", "dmes_WGS", "POS_NC_085326", "NEG_Phaeostrophion", "NEG_Protohalopteris", "NEG_Saccharina", "NEG_Petalonia", "NEG_Sargassum"]

print("== consensus-M sites: total M columns and all-codon breakdown, per target; aligned consensus codons")
for t in T:
    x = C.filter(pl.col("target") == t)
    m = x.filter(pl.col("cres") == "M"); i = x.filter(pl.col("cres") == "I")
    print(t, "consensus codons", x.height, "M cols", m.height, f"({m.height / x.height:.4f})", "I cols", i.height,
          f"({i.height / x.height:.4f})", "M codons", dict(Counter(m["codon"].to_list()).most_common(6)))

print("\n== relative position in Desmarestia protein (dpos/dlen) of ATA@M (cten) vs ATG@M (each target)")
def qdist(x):
    r = (x["dpos"] / x["dlen"]).drop_nulls().to_list()
    if not r: return "n=0"
    bins = Counter(min(int(v * 5), 4) for v in r)
    return f"n={len(r)} quintiles " + " ".join(str(bins[b]) for b in range(5)) + f" | dpos<=30: {sum(1 for v in x['dpos'].drop_nulls() if v <= 30)}"
cm = C.filter((pl.col("target") == "cten_WGS") & (pl.col("cres") == "M"))
print("cten ATA@M", qdist(cm.filter(pl.col("codon") == "ATA")))
print("cten ATG@M", qdist(cm.filter(pl.col("codon") == "ATG")))
for t in T[1:]:
    print(t, "ATG@M", qdist(C.filter((pl.col("target") == t) & (pl.col("cres") == "M") & (pl.col("codon") == "ATG"))))
print("cten ATT@I", qdist(C.filter((pl.col("target") == "cten_WGS") & (pl.col("cres") == "I") & (pl.col("codon") == "ATT"))))

print("\n== per gene, cten: ATA@M, ATG@M  (genes by Desmarestia name)")
g = cm.filter(pl.col("codon").is_in(["ATA", "ATG"])).group_by("dgene", "codon").len().pivot(on="codon", index="dgene", values="len").fill_null(0).sort("ATA", descending=True)
print(g.to_pandas().to_string(index=False) if False else "\n".join(f"{r['dgene']}\tATA={r['ATA']}\tATG={r['ATG']}" for r in g.iter_rows(named=True)))
print("genes with >=1 ATA@M:", g.filter(pl.col("ATA") > 0).height, "of", g.height, "; max share of ATA@M in one gene:",
      round(g["ATA"].max() / g["ATA"].sum(), 3), "; top-2 share:", round(sum(sorted(g["ATA"].to_list())[-2:]) / g["ATA"].sum(), 3))

print("\n== at cten ATA@M columns (Desmarestia key), codon in other genomes (consensus M or I there)")
key = cm.filter(pl.col("codon") == "ATA").select("dgene", "dpos").drop_nulls()
print("cten ATA@M sites with a Desmarestia key:", key.height)
for t in T[1:]:
    y = C.filter(pl.col("target") == t).join(key, on=["dgene", "dpos"], how="inner")
    print(t, "matched", y.height, "codon", dict(Counter(y["codon"].to_list()).most_common(5)), "cons", dict(Counter(y["cres"].to_list())))
print("\n== at cten ATT@I columns, same view (control)")
key2 = C.filter((pl.col("target") == "cten_WGS") & (pl.col("cres") == "I") & (pl.col("codon") == "ATT")).select("dgene", "dpos").drop_nulls()
for t in ["dmes_WGS", "NEG_Phaeostrophion", "NEG_Saccharina"]:
    y = C.filter(pl.col("target") == t).join(key2, on=["dgene", "dpos"], how="inner")
    print(t, "matched", y.height, "codon", dict(Counter(y["codon"].to_list()).most_common(5)), "cons", dict(Counter(y["cres"].to_list())))
print("\n== Dmes ATG@M columns: what does cten have there?")
key3 = C.filter((pl.col("target") == "dmes_WGS") & (pl.col("cres") == "M") & (pl.col("codon") == "ATG")).select("dgene", "dpos").drop_nulls()
y = C.filter(pl.col("target") == "cten_WGS").join(key3, on=["dgene", "dpos"], how="inner")
print("dmes ATG@M keys", key3.height, "cten matched", y.height, "codon", dict(Counter(y["codon"].to_list()).most_common(6)), "cons", dict(Counter(y["cres"].to_list())))
print("\n== Dmes ATA@I columns: what does cten have there?")
key4 = C.filter((pl.col("target") == "dmes_WGS") & (pl.col("cres") == "I") & (pl.col("codon") == "ATA")).select("dgene", "dpos").drop_nulls()
y = C.filter(pl.col("target") == "cten_WGS").join(key4, on=["dgene", "dpos"], how="inner")
print("dmes ATA@I keys", key4.height, "cten matched", y.height, "codon", dict(Counter(y["codon"].to_list()).most_common(6)), "cons", dict(Counter(y["cres"].to_list())))
