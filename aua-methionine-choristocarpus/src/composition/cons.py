"""Brown-alga consensus residue (Lyra kill consensus.parquet, cres) at Choristocarpus codons, by gene class.
Respiratory = nad/cox/cob/atp; ribo+tatC = rps/rpl/tatC. Controls (NEG_*) pooled for the same codons."""
import polars as pl
d = pl.read_parquet("lab/work/organelle-reading/checks/chor_wgs/kill/consensus.parquet")
d = d.with_columns(pl.when(pl.col("dgene").str.to_lowercase().str.contains("^(rps|rpl|tatc)")).then(pl.lit("ribo+tatC"))
                   .otherwise(pl.lit("resp")).alias("cls"),
                   pl.when(pl.col("target").str.starts_with("NEG")).then(pl.lit("NEG")).otherwise(pl.col("target")).alias("t"))
print(d["target"].unique().to_list())
for t in ["POS_NC_085326", "NEG"]:
    for cls in ["resp", "ribo+tatC"]:
        for cod in ["ATA", "ATG", "ATT"]:
            s = d.filter((pl.col("t") == t) & (pl.col("cls") == cls) & (pl.col("codon") == cod))
            n = s.height
            top = s.group_by("cres").len().sort("len", descending=True).head(6)
            print(t, cls, cod, n, " ".join(f"{r}:{c}({c/n:.2f})" for r, c in top.iter_rows()) if n else "")
