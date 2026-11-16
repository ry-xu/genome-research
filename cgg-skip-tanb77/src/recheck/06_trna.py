"""tRNA context from the lab's ARAGORN calls (codon-atlas trna product): Arg and Gly anticodons per genome,
plus any tRNA with an unusual (non-7-nt) anticodon loop if the column exists."""
import duckdb, sys
W = "recheck_work"
FOCAL = {"G1": "GCA_946496695.1", "G2": "GCF_012317385.1", "G3": "GCA_900553755.1", "G4": "GCA_963931445.1",
         "G5": "GCF_005473905.2", "G6": "GCA_964631135.1", "G7": "GCA_009785865.1", "G8": "GCF_900187165.1",
         "G9": "GCA_949283575.1"}
con = duckdb.connect()
P = lambda acc: f"data/atlas/products/{acc[:8]}/{acc}/trna.parquet"
print(con.sql(f"describe select * from '{P(FOCAL['G1'])}'").fetchall())
for lab, acc in FOCAL.items():


    n = con.sql(f"select count(*) from '{P(acc)}'").fetchall()[0][0]
    arg = con.sql(f"select anticodon, count(*) from '{P(acc)}' where aa ilike 'Arg%' group by 1 order by 1").fetchall()
    gly = con.sql(f"select anticodon, count(*) from '{P(acc)}' where aa ilike 'Gly%' group by 1 order by 1").fetchall()
    other = con.sql(f"select aa, anticodon, count(*) from '{P(acc)}' where aa not ilike 'Arg%' and aa not ilike 'Gly%' and (anticodon ilike '%cg' or length(anticodon) <> 3) group by all order by 1").fetchall()
    print(lab, acc, "n_tRNA", n, "| Arg:", arg, "| Gly:", gly, "| other *CG / non-3 anticodons:", other)
