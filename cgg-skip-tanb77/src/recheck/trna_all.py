import duckdb, sys
acc=sys.argv[1]
con=duckdb.connect()
print(con.sql(f"select aa, anticodon, count(*) n from 'data/atlas/products/{acc[:8]}/{acc}/trna.parquet' group by all order by aa, anticodon").fetchall())
