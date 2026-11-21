"""List 23S census candidates under the rules in leads/lsu_23s.md (written before the bacterial census was read).

Usage: python -I src/lsu_candidates.py FIXED.tsv DEVIATIONS.parquet OUT.tsv [--bar 0.999]
A candidate is a fixed.tsv row whose column is >= bar conserved (f_cons) or that sits at a named site (site column).
Per candidate: the carrier genomes (deviation rows with that base in that clade), their distinct NCBI BioProjects, and
whether it is a macrolide site (2058/2059/2611; reported as resistance biology). The row index is fixed.tsv's 0-based
data row, which is what lsu_copies.py --rows takes. Rules 1-2 (raw base, all copies) are lsu_copies.py's job.
"""
import os
import argparse

import duckdb

META = os.path.join(os.environ.get("FROZEN_CORE", "frozen-core"), "results/meta.parquet")
MACROLIDE = {"2058", "2059", "2611"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fixed"); ap.add_argument("devs"); ap.add_argument("out")
    ap.add_argument("--bar", type=float, default=0.999)
    a = ap.parse_args()
    fx = duckdb.sql(f"select * from read_csv('{a.fixed}', delim='\t', header=true, all_varchar=true)").df()
    n = 0
    with open(a.out, "w") as out:
        out.write("row\tclade\tecoli_pos\tnt\tcons\tf_cons\tk\tn\tsite\tkind\tmacrolide\tcarriers\tbioprojects\n")
        for i, r in fx.iterrows():
            site = r.site if isinstance(r.site, str) else ""
            if float(r.f_cons) < a.bar and not site:
                continue
            accs = duckdb.sql(
                f"select d.acc, m.ncbi_bioproject from '{a.devs}' d join '{META}' m using(acc) "
                f"where d.col = {int(r.col)} and d.nt = '{r.nt}' "
                f"and (m.gtdb_taxonomy like '%{r.clade};%' or m.gtdb_taxonomy like '%{r.clade}')").fetchall()
            bp = {b for _, b in accs if b and b != "none"}
            out.write(f"{i}\t{r.clade}\t{r.ecoli_pos}\t{r.nt}\t{r.cons}\t{r.f_cons}\t{r.k}\t{r.n}\t{site}\t{r.kind}\t"
                      f"{'yes' if r.ecoli_pos in MACROLIDE else ''}\t{len(accs)}\t{len(bp)}\n")
            n += 1
    print("candidates", n, "of", len(fx), "fixed rows")


if __name__ == "__main__":
    main()
