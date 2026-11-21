"""Do bac120 residues change with a 23S variant across independent origins? (plan: leads/lsu_23s.md, 14 Nov)

Usage: python -I src/srl_covar.py SET OUT.tsv [--perms 10000] [--seed 1]
SET is json:KEY (a final list in results/lsu/bac/kill/work/recount.json), a kill-agent set name from results/lsu/bac/kill/carrier_sets.tsv (rows whose status starts "all copies"), or
lsu:POS:NT for every census genome with base NT at E. coli position POS (from results/lsu/bac/deviations.parquet).
Genomes: the bacterial 23S census (data/lsu/bac_lsu*.fna.tsv). Background: census genomes with no deviation at POS
(and, for the SRL sets, none at 2659-2662). Per bac120 column (f_cons >= 0.99): k = distinct GTDB orders with >= 1
set genome carrying a non-consensus residue (gaps ignored). Null: within strata (order x above/below that order's
median genome-wide deviation count) draw as many genomes as the set has there; p = (1 + #null >= obs) / (1 + perms).
Prints Bonferroni passers (over columns with obs >= 2) and how many fall in the three SRL-binding GTPases.
"""
import os
import argparse
import json

import duckdb
import numpy as np

W = os.environ.get("FROZEN_CORE", "frozen-core")  # lab working tree of the census
GTP = {"TIGR00487": "infB", "TIGR01393": "lepA", "TIGR01394": "bipA"}
SRL = (2659, 2660, 2661, 2662)


def lsu_col(pos):
    return duckdb.sql(f"select col from read_csv('{W}/results/lsu/bac/columns.tsv', delim='\t', header=true) "
                      f"where ecoli_pos = {pos} and ecoli_nt != '-'").fetchall()[0][0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("set"); ap.add_argument("out")
    ap.add_argument("--perms", type=int, default=10000); ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    census = [r[0] for r in duckdb.sql(f"select acc from read_csv('{W}/data/lsu/bac_lsu*.fna.tsv', delim='\t', "
                                       f"header=true, union_by_name=true)").fetchall()]
    ldev = f"{W}/results/lsu/bac/deviations.parquet"
    if a.set.startswith("lsu:"):
        _, pos, nt = a.set.split(":"); pos = int(pos)
        target = {r[0] for r in duckdb.sql(f"select acc from '{ldev}' where col = {lsu_col(pos)} and nt = '{nt}'").fetchall()}
        excl = {r[0] for r in duckdb.sql(f"select acc from '{ldev}' where col = {lsu_col(pos)}").fetchall()}
    elif a.set.startswith("json:"):  # json:KEY = a final list in kill/work/recount.json (e.g. json:GAAA tetraloop (2661A list))
        target = set(json.load(open(f"{W}/results/lsu/bac/kill/work/recount.json"))[a.set[5:]]["full"])
        pos = 2661
        excl = set()
    else:
        target = {r[0] for r in duckdb.sql(f"select acc from read_csv('{W}/results/lsu/bac/kill/carrier_sets.tsv', "
                                           f"delim='\t', header=true) where \"set\" = '{a.set}' "
                                           f"and status like 'all copies%'").fetchall()}
        pos = 2661
        excl = set()
    if pos in SRL or a.set.startswith("2661") or a.set.startswith("json:") or a.set.startswith("2660"):
        cols = ",".join(str(lsu_col(p)) for p in SRL)
        excl |= {r[0] for r in duckdb.sql(f"select acc from '{ldev}' where col in ({cols})").fetchall()}
    genomes = sorted(set(census) - (excl - target))
    target &= set(genomes)
    gi = {g: i for i, g in enumerate(genomes)}
    meta = dict(duckdb.sql(f"select acc, \"order\" from '{W}/results/meta.parquet'").fetchall())
    cols = duckdb.sql(f"select marker, col, cons from '{W}/results/bac120/columns.parquet' where f_cons >= 0.99 "
                      f"order by marker, col").fetchall()
    ci = {(m, c): j for j, (m, c, _) in enumerate(cols)}
    dev = duckdb.sql(f"select acc, marker, col, res from '{W}/results/bac120/deviations.parquet' "
                     f"where res != '-' and acc in (select unnest(?::varchar[]))", params=[genomes]).fetchall()
    M = np.zeros((len(genomes), len(cols)), dtype=bool)
    resA = {}
    for acc, m, c, res in dev:
        j = ci.get((m, c))
        if j is None:
            continue
        M[gi[acc], j] = True
        if acc in target:
            resA.setdefault(j, []).append(res)
    devcount = M.sum(1)
    order = np.array([meta.get(g) or "none" for g in genomes])
    strata = {}
    for o in set(order[[gi[g] for g in target]]):
        idx = np.where(order == o)[0]
        med = np.median(devcount[idx])
        for hi in (False, True):
            s = idx[(devcount[idx] > med) == hi]
            t = [i for i in s if genomes[i] in target]
            if t:
                strata[(o, hi)] = (s, len(t))
    tidx = np.array(sorted(gi[g] for g in target))

    def k_of(rows):
        rows = rows[np.argsort(order[rows], kind="stable")]
        o = order[rows]
        starts = np.r_[0, np.where(o[1:] != o[:-1])[0] + 1]
        return np.logical_or.reduceat(M[rows], starts, axis=0).sum(0)

    obs = k_of(tidx)
    rng = np.random.default_rng(a.seed)
    ge = np.zeros(len(cols), dtype=np.int64); tot = np.zeros(len(cols))
    for _ in range(a.perms):
        rows = np.concatenate([rng.choice(s, n, replace=False) for s, n in strata.values()])
        k = k_of(rows)
        ge += k >= obs; tot += k
    p = (1 + ge) / (1 + a.perms)
    colmap = {(m, c): (g, rp, rr) for m, c, g, rp, rr in duckdb.sql(
        f"select marker, col, gene, ref_pos, ref_res from '{W}/results/bac120/colmap.parquet'").fetchall()}
    tested = obs >= 2
    bonf = 0.05 / max(1, tested.sum())
    with open(a.out, "w") as out:
        out.write("marker\tgene\tcol\tref_pos\tref_res\tcons\tobs_orders\tnull_mean\tp\tbonferroni\tset_residues\n")
        for j in np.argsort(p):
            if not tested[j]:
                continue
            m, c, cons = cols[j]
            g, rp, rr = colmap.get((m, c), ("", "", ""))
            r = resA.get(j, [])
            rs = ",".join(f"{x}{r.count(x)}" for x in sorted(set(r)))
            out.write(f"{m}\t{g}\t{c}\t{rp}\t{rr}\t{cons}\t{obs[j]}\t{tot[j] / a.perms:.2f}\t{p[j]:.5f}\t"
                      f"{'yes' if p[j] < bonf else ''}\t{rs}\n")
    passers = [j for j in np.where(tested & (p < bonf))[0]]
    gtp_tested = sum(cols[j][0] in GTP for j in np.where(tested)[0])
    print(f"set {a.set}: {len(tidx)} genomes in {len({order[i] for i in tidx})} orders; strata {len(strata)}; "
          f"background {len(genomes) - len(tidx)}")
    print(f"columns tested (obs>=2) {tested.sum()} (GTPase {gtp_tested}); Bonferroni p < {bonf:.2e}: "
          f"{len(passers)} (GTPase {sum(cols[j][0] in GTP for j in passers)})")
    print(f"set genomes' mean deviations {devcount[tidx].mean():.1f} vs all {devcount.mean():.1f}")


if __name__ == "__main__":
    main()
