"""Where in GTDB does CGG have no tRNA reader? From codon-hunt's ARAGORN screen of every R232 species rep (archive,
read-only: project-13/results/trna_screen/genomes.tsv, 191,958 genomes, per-genome anticodon counts).
Rules, fixed before the first run:
  usable genome: CheckM2 completeness >= 90, contamination <= 5, ARAGORN n_trna >= 30, and at least one other Arg tRNA
    (anticodon ACG, CCT or TCT) so that missing Arg tRNAs are not just missing tRNAs;
  no CGG reader: no tRNA with anticodon CCG or TCG (U34 reads CGA and CGG; ACG -> ICG reads CGU/CGC/CGA only).
usage: python -I src/analysis/cgg_reader_census.py <genomes.tsv> <gtdb.parquet> <outdir>
  -> <outdir>/cgg_reader_genomes.tsv (usable genomes with no CGG reader), cgg_reader_order.tsv, cgg_reader_family.tsv
"""
import csv
import sys
from collections import defaultdict
from pathlib import Path

import duckdb

OTHER_ARG = {"ACG", "CCT", "TCT"}


def main():
    gpath, tax, out = sys.argv[1], sys.argv[2], Path(sys.argv[3])
    meta = {a: r for a, *r in duckdb.sql(
        f"select substr(gtdb_acc, 4), phylum, class, \"order\", family, genus, checkm2_completeness, checkm2_contamination, "
        f"gc_percentage, ncbi_translation_table from '{tax}' where gtdb_representative = 't'").fetchall()}
    n = defaultdict(lambda: [0, 0])          # (level, taxon) -> [usable, no reader]
    rows, seen, skipped = [], 0, defaultdict(int)
    for r in csv.DictReader(open(gpath), delimiter="\t"):
        seen += 1
        m = meta.get(r["acc"])
        if m is None:
            skipped["no meta"] += 1
            continue
        phylum, cls, order, family, genus, comp, cont, gc, tt = m
        ac = {}
        for x in (r["anticodons"] or "").split(","):
            if x:
                k, v = x.split(":")
                ac[k.upper()] = int(v)
        if r["status"] != "ok" or comp is None or comp < 90 or cont > 5 or int(r["n_trna"]) < 30:
            skipped["quality"] += 1
            continue
        if not OTHER_ARG & ac.keys():
            skipped["no other Arg tRNA"] += 1
            continue
        none = "CCG" not in ac and "TCG" not in ac
        for lev, t in (("phylum", phylum), ("class", cls), ("order", order), ("family", family)):
            n[lev, t][0] += 1
            n[lev, t][1] += none
        if none:
            rows.append([r["acc"], phylum, cls, order, family, genus, f"{gc:.1f}", tt, r["n_trna"],
                         ",".join(f"{k}:{ac[k]}" for k in sorted(ac) if k in OTHER_ARG)])
    out.mkdir(parents=True, exist_ok=True)
    w = csv.writer(open(out / "cgg_reader_genomes.tsv", "w"), delimiter="\t", lineterminator="\n")
    w.writerow(["acc", "phylum", "class", "order", "family", "genus", "gc", "tt", "n_trna", "arg_anticodons"])
    w.writerows(rows)
    for lev in ("order", "family"):
        w = csv.writer(open(out / f"cgg_reader_{lev}.tsv", "w"), delimiter="\t", lineterminator="\n")
        w.writerow([lev, "usable", "no_reader", "share"])
        for (l, t), (u, k) in sorted(n.items(), key=lambda x: -x[1][1]):
            if l == lev and k:
                w.writerow([t, u, k, f"{k / u:.3f}"])
    usable = sum(u for (l, _), (u, _) in n.items() if l == "phylum")
    print(f"{seen} screened, {usable} usable, {len(rows)} with no CGG reader; skipped {dict(skipped)}")


if __name__ == "__main__":
    main()
