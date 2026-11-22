"""Label every call and watch from call.py against the known-codes table and NCBI's genetic code.

usage: label.py CALLS.tsv UNITS.tsv LINEAGE.tsv KNOWN.tsv > labelled.tsv
KNOWN.tsv (data/known_codes.tsv, from notes/known-nuclear-stop-codes.md) has clade, status, pairs, basis. A clade is an
NCBI lineage name or a genus (the first word of the unit's name); status is known, candidate (published but unresolved)
or group (a clade holding a known code, with no pairs of its own).
label  known      a known row in the unit's lineage lists this codon:meaning
       candidate  only a candidate row lists it
       lineage    the lineage holds a row (any status), but none lists this codon:meaning
       new        no row in the lineage
ncbi   yes if NCBI's genetic code for the taxid translates the codon as the meaning, no otherwise (NCBI's codes are often
       wrong for ciliates; the label, not this column, is the answer)
Only call/watch rows are printed, with the clades that matched and the lineage's last six ranks.
"""
import csv
import sys

from Bio.Data import CodonTable


def main(calls, units_tsv, lineage_tsv, known_tsv):
    units = {r["unit"]: r for r in csv.DictReader(open(units_tsv), delimiter="\t")}
    lin = {r["taxid"]: r["lineage"].split("; ") for r in csv.DictReader(open(lineage_tsv), delimiter="\t")}
    known = {r["clade"]: (r["status"], set(filter(None, r["pairs"].split(","))))
             for r in csv.DictReader(open(known_tsv), delimiter="\t")}
    cols = ["unit", "stratum", "name", "group", "taxid", "gc", "duplicate_of", "codon", "meaning", "tier", "n", "prot",
            "prot_meaning", "rate_al", "rel", "inner", "powered", "label", "clade", "ncbi", "lineage"]
    print("\t".join(cols))
    for r in csv.DictReader(open(calls), delimiter="\t"):
        if r["tier"] not in ("call", "watch"):
            continue
        u = units[r["unit"]]
        L = lin.get(u["taxid"], [])
        names = set(L) | {u["name"].split(" ")[0]}
        hit = [c for c in known if c in names]
        pair = f"{r['codon']}:{r['meaning']}"
        st = {known[c][0] for c in hit if pair in known[c][1]}
        label = "known" if "known" in st else "candidate" if st else "lineage" if hit else "new"
        ncbi = ""
        if u["gc"]:
            t = CodonTable.unambiguous_dna_by_id.get(int(u["gc"]))
            ncbi = "yes" if t and t.forward_table.get(r["codon"]) == r["meaning"] else "no"
        r.update(taxid=u["taxid"], gc=u["gc"], duplicate_of=u["duplicate_of"], label=label, clade=",".join(hit),
                 ncbi=ncbi, lineage="; ".join(L[-6:]))
        print("\t".join(str(r[c]) for c in cols))


if __name__ == "__main__":
    main(*sys.argv[1:5])
