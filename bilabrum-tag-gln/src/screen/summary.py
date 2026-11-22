"""Summarise a full call table: coverage, wider recall on known lineages outside the controls, and every lead.

usage: summary.py CALLS.tsv UNITS.tsv LINEAGE.tsv KNOWN.tsv CONTROLS.tsv > report.txt
CALLS.tsv is call.py output over all units. A unit's lineage key is duplicate_of if set, else the unit itself.
Wider recall: every non-control unit whose lineage holds a known row. The expectation is the deepest known row in
its lineage (the genus, the first word of the name, is deepest), so a peritrich expects E, not the class's Q. Clade
rows that the ledger marks as clade statements (Colpodea, Nassophorea) are reported but are weak expectations.
Leads: every call/watch whose label is not known (label.py's rules), with its numbers.
"""
import csv
import sys
from collections import Counter, defaultdict


def main(calls, units_tsv, lineage_tsv, known_tsv, controls_tsv):
    units = {r["unit"]: r for r in csv.DictReader(open(units_tsv), delimiter="\t")}
    lin = {r["taxid"]: r["lineage"].split("; ") for r in csv.DictReader(open(lineage_tsv), delimiter="\t")}
    known = {r["clade"]: (r["status"], set(filter(None, r["pairs"].split(","))))
             for r in csv.DictReader(open(known_tsv), delimiter="\t")}
    ctl = {r["unit"] for r in csv.DictReader(open(controls_tsv), delimiter="\t")}
    rows = defaultdict(dict)
    for r in csv.DictReader(open(calls), delimiter="\t"):
        rows[r["unit"]][r["codon"]] = r

    def path(u):
        x = units[u]
        return lin.get(x["taxid"], []) + [x["name"].split(" ")[0]]

    def label(u, pair):
        hit = [c for c in path(u) if c in known]
        st = {known[c][0] for c in hit if pair in known[c][1]}
        return "known" if "known" in st else "candidate" if st else "lineage" if hit else "new"

    print("## Coverage (units with a score row; powered = >= 100,000 inner columns)")
    cov = defaultdict(Counter)
    for u, x in units.items():
        s = x["stratum"]
        cov[s]["units"] += 1
        if u in rows:
            cov[s]["scored"] += 1
            cov[s]["powered"] += next(iter(rows[u].values()))["powered"] == "yes"
    print("stratum\tunits\tscored\tpowered")
    for s in sorted(cov):
        print(f"{s}\t{cov[s]['units']}\t{cov[s]['scored']}\t{cov[s]['powered']}")
    print(f"all\t{sum(c['units'] for c in cov.values())}\t{sum(c['scored'] for c in cov.values())}\t"
          f"{sum(c['powered'] for c in cov.values())}")

    print("\n## Calls and watches by label (codon rows; distinct lineages in brackets)")
    tally = defaultdict(set)
    tc = Counter()
    for u, cc in rows.items():
        for codon, r in cc.items():
            if r["tier"] in ("call", "watch"):
                k = (r["tier"], label(u, f"{codon}:{r['meaning']}"))
                tc[k] += 1
                tally[k].add(units[u]["duplicate_of"] or u)
    for k in sorted(tc):
        print(f"{k[0]}\t{k[1]}\t{tc[k]}\t[{len(tally[k])}]")

    print("\n## Wider recall: non-control units in a lineage with a known row")
    print("unit\tstratum\tname\trow\texpected\tpowered\tinner\tresult")
    rec = defaultdict(Counter)
    for u in sorted(rows):
        if u in ctl:
            continue
        p = path(u)
        deep = [c for c in p if c in known and known[c][0] == "known"]
        if not deep:
            continue
        row = deep[-1]
        exp = defaultdict(set)
        for pr in known[row][1]:
            c, m = pr.split(":")
            exp[c].add(m)
        res = []
        any_r = next(iter(rows[u].values()))
        for c in sorted(exp):
            r = rows[u].get(c)
            if r and r["tier"] in ("call", "watch") and r["meaning"] in exp[c]:
                out = r["tier"]
            elif r and r["tier"] in ("call", "watch"):
                out = f"{r['tier']}-other:{r['meaning']}"
            else:
                out = "miss"
            res.append(f"{c}:{'/'.join(sorted(exp[c]))}={out}(n={r['n'] if r else 0})")
            rec[(row, any_r["powered"])][out.split("-")[0] if out != "miss" else "miss"] += 1
        x = units[u]
        print(f"{u}\t{x['stratum']}\t{x['name']}\t{row}\t{','.join(sorted(known[row][1]))}\t{any_r['powered']}\t"
              f"{any_r['inner']}\t{' '.join(res)}")
    print("\nrow\tpowered\tcall\twatch\tmiss")
    for (row, pw), c in sorted(rec.items()):
        print(f"{row}\t{pw}\t{c['call']}\t{c['watch']}\t{c['miss']}")
    tot = Counter()
    for (row, pw), c in rec.items():
        if pw == "yes":
            tot.update(c)
    print(f"ALL powered\t\t{tot['call']}\t{tot['watch']}\t{tot['miss']}")

    print("\n## Leads: calls and watches not labelled known")
    keys = ["n", "prot", "prot_meaning", "rate_al", "rel", "inner", "powered"]
    print("\t".join(["unit", "stratum", "name", "duplicate_of", "codon", "meaning", "tier", "label"] + keys +
                    ["lineage"]))
    for u in sorted(rows):
        for codon, r in sorted(rows[u].items()):
            if r["tier"] not in ("call", "watch"):
                continue
            lb = label(u, f"{codon}:{r['meaning']}")
            if lb == "known":
                continue
            x = units[u]
            print("\t".join([u, x["stratum"], x["name"], x["duplicate_of"], codon, r["meaning"], r["tier"], lb] +
                            [r[k] for k in keys] + ["; ".join(path(u)[-6:-1])]))


if __name__ == "__main__":
    main(*sys.argv[1:6])
