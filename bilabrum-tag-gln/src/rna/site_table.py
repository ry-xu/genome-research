"""Join per-run read codons (RUN.sites.tsv from site_reads.py) with the reviewer's stop sites (bil.sites.tsv).
usage: site_table.py BIL_SITES.tsv RUN.sites.tsv [...]
Prints, per genome codon and site class (all aligned / clean), sites covered by >= 2 reads and their majority read codon."""
import sys
from collections import Counter, defaultdict


def rows(path):
    with open(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            yield dict(zip(head, line.rstrip("\n").split("\t")))


def main(bil, runs):
    site = {}
    for r in rows(bil):
        aligned = r["qres"] not in ("-", "None")
        clean = aligned and int(r["edge_dist"]) >= 10 and int(r["qgaps6"]) == 0 and int(r["win_ident"]) >= 8
        site[(r["contig"], r["pos1"], r["strand"])] = (r["codon"], aligned, clean, r["qres"])
    reads = defaultdict(Counter)
    for run in runs:
        for r in rows(run):
            k = (r["contig"], r["pos1"], r["strand"])
            reads[k]["TAG"] += int(r["TAG"])
            reads[k]["CAG"] += int(r["CAG"])
            for x in filter(None, r["other"].split(",")):
                c, n = x.split(":")
                if c == "indel":
                    continue
                reads[k][c] += int(n)
    for cls in ("aligned", "clean"):
        for cod in ("TAG", "TAA", "TGA"):
            ks = [k for k, v in site.items() if v[0] == cod and v[1 if cls == "aligned" else 2]]
            cov = [k for k in ks if sum(reads[k].values()) >= 2]
            maj = Counter(reads[k].most_common(1)[0][0] for k in cov)
            tot = sum((reads[k] for k in cov), Counter())
            print(f"{cls:7s} {cod}: sites {len(ks)}, covered (>=2 reads) {len(cov)}; majority read codon {dict(maj.most_common())}; "
                  f"reads {dict(tot.most_common(6))}")
    ks = [k for k, v in site.items() if v[0] == "TAG" and v[2] and sum(reads[k].values()) >= 2]
    for res in ("Q",):
        sub = [k for k in ks if site[k][3] == res]
        print(f"clean TAG opposite {res}: covered {len(sub)}; majority {dict(Counter(reads[k].most_common(1)[0][0] for k in sub))}")
    sub = [k for k in ks if site[k][3] != "Q"]
    print(f"clean TAG opposite other residues: covered {len(sub)}; majority {dict(Counter(reads[k].most_common(1)[0][0] for k in sub))}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2:])
