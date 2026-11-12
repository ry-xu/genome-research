"""Check 1 (clade) of PREREG.md "From an Arm B pass to a result", as made operational in Changes (2026-11-11 06:21).

For each candidate whose genus has no other assembly passing X in GTDB Arm B: every other GTDB R232 assembly of its
species, then the genus's not-yet-run assemblies (ranked by CheckM2 completeness - 5 x contamination, ties by
accession) up to 20 added genomes per genus.

usage: uv run --no-project python -I src/check1_plan.py  -> data/check1/plan.tsv, data/check1/fetch.txt,
                                                           data/check1/args.txt (<fasta>:g per genome)
"""
import csv
import gzip
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
G = ROOT / "results" / "gtdb"
OUT = ROOT / "data" / "check1"
GEN = ROOT / "data" / "genomes"
CENSUS = ROOT.parent / "tag-bacteria" / "data" / "genomes"
CAP = 20


def tsv(fn):
    return list(csv.DictReader(open(fn), delimiter="\t"))


def path(acc):
    for d in (GEN, CENSUS):
        p = d / f"{acc}.fna.gz"
        if p.exists():
            return p
    return None


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    cand = tsv(G / "candidates.tsv")
    passing = {(r["gtdb_acc"], r["codon"]) for r in tsv(G / "armB_passes.tsv")}
    run = {r["gtdb_acc"] for r in tsv(G / "armB_select.tsv")}
    need = {}
    for c in cand:
        if int(c["genus_others_pass"]) == 0:
            need[c["genus"]] = c
    meta = {}
    for fn in ("ar53_metadata.tsv.gz", "bac120_metadata.tsv.gz"):
        with gzip.open(ROOT / "data" / "gtdb_meta" / fn, "rt") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                tax = r["gtdb_taxonomy"].split(";")
                if tax[5] in need:
                    meta[r["accession"]] = r
    rows = []
    for g, c in sorted(need.items()):
        sp = c["taxonomy"].split(";")[6]
        pool = [a for a, r in meta.items() if r["gtdb_taxonomy"].split(";")[5] == g and a != c["gtdb_acc"]
                and a not in run]
        score = lambda a: (-(float(meta[a]["checkm2_completeness"]) - 5 * float(meta[a]["checkm2_contamination"])), a)
        same = sorted((a for a in pool if meta[a]["gtdb_taxonomy"].split(";")[6] == sp), key=score)
        rest = sorted((a for a in pool if a not in same), key=score)
        chosen = same + rest[:max(0, CAP - len(same))]
        print(f"{g} {c['codon']} ({c['gtdb_acc']}, {sp}): pool {len(pool)}, same species {len(same)}, "
              f"chosen {len(chosen)}")
        for a in chosen:
            m = meta[a]
            rows.append(dict(gtdb_acc=a, acc=a[3:], genus=g, codon=c["codon"], candidate=c["gtdb_acc"],
                             same_species=int(a in same), species=m["gtdb_taxonomy"].split(";")[6],
                             completeness=m["checkm2_completeness"], contamination=m["checkm2_contamination"],
                             representative=m["gtdb_representative"]))
    with open(OUT / "plan.tsv", "w") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    missing = [r["acc"] for r in rows if path(r["acc"]) is None]
    (OUT / "fetch.txt").write_text("".join(a + "\n" for a in missing))
    (OUT / "args.txt").write_text("".join(f"{path(r['acc']) or GEN / (r['acc'] + '.fna.gz')}:g\n" for r in rows))
    print(f"{len(rows)} genomes, {len(missing)} to fetch -> {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
