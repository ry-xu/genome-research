"""Arm B run list for the GTDB genomes chosen by src/armB_select.py (PREREG.md, Changes 2026-11-10 23:01).

Each genome is run with the table armB_select.tsv gives it (11 for a suspect table-4 call, else the Arm A calling
table). A genome already run in results/corio/ with that same table is not rerun; one run there with another table
is. Run order: flag, actino, strat, armc (a genome goes by its first reason in that order), split round-robin.
A GenBank assembly NCBI cannot serve is run from its RefSeq twin (GCA_x -> GCF_x, same assembly); labels.tsv maps
each Arm B label back to the GTDB accession.

usage: python -I src/armB_plan.py fetch          -> data/gtdb_armB/fetch.txt   accessions not on disk
       python -I src/armB_plan.py twins          -> data/gtdb_armB/fetch_twins.txt  GCF twins of unfetched GCA
       python -I src/armB_plan.py args <nparts>  -> data/gtdb_armB/args_<i>.txt, labels.tsv, missing.txt, reuse.tsv
"""
import collections
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SEL = ROOT / "results" / "gtdb" / "armB_select.tsv"
OUT = ROOT / "data" / "gtdb_armB"
CORIO = ROOT / "results" / "corio"
ORDER = ("flag", "actino", "strat", "armc")


def twin(acc):
    return "GCF_" + acc[4:] if acc.startswith("GCA_") else None


def path(acc, twins=True):
    for a in (acc, twin(acc)) if twins else (acc,):
        for p in (ROOT.parent / "tag-bacteria" / "data" / "genomes" / f"{a}.fna.gz",
                  ROOT / "data" / "genomes" / f"{a}.fna.gz"):
            if a and p.exists():
                return p
    return None


def selected():
    rows = list(csv.DictReader(open(SEL), delimiter="\t"))
    pri = lambda r: min(ORDER.index(x.split(":")[0]) for x in r["reasons"].split(","))
    return sorted(rows, key=lambda r: (pri(r), r["acc"]))


def corio_done():
    """label -> table used, over every Coriobacteriales Arm B run."""
    used = {}
    for fn in sorted(CORIO.glob("armB_*.tsv")):
        if fn.name.endswith("_events.tsv"):
            continue
        for r in csv.DictReader(open(fn), delimiter="\t"):
            used[r["label"]] = r["used"]
    return used


def main(mode, nparts=2):
    OUT.mkdir(parents=True, exist_ok=True)
    rows = selected()
    if mode == "fetch":
        need = [r["acc"] for r in rows if path(r["acc"]) is None]
        (OUT / "fetch.txt").write_text("".join(a + "\n" for a in need))
        print(f"{len(rows)} selected, {len(need)} to fetch")
        return
    if mode == "twins":
        need = [twin(r["acc"]) for r in rows if path(r["acc"]) is None and twin(r["acc"])]
        (OUT / "fetch_twins.txt").write_text("".join(a + "\n" for a in need))
        print(f"{len(need)} GCF twins to fetch")
        return
    done = corio_done()
    args, missing, reuse, labels = [], [], [], []
    for r in rows:
        prev = [a for a in (r["acc"], twin(r["acc"])) if a and done.get(a) == r["table"]]
        if prev:
            reuse.append((r["acc"], r["table"], f"results/corio:{prev[0]}"))
            continue
        p = path(r["acc"])
        if p is None:
            missing.append(r["acc"])
        else:
            args.append(f"{p}:{r['table']}")
            labels.append((p.name.split(".fna")[0], r["gtdb_acc"], r["table"], r["reasons"]))
    for i in range(nparts):
        (OUT / f"args_{i}.txt").write_text("".join(a + "\n" for a in args[i::nparts]))
    (OUT / "missing.txt").write_text("".join(a + "\n" for a in missing))
    with open(OUT / "labels.tsv", "w") as f:
        f.write("label\tgtdb_acc\ttable\treasons\n")
        f.writelines("\t".join(x) + "\n" for x in labels)
    with open(OUT / "reuse.tsv", "w") as f:
        f.write("acc\ttable\tfrom\n")
        f.writelines("\t".join(x) + "\n" for x in reuse)
    tabs = collections.Counter(a.rsplit(":", 1)[1] for a in args)
    print(f"{len(rows)} selected: {len(args)} to run (tables {dict(tabs)}), {len(reuse)} reused from Corio runs, "
          f"{len(missing)} missing")


if __name__ == "__main__":
    main(sys.argv[1], *(int(x) for x in sys.argv[2:]))
