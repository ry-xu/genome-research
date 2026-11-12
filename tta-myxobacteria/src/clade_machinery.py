"""Clade view of a candidate codon: every GTDB species representative in the given families, with Arm A (row for X
if lr >= 0.5, flag), Arm B (X's O, O_fam, ratio, log10p, events; pass) where run, RF motifs (rf_motifs.tsv) and tRNA
readers of X from the seed-archive ARAGORN screen (project-13 results/trna_screen/genomes.tsv, read-only; covers
191,958 of 199,923 representatives).

Readers counted per codon (DNA anticodons, standard wobble; archaea have no A34 inosine):
  TTA: TAA (Leu)          CGA: TCG, and ACG in bacteria (I34)          others: the exact complement only
usage: uv run --no-project python -I src/clade_machinery.py <codon> <family> [<family> ...]
       -> results/gtdb/clade_<codon>_<first family>.tsv
"""
import collections
import csv
import gzip
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
G = ROOT / "results" / "gtdb"
TRNA = Path("data/trna_screen/genomes.tsv")
COMP = str.maketrans("ACGT", "TGCA")


def readers(x, domain):
    if x == "CGA":
        return ["TCG"] + (["ACG"] if domain == "bacteria" else [])
    return [x.translate(COMP)[::-1]]


def tsv(fn):
    return csv.DictReader(open(fn), delimiter="\t")


def main(x, fams):
    gen = {r["acc"]: r for r in tsv(G / "armA_genomes.tsv") if any(f"f__{f};" in r["taxonomy"] for f in fams)}
    arma = {}
    with gzip.open(G / "armA_rows.tsv.gz", "rt") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["codon"] == x and r["acc"] in gen:
                arma[r["acc"]] = r
    label_of = {r["gtdb_acc"]: r["label"] for r in tsv(ROOT / "data" / "gtdb_armB" / "labels.tsv")}
    want = {label_of[a]: a for a in gen if a in label_of}
    armb = {}
    for i in (0, 1):
        for r in tsv(G / f"armB_{i}.tsv"):
            if r["codon"] == x and r["label"] in want:
                armb[want[r["label"]]] = r
    rf = collections.defaultdict(list)
    for r in tsv(G / "rf_motifs.tsv"):
        if r["acc"] in gen and r["class"] in ("RF1", "RF2", "aRF1"):
            rf[r["acc"]].append(f"{r['class']}:{r['residues']}" + ("" if r["canonical"] == "1" else "*"))
    trna = {}
    for r in tsv(TRNA):
        trna[r["acc"]] = r
    out = []
    for a, g in sorted(gen.items(), key=lambda kv: kv[1]["taxonomy"]):
        t = trna.get(a[3:])
        ac = dict(kv.split(":") for kv in t["anticodons"].split(",") if kv) if t else None
        rd = "na" if ac is None else ",".join(f"{k}:{ac.get(k, 0)}" for k in readers(x, g["domain"]))
        A, B = arma.get(a, {}), armb.get(a, {})
        bpass = bool(B) and int(B["O_fam"]) >= 10 and float(B["ratio"]) >= 3 and float(B["log10p"]) <= -10
        tax = g["taxonomy"].split(";")
        out.append(dict(acc=a, family=tax[4], genus=tax[5], species=tax[6], gc=g["gc"], ncomplete=g["ncomplete"],
                        armA_flags=g["flags"], armA_lr=A.get("lr", "<0.5"), armA_n3=A.get("n3", ""),
                        armA_nI=A.get("nI", ""), armB_events=B.get("events", ""), armB_O=B.get("O", ""),
                        armB_Ofam=B.get("O_fam", ""), armB_ratio=B.get("ratio", ""), armB_log10p=B.get("log10p", ""),
                        armB_pass=("yes" if bpass else "no") if B else "", rf=" ".join(sorted(rf[a])),
                        trna_n=t["n_trna"] if t else "na", readers=rd))
    fn = G / f"clade_{x}_{fams[0]}.tsv"
    with open(fn, "w") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0]), delimiter="\t")
        w.writeheader()
        w.writerows(out)
    print(f"{len(out)} representatives in {', '.join(fams)} -> {fn.relative_to(ROOT)}")
    by = collections.defaultdict(list)
    for r in out:
        by[r["genus"]].append(r)
    for gname, rs in by.items():
        flag = sum(x in r["armA_flags"] for r in rs)
        run = [r for r in rs if r["armB_pass"]]
        npass = sum(r["armB_pass"] == "yes" for r in rs)
        rfs = collections.Counter(m for r in rs for m in r["rf"].split() if m.startswith(("RF2", "aRF1:NIKS", "aRF1")))
        rd = collections.Counter(r["readers"] for r in rs)
        print(f"{rs[0]['family']:28s} {gname:24s} n={len(rs):3d} GC~{sorted(float(r['gc']) for r in rs)[len(rs)//2]:.2f}"
              f"  armA {x} flag {flag:3d}  armB pass {npass}/{len(run)}  RF {dict(rfs.most_common(4))}"
              f"  readers {dict(rd.most_common(3))}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2:])
