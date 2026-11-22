"""Shared helpers: read cmalign --matchonly Stockholm, map columns to E. coli numbering via census columns.tsv."""
import csv

COLS = {"23S": "lab/work/frozen-core/results/lsu/bac/columns.tsv",
        "16S": "lab/work/frozen-core/results/ssu/bac/columns.tsv"}


def read_sto(fn):
    seqs, order = {}, []
    for line in open(fn):
        if not line.strip() or line.startswith("#") or line.startswith("//"):
            continue
        name, s = line.split()
        if name not in seqs:
            order.append(name); seqs[name] = []
        seqs[name].append(s)
    out = {}
    for n in order:
        s = "".join(seqs[n]).upper().replace("T", "U").replace(".", "-")
        out[n] = s
    return out


def read_cols(kind):
    rows = list(csv.DictReader(open(COLS[kind]), delimiter="\t"))
    for r in rows:
        r["col"] = int(r["col"]); r["ecoli_pos"] = int(r["ecoli_pos"]); r["f_cons"] = float(r["f_cons"])
        r["f_gap"] = float(r["f_gap"]); r["partner"] = int(r["partner"])
    return rows


def short(name):
    """23S_004|GCA_002782845.1|Aquicultor_A|div -> 23S_004"""
    return name.split("|")[0]
