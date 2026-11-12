"""Calibration: pyskani ANI among the 10 clade references and from 13 non-clade Polyangiaceae genomes on disk.
usage: uv run --no-project --with pyskani python -I src/catalogues/skani_calib.py
"""
import gzip, time
from pathlib import Path
import pyskani
ROOT = Path(__file__).resolve().parents[2]
G = ROOT / "data/genomes"
REFS = {"GCA_020849655.1": "JADLFU01", "GCA_963169165.1": "JADLFU01", "GCA_963169955.1": "JADLFU01",
        "GCA_964605685.1": "JANWXX01", "GCA_025453595.1": "JANWXX01", "GCA_025057345.1": "JANWXX01",
        "GCA_033578855.1": "JANWXX01", "GCA_035762175.1": "DASTPP01", "GCA_046351315.1": "JAXLYQ01",
        "GCA_964396175.1": "CAYNEK01"}
NEG = [l.split("\t") for l in open(ROOT / "results/catalogues/neg_controls.tsv").read().split("\n")[1:] if l]
def seqs(p):
    out, cur = [], []
    with gzip.open(p, "rt") as f:
        for line in f:
            if line.startswith(">"):
                if cur: out.append("".join(cur).encode()); cur = []
            else: cur.append(line.strip())
    if cur: out.append("".join(cur).encode())
    return out
t0 = time.time()
db = pyskani.Database()
S = {}
for a in REFS:
    S[a] = seqs(G / f"{a}.fna.gz"); db.sketch(a, *S[a])
print("sketch s", round(time.time() - t0, 1))
print("query\tq_genus\tref\tref_genus\tani\taf_q\taf_r")
for a, g in REFS.items():
    for h in db.query(a, *S[a]):
        if h.reference_name != a:
            print(a, g, h.reference_name, REFS[h.reference_name], round(h.identity * 100, 2), round(h.query_fraction * 100, 1), round(h.reference_fraction * 100, 1), sep="\t")
for a, g in NEG:
    t1 = time.time()
    hs = db.query(a, *seqs(G / f"{a}.fna.gz"))
    print(a, g, "NEG", len(hs), [(h.reference_name, round(h.identity * 100, 2), round(h.query_fraction * 100, 1)) for h in hs], round(time.time() - t1, 1), sep="\t")
