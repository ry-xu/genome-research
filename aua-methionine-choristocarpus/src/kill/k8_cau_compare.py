"""Kill check 4c: the CAU tRNAs of Choristocarpus and Discosporangium compared with trnW, Met-e, Met-f and Ile-CAU of
standard-code brown algae (global identity, anticodon loops). Discosporangium loci (unannotated) are cut by local
alignment of the closest annotated tRNA to the k7 locus +-40 nt.
usage: uv run python -I checks/chor_wgs/kill/k8_cau_compare.py GBDIR DMES_FA"""
import re, sys
from Bio import SeqIO
from Bio.Align import PairwiseAligner

GB, DMES = sys.argv[1:3]
SP = {"NC_007683.1": "Fucus", "NC_007684.1": "Desmarestia", "NC_085327.1": "Cutleria", "NC_013473.1": "Saccharina",
      "ON552993.1": "Phaeostrophion", "MZ156064.1": "Protohalopteris", "NC_085326.1": "Chor"}
PICK = {("Fucus", 3628): "MetE", ("Fucus", 8394): "IleCAU", ("Fucus", 34901): "MetF", ("Desmarestia", 4412): "MetE",
        ("Desmarestia", 9300): "IleCAU", ("Desmarestia", 37536): "MetF", ("Cutleria", 3287): "MetE", ("Cutleria", 8055): "IleCAU",
        ("Cutleria", 35861): "MetF", ("Saccharina", 3464): "MetE", ("Saccharina", 8353): "IleCAU", ("Saccharina", 36010): "MetF",
        ("Phaeostrophion", 17294): "MetE", ("Phaeostrophion", 22257): "IleCAU", ("Phaeostrophion", 46360): "MetF",
        ("Protohalopteris", 2762): "MetE", ("Protohalopteris", 7793): "IleCAU", ("Protohalopteris", 40915): "MetF",
        ("Chor", 2625): "Chor_2625", ("Chor", 2859): "Chor_2859"}
T = {}
for acc, sp in SP.items():
    r = SeqIO.read(f"{GB}/{acc}.gb", "genbank")
    for f in r.features:
        if f.type != "tRNA": continue
        st = int(f.location.start) + 1; s = str(f.extract(r.seq)).upper()
        lab = (f.qualifiers.get("gene") or f.qualifiers.get("product") or [""])[0]
        if (sp, st) in PICK: T[f"{sp}:{PICK[(sp, st)]}"] = s
        elif re.search(r"trnW|Trp", lab): T[f"{sp}:trnW"] = s
dm = "".join(l.strip().upper() for l in open(DMES) if not l.startswith(">"))
loc = PairwiseAligner(); loc.mode = "local"; loc.match_score, loc.mismatch_score = 2, -1; loc.open_gap_score, loc.extend_gap_score = -3, -1
def cut(query, lo, hi):
    w = dm[lo - 40:hi + 40]; a = loc.align(query, w)[0]
    (q0, q1), (s0, s1) = a.aligned[0][0], a.aligned[1][0]
    q0e, q1e = a.aligned[0][0][0], a.aligned[0][-1][1]; s0e, s1e = a.aligned[1][0][0], a.aligned[1][-1][1]
    # extend to full query length where the alignment stops short
    return w[max(0, s0e - q0e): s1e + (len(query) - q1e)]
T["Dmes:X_33026"] = cut(T["Chor:Chor_2859"], 33026, 33095)
T["Dmes:MetE_33179"] = cut(T["Cutleria:MetE"], 33179, 33245)
T["Dmes:MetF_16955"] = cut(T["Chor:Chor_2625"], 16955, 17028)
T["Dmes:trnW_21714"] = cut(T["Chor:trnW"], 21714, 21788)
g = PairwiseAligner(); g.mode = "global"; g.match_score, g.mismatch_score = 1, 0; g.open_gap_score, g.extend_gap_score = -1, -0.5
def ident(a, b):
    x = g.align(a, b)[0]; return sum(1 for i, j in zip(*x) if i == j and i != "-") / max(len(a), len(b))
def loopstr(s):
    c = [m.start() + 2 for m in re.finditer(r"(?=[ACGT]T(CAT|CCA)[AG])", s) if 30 <= m.start() + 2 <= 38]
    return s[c[0] - 2:c[0] + 5] if c else "?"
print("== sequences (len, anticodon loop)")
for k, s in T.items(): print(f"{k:28s} {len(s):3d} {loopstr(s)} {s}")
cls = ["MetE", "MetF", "IleCAU", "trnW"]
print("\n== mean / max global identity to standard-algae classes (Fucus, Desmarestia, Cutleria, Saccharina, Phaeostrophion, Protohalopteris)")
std = ["Fucus", "Desmarestia", "Cutleria", "Saccharina", "Phaeostrophion", "Protohalopteris"]
for k in ["Chor:Chor_2625", "Chor:Chor_2859", "Chor:trnW", "Dmes:X_33026", "Dmes:MetE_33179", "Dmes:MetF_16955", "Dmes:trnW_21714", "Fucus:IleCAU", "Fucus:MetE"]:
    out = []
    for c in cls:
        v = [ident(T[k], T[f"{sp}:{c}"]) for sp in std if f"{sp}:{c}" in T and f"{sp}:{c}" != k]
        out.append(f"{c} mean={sum(v) / len(v):.3f} max={max(v):.3f} (n={len(v)})")
    print(f"{k:22s} " + " | ".join(out))
print("\n== Chor_2859 vs Dmes X_33026:", round(ident(T["Chor:Chor_2859"], T["Dmes:X_33026"]), 3),
      "| Chor_2625 vs Dmes MetF:", round(ident(T["Chor:Chor_2625"], T["Dmes:MetF_16955"]), 3),
      "| Chor_2859 vs Dmes MetE:", round(ident(T["Chor:Chor_2859"], T["Dmes:MetE_33179"]), 3),
      "| Chor trnW vs Dmes trnW:", round(ident(T["Chor:trnW"], T["Dmes:trnW_21714"]), 3))
