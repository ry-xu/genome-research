"""Kill check 4: CAU tRNAs. Extract annotated tRNAs from GenBank flat files, read anticodon loops, classify each CAU tRNA
by best global identity to the labelled Fucus/Desmarestia trnM(cau)e (elongator), trnI(cau) (lysidine Ile) and trnM(cau)f
(initiator); search Choristocarpus (both sequences) and Discosporangium contig by blastn for CAU tRNA loci, annotated or not.
usage: uv run python -I checks/chor_wgs/kill/k6_trna.py GBDIR DMES_FA"""
import os, re, subprocess, sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner

GB, DMES = sys.argv[1:3]
tmp = "checks/chor_wgs/kill/tmp"; os.makedirs(tmp, exist_ok=True)
FILES = ["NC_007683.1", "NC_007684.1", "NC_085327.1", "NC_013473.1", "ON552993.1", "MZ156064.1", "NC_085326.1", "CAXNAD020000071.1"]
SP = {"NC_007683.1": "Fucus", "NC_007684.1": "Desmarestia", "NC_085327.1": "Cutleria", "NC_013473.1": "Saccharina",
      "ON552993.1": "Phaeostrophion", "MZ156064.1": "Protohalopteris", "NC_085326.1": "Chor_NC", "CAXNAD020000071.1": "Chor_WGS"}
al = PairwiseAligner(); al.mode = "global"; al.match_score, al.mismatch_score = 1, 0; al.open_gap_score, al.extend_gap_score = -1, -0.5

def loop(s):
    """anticodon candidates: 7-nt loop pattern NT NNN [AG] with the anticodon at 0-based 30..38"""
    c = [(m.start() + 2, s[m.start() + 2:m.start() + 5]) for m in re.finditer(r"(?=[ACGT]T[ACGT]{3}[AG])", s) if 30 <= m.start() + 2 <= 38]
    cat = [x for x in c if x[1] == "CAT"]
    return (cat or c or [(None, "?")])[0]

trna = []  # (species, acc, start, end, strand, label, seq, anticodon_pos, anticodon)
for acc in FILES:
    r = SeqIO.read(f"{GB}/{acc}.gb", "genbank")
    for f in r.features:
        if f.type != "tRNA": continue
        s = str(f.extract(r.seq)).upper()
        lab = (f.qualifiers.get("gene") or f.qualifiers.get("locus_tag") or f.qualifiers.get("product") or ["?"])[0]
        note = f.qualifiers.get("note", [""])[0]
        ac = f.qualifiers.get("anticodon", [""])[0]
        p, a = loop(s)
        trna.append((SP[acc], acc, int(f.location.start) + 1, int(f.location.end), f.location.strand, lab + (f" [{note}]" if note else ""), s, p, a, ac))
print("== annotated tRNAs with a CAT anticodon (loop-pattern read) or annotated Met/Ile")
cau = []
for t in trna:
    if (t[8] == "CAT" and not re.search(r"trnD|TRND", t[5])) or re.search(r"M|I-|Ile|Met|TRNI|trnI", t[5]):
        print(f"{t[0]:16s} {t[2]:>6}-{t[3]:<6} {'+' if t[4] == 1 else '-'} {t[5]:38s} len={len(t[6])} anticodon={t[8]} loop={t[6][t[7] - 2:t[7] + 5] if t[7] else '?'} {t[9]}")
        if t[8] == "CAT" and not re.search(r"trnD|TRND", t[5]): cau.append(t)
lab = {("Fucus", 3628): "Met-e", ("Fucus", 8394): "Ile-CAU", ("Fucus", 34901): "Met-f",
       ("Desmarestia", 4412): "Met-e", ("Desmarestia", 9300): "Ile-CAU", ("Desmarestia", 37536): "Met-f"}
refs = [(lab[(t[0], t[2])], t) for t in cau if (t[0], t[2]) in lab]
print("\n== each CAU tRNA: global identity to labelled Fucus/Desmarestia CAU tRNAs (best per class)")
def ident(a, b):
    x = al.align(a, b)[0]; m = sum(1 for i, j in zip(*x) if i == j and i != "-")
    return m / max(len(a), len(b))
for t in cau:
    best = {}
    for c, rt in refs:
        if rt is t: continue
        v = ident(t[6], rt[6]); best[c] = max(best.get(c, 0), v)
    call = max(best, key=best.get)
    print(f"{t[0]:16s} {t[2]:>6} {t[5]:38s} " + " ".join(f"{c}={best[c]:.3f}" for c in ["Met-e", "Ile-CAU", "Met-f"]) + f"  -> {call}")
# blastn search for CAU tRNA loci (annotated or not) in Choristocarpus and Discosporangium
qf = os.path.join(tmp, "cau_refs.fa")
with open(qf, "w") as fh:
    for c, t in refs: fh.write(f">{t[0]}_{c}\n{t[6]}\n")
    for t in cau:
        if t[0] in ("Chor_NC", "Cutleria", "Phaeostrophion", "Protohalopteris", "Saccharina"): fh.write(f">{t[0]}_{t[2]}\n{t[6]}\n")
RD = "scratchpad/chor_work"
for name, fa in [("Chor_NC", f"{RD}/NC_085326.fa"), ("Chor_WGS", f"{RD}/cten_mito.fa"), ("Dmes", DMES)]:
    out = subprocess.run(["blastn", "-task", "blastn", "-word_size", "7", "-evalue", "1e-5", "-dust", "no", "-query", qf, "-subject", fa,
                          "-outfmt", "6 qseqid sstart send pident length evalue bitscore"], check=True, capture_output=True, text=True).stdout
    loci = {}
    for l in out.splitlines():
        q, ss, se, pid, L, ev, bs = l.split("\t"); ss, se = int(ss), int(se)
        k = (min(ss, se) // 100)
        key = next((x for x in loci if abs(x[0] - min(ss, se)) < 40), (min(ss, se), max(ss, se), "+" if se > ss else "-"))
        loci.setdefault(key, []).append((float(bs), q, pid, L))
    g = "".join(l.strip().upper() for l in open(fa) if not l.startswith(">"))
    print(f"\n== {name}: CAU-tRNA blastn loci (merged within 40 nt); best 3 queries; anticodon read from extracted locus")
    for (a, b, st), hits in sorted(loci.items()):
        hits.sort(reverse=True)
        s = g[max(0, a - 6):b + 6]
        if st == "-": s = str(Seq(s).reverse_complement())
        cats = [m.start() for m in re.finditer(r"(?=[ACGT]TCAT[AG])", s)]
        print(f"  {a}-{b} {st} anticodon-loop CAT at {cats} | " + "; ".join(f"{q} {pid}%/{L} bits={bs:.0f}" for bs, q, pid, L in hits[:3]))
