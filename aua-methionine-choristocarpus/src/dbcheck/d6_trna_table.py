"""Tabulate tRNA-finder results: all tRNAs with anticodon CAT/TAT/AAT/GAT (the ones that could pair with AUA or are Ile),
per genome and tool; isotype scores (tRNAscan -B); identity of every CAT tRNA to the curated Desmarestia
trnI(cau), trnM(cau)e, trnM(cau)f (NC_007684.1 GenBank annotation) and to Saccharina's three CAT tRNAs."""
import re, glob, os
from Bio import SeqIO
from Bio.Align import PairwiseAligner
B = "lab/work/organelle-reading/checks/chor_wgs/dbcheck/"
G = ["NC_085326", "CAXNAD020000071.1", "CAXNAE020000396.1", "NC_013473.1", "ON552993.1", "NC_007684.1"]
NAME = {"NC_085326": "Choristocarpus (RefSeq)", "CAXNAD020000071.1": "Choristocarpus (Genoscope WGS)",
        "CAXNAE020000396.1": "Discosporangium (WGS contig)", "NC_013473.1": "Saccharina", "ON552993.1": "Phaeostrophion",
        "NC_007684.1": "Desmarestia"}
WATCH = {"CAT", "TAT", "AAT", "GAT"}

def parse_tscan(path):
    rows = []
    for line in open(path):
        p = [x.strip() for x in line.rstrip("\n").split("\t")]
        if len(p) < 9 or not p[1].isdigit():
            continue
        rows.append(p)
    return rows

def parse_ss(path):
    out = {}
    txt = open(path).read()
    for blk in txt.strip().split("\n\n"):
        m = re.match(r"(\S+)\.trna(\d+) \((\d+)-(\d+)\)", blk)
        s = re.search(r"^Seq: (\S+)", blk, re.M)
        if m and s:
            out[int(m.group(2))] = (int(m.group(3)), int(m.group(4)), s.group(1).upper().replace("U", "T"))
    return out

def parse_iso(path):
    out = {}
    lines = open(path).read().strip().split("\n")
    hdr = lines[0].split("\t")
    for l in lines[1:]:
        p = l.split("\t"); out[int(p[0].rsplit("trna", 1)[1])] = dict(zip(hdr, p))
    return out

def parse_aragorn(path):
    rows = []
    for line in open(path):
        m = re.match(r"\s*\d+\s+tRNA-(\S+)\s+(c?)\[(\d+),(\d+)\]\s+\S+\s+\((\w+)\)", line)
        if m:
            rows.append((m.group(1), m.group(2) == "c", int(m.group(3)), int(m.group(4)), m.group(5).upper()))
    return rows

print("== Counts of predicted tRNAs per tool")
for g in G:
    o = parse_tscan(B + f"trna/{g}.tscanO.txt"); gg = parse_tscan(B + f"trna/{g}.tscanG.txt")
    bb = parse_tscan(B + f"trna/{g}.tscanB.txt"); a = parse_aragorn(B + f"trna/{g}.aragorn.txt"); aps = parse_aragorn(B + f"trna/{g}.aragorn_ps.txt")
    print(f"{NAME[g]:32s} tRNAscan-O={len(o)} tRNAscan-G={len(gg)} tRNAscan-B={len(bb)} ARAGORN={len(a)} ARAGORN-ps={len(aps)}")

print("\n== tRNAs with anticodon CAT/TAT/AAT/GAT, all tools (begin, end, type, anticodon, score, note)")
for g in G:
    print(f"-- {NAME[g]} ({g})")
    for tool in ["tscanO", "tscanG", "tscanB"]:
        for p in parse_tscan(B + f"trna/{g}.{tool}.txt"):
            if p[5] in WATCH:
                print(f"   {tool:7s} #{p[1]:>2} {p[2]:>6}-{p[3]:<6} {p[4]:5s} {p[5]}  score={p[8]}  " + "  ".join(p[9:]))
    for tool in ["aragorn", "aragorn_ps"]:
        for t in parse_aragorn(B + f"trna/{g}.{tool}.txt"):
            if t[4] in WATCH:
                print(f"   {tool:10s} {t[2]:>6}-{t[3]:<6} {t[0]:5s} {t[4]} {'minus' if t[1] else 'plus'}")

# reference CAT tRNAs from curated Desmarestia annotation
des = SeqIO.read(B + "data/NC_007684.1.gb", "genbank")
REF = {}
for f in des.features:
    if f.type == "tRNA":
        gname = f.qualifiers.get("gene", [""])[0]
        if gname in ("trnI(cau)", "trnM(cau)e", "trnM(cau)f"):
            REF["Desm " + gname] = str(f.extract(des.seq)).upper()
sac = SeqIO.read(B + "data/NC_013473.1.gb", "genbank")
for f in sac.features:
    if f.type == "tRNA":
        lt = f.qualifiers.get("locus_tag", [""])[0]
        if lt in ("SaanoM_t04", "SaanoM_t11", "SaanoM_t24"):
            REF[f"Sacc {lt} {f.qualifiers.get('product',[''])[0]}"] = str(f.extract(sac.seq)).upper()
print("\n== Reference CAT tRNAs:", {k: len(v) for k, v in REF.items()})

al = PairwiseAligner(); al.mode = "global"; al.match_score = 1; al.mismatch_score = -1
al.open_gap_score = -2; al.extend_gap_score = -1; al.end_gap_score = 0
def ident(a, b):
    aln = al.align(a, b)[0]
    m = cols = 0
    for (s1, e1), (s2, e2) in zip(*aln.aligned):
        for i, j in zip(range(s1, e1), range(s2, e2)):
            cols += 1; m += a[i] == b[j]
    return m, cols

print("\n== Every CAT tRNA (tRNAscan -B) vs reference CAT tRNAs: identity = matches/aligned (ungapped-pair) columns")
print("   isotype scores from tRNAscan -B isotype-specific CMs: Met, fMet, Ile, Ile2 and best isotype")
for g in G:
    ss = parse_ss(B + f"trna/{g}.tscanB.ss"); iso = parse_iso(B + f"trna/{g}.tscanB.iso.txt")
    for p in parse_tscan(B + f"trna/{g}.tscanB.txt"):
        if p[5] != "CAT":
            continue
        n = int(p[1]); b, e, s = ss[n]; d = iso[n]
        scores = {k: float(v) for k, v in d.items() if k not in ("tRNAscanID", "Anticodon_predicted_isotype")}
        best = max(scores, key=scores.get)
        ids = []
        for rk, rs in REF.items():
            m, c = ident(s, rs); ids.append((m / c, rk, m, c))
        ids.sort(reverse=True)
        print(f"-- {NAME[g]:30s} {b}-{e} type={p[4]} Met={scores['Met']} fMet={scores['fMet']} Ile={scores['Ile']} "
              f"Ile2={scores['Ile2']} best={best}({scores[best]})")
        for fr, rk, m, c in ids:
            print(f"      {fr:.3f} ({m}/{c})  {rk}")
