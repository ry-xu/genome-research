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
al = PairwiseAligner(); al.mode = "global"; al.match_score = 1; al.mismatch_score = -1
al.open_gap_score = -2; al.extend_gap_score = -1; al.end_gap_score = 0
