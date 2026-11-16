"""Positive control 2: plant frameshifts into a genome and check that the same pipeline finds them.

usage: 07_plant.py SRC_LABEL N  -> genomes/<SRC>planted.fna, out/<SRC>planted.sites.tsv
Site types (N each), at aligned codons with window identity >= 0.5, >= 30 columns from either alignment end,
sites >= 600 nt apart:
  cgg_gly : at a Gly codon GGN insert a C in front of it, so the genome reads ..CGG N.. with the CGG in the old frame
            and a +1 shift after it -- exactly the pattern tested at natural CGG; the expected residue is Gly
  plus1   : insert one random nt after a random codon (needs a +1 shift)
  minus1  : delete the first nt after a random codon (needs a -1 shift)
Then run miniprot + 03_analyze.py + 04_fs_detail.py on <SRC>planted and score with 08_score_plant.py.
"""
import sys, random
sys.path.insert(0, "src/recheck")
from mp import read_fasta, parse_paf, best_per_ref, STD

W = "recheck_work"
src, N = sys.argv[1], int(sys.argv[2])
TYPES = sys.argv[3].split(",") if len(sys.argv) > 3 else ["cgg_gly", "plus1", "minus1"]
OUT = sys.argv[4] if len(sys.argv) > 4 else f"{src}planted"
random.seed(20261111)
genome = read_fasta(f"{W}/genomes/{src}.fna")
prots = read_fasta(f"{W}/refs/refs_all.faa")
alns = [a for a in parse_paf(f"{W}/out/{src}.paf", genome, prots) if a["nM"] >= 30 and a["AS"] >= 60]
best = best_per_ref(alns)
cand = []
for (ctg, st, s, ref), (AS, ai, o, codon, raa, wid, neargap, kf, kt) in best.items():
    if wid >= 0.5 and not neargap and kf >= 30 and kt >= 30 and not alns[ai]["events"]:
        cand.append((ctg, st, s, codon, raa))
random.shuffle(cand)
chosen, used = [], {}
want = {t: (N if t in TYPES else 0) for t in ["cgg_gly", "cgg_arg", "plus1", "minus1"]}
for (ctg, st, s, codon, raa) in cand:
    if all(v == 0 for v in want.values()): break
    if any(abs(s - x) < 600 for x in used.get(ctg, [])): continue
    if want["cgg_gly"] and codon.startswith("GG") and raa == "G":
        t = "cgg_gly"
    elif want["cgg_arg"] and STD.get(codon) == "R" and raa == "R":
        t = "cgg_arg"
    elif want["plus1"] and want["plus1"] >= want["minus1"]:
        t = "plus1"
    elif want["minus1"]:
        t = "minus1"
    else:
        continue
    want[t] -= 1
    used.setdefault(ctg, []).append(s)
    chosen.append((ctg, st, s, codon, raa, t))
COMPL = {"A": "T", "C": "G", "G": "C", "T": "A"}
edits = {}  # ctg -> list of (pos, kind, base): kind ins (insert base before pos) / del (delete base at pos)
rows = []
for (ctg, st, s, codon, raa, t) in chosen:
    b = random.choice("ACGT")
    if st == "+":
        if t == "cgg_arg": edits.setdefault(ctg, []).append((s, "sub3", "CGG"))
        e = {"cgg_gly": (s, "ins", "C"), "cgg_arg": (s + 3, "ins", b), "plus1": (s + 3, "ins", b), "minus1": (s + 3, "del", "")}[t]
    else:
        if t == "cgg_arg": edits.setdefault(ctg, []).append((s, "sub3", "CCG"))
        e = {"cgg_gly": (s + 3, "ins", "G"), "cgg_arg": (s, "ins", COMPL[b]), "plus1": (s, "ins", COMPL[b]), "minus1": (s - 1, "del", "")}[t]
    edits.setdefault(ctg, []).append(e)
    rows.append((ctg, st, s, codon, raa, t, e[0]))
newpos = {}
out = open(f"{W}/genomes/{OUT}.fna", "w")
for ctg, seq in genome.items():
    es = sorted(edits.get(ctg, []))
    parts, last, shift, shifts = [], 0, 0, []
    for (p, kind, base) in es:
        parts.append(seq[last:p])
        shifts.append((p, shift))
        if kind == "ins":
            parts.append(base); last = p; shift += 1
        elif kind == "sub3":
            parts.append(base); last = p + 3
        else:
            last = p + 1; shift -= 1
    parts.append(seq[last:])
    out.write(f">{ctg}\n{''.join(parts)}\n")
    for (p, sh) in shifts:
        newpos[(ctg, p)] = p + sh
out.close()
with open(f"{W}/out/{OUT}.sites.tsv", "w") as f:
    f.write("contig\tstrand\torig_codon_start0\tcodon\tref_res\ttype\tedit_pos_orig\tedit_pos_new\n")
    for (ctg, st, s, codon, raa, t, p) in rows:
        f.write(f"{ctg}\t{st}\t{s}\t{codon}\t{raa}\t{t}\t{p}\t{newpos[(ctg, p)]}\n")
from collections import Counter
print(src, "planted", Counter(r[5] for r in rows))
