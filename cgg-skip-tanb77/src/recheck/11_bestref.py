"""Identity-stratified residue test (null comparison that does not depend on how many references cover a locus).
For every unique in-frame locus (03 definition: best alignment per reference, max window identity >= 0.3,
repeated copies collapsed by +-30 nt context), take the residue of the single reference with the highest window
identity among clean columns (no gap within 3 columns, residue not X). Report residue counts per codon, split by
that identity (<0.6, 0.6-0.8, >=0.8). Near-identical windows (>=0.8) flag recently shared / transferred DNA.
usage: 11_bestref.py LABEL [LABEL...] -> out/bestref.txt (appended)"""
import sys
from collections import Counter, defaultdict
sys.path.insert(0, "src/recheck")
from mp import read_fasta, parse_paf, best_per_ref, rc
W = "recheck_work"
prots = read_fasta(f"{W}/refs/refs_all.faa")
COD = ["CGG", "CGA", "CGC", "CGT", "AGG", "AGA"]
for lab in sys.argv[1:]:
    genome = read_fasta(f"{W}/genomes/{lab}.fna")
    alns = [a for a in parse_paf(f"{W}/out/{lab}.paf", genome, prots) if a["nM"] >= 30 and a["AS"] >= 60]
    best = best_per_ref(alns)
    loci = defaultdict(list)
    for (ctg, st, s, ref), (AS, ai, o, codon, raa, wid, neargap, kf, kt) in best.items():
        if codon in COD: loci[(ctg, st, s, codon)].append((ref, raa, wid, neargap))
    seen = set(); tab = {c: [Counter(), Counter(), Counter()] for c in COD}
    for (ctg, st, s, codon), v in sorted(loci.items()):
        if max(x[2] for x in v) < 0.3: continue
        x = genome[ctg][max(0, s - 30):s + 33]; x = x if st == "+" else rc(x)
        if x in seen: continue
        seen.add(x)
        cl = [y for y in v if not y[3] and y[1] not in "X*-"]
        if not cl: continue
        ref, raa, wid, _ = max(cl, key=lambda y: y[2])
        tab[codon][0 if wid < 0.6 else (1 if wid < 0.8 else 2)][raa] += 1
    with open(f"{W}/out/bestref.txt", "a") as f:
        for c in COD:
            parts = []
            for i, nm in enumerate(["id<0.6", "0.6-0.8", ">=0.8"]):
                t = tab[c][i]; n = sum(t.values())
                parts.append(f"{nm}: n={n} R={t['R']/n if n else 0:.2f} [{' '.join(f'{a}{b}' for a, b in t.most_common(4))}]")
            line = f"{lab} {c} " + " | ".join(parts)
            print(line, flush=True); f.write(line + "\n")
