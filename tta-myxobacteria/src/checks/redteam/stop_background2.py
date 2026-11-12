"""Add premature canonical stops found from the 3' side, to remove an ascertainment asymmetry in stop_background.py.

usage: stop_background2.py <panel.tsv> <genes_dir> <genomes_dir> <ortho_dir> <stop_background.tsv> <threshold> <out.tsv>

stop_background.py finds a premature canonical stop only when the 5' fragment is a called gene whose best sister RBH
partner is the full-length ortholog. Pyrodigal reads through TTA, so interior TTAs never need that. Here, for every
gene whose DIAMOND alignment to its best sister partner starts within 10 residues of the gene's N-terminus while the
partner has >= 30 residues before the aligned start: walk upstream from the start codon in frame to the nearest
canonical stop X. Translate the W + 10 codons before X (W = min(100, partner residues before the one aligned with X)),
local-align (BLOSUM62, open -11, extend -1) to the W partner residues before X's partner position. Homologous =
score >= threshold (the 99th percentile of shuffled scores from stop_background.py) and the alignment ends within 10
residues of both ends (colinear).
Union with the 5'-side STOP events, deduplicated by the stop codon's genomic position. Reports per genome: TTA_hom
(from stop_background.tsv), STOP_hom union, STOP_hom without TGA (selenocysteine possible), and expected TTA_hom if TTA
behaved like a canonical stop (STOP_hom x opportunity ratio, with and without TGA).
"""
import gzip
import random
import sys
import warnings
from collections import Counter, defaultdict

from Bio import Align
from Bio.Align import substitution_matrices

warnings.filterwarnings("ignore")
COMP = str.maketrans("ACGTN", "TGCAN")
STOPS = {"TAA", "TAG", "TGA"}
B = "TCAG"
AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODE = {}
k = 0
for x in B:
    for y in B:
        for z in B:
            CODE[x + y + z] = AAS[k]
            k += 1


def rc(s):
    return s.translate(COMP)[::-1]


def read_fasta(fn):
    op = gzip.open if fn.endswith(".gz") else open
    name, buf = None, []
    with op(fn, "rt") as f:
        for line in f:
            if line.startswith(">"):
                if name:
                    yield name, "".join(buf).upper()
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name:
        yield name, "".join(buf).upper()


def tr(s):
    return "".join(CODE.get(s[i:i + 3], "X") for i in range(0, len(s) - 2, 3))


def one_step(c, targets):
    n = 0
    for i in range(3):
        for b in "ACGT":
            if b != c[i] and c[:i] + b + c[i + 1:] in targets:
                n += 1
    return n


def main(panel_fn, gdir, gendir, odir, sb_fn, thr, out_fn):
    thr = float(thr)
    panel = [l.rstrip("\n").split("\t") for l in open(panel_fn)]
    labs = [r[0] for r in panel]
    acc = {r[0]: r[1] for r in panel}
    role = {r[0]: ("clade" if r[2] == "clade" else "sister") for r in panel}
    genes = {}
    for lab in labs:
        with open(f"{gdir}/{lab}.genes.tsv") as f:
            next(f)
            for line in f:
                gid, ctg, b, e, st, pb, pe, nt = line.rstrip("\n").split("\t")
                genes[f"{lab}|{gid}"] = (ctg, int(b), int(e), int(st), pb == "0", pe == "0", nt)
    best = {}
    with open(f"{odir}/rbh.tsv") as f:
        next(f)
        for line in f:
            q, s, pid, qs, qe, ss, se, ql, sl, bits = line.rstrip("\n").split("\t")
            if role[s.split("|")[0]] != "sister":
                continue
            bits = float(bits)
            if q not in best or bits > best[q][0]:
                best[q] = (bits, s, float(pid), int(qs), int(qe), int(ss), int(se), int(ql), int(sl))
    # 5'-side events from stop_background.tsv
    ev = defaultdict(dict)  # lab -> key -> (codon, side, score)
    tta_hom = Counter()
    with open(sb_fn) as f:
        h = next(f).rstrip("\n").split("\t")
        for line in f:
            r = dict(zip(h, line.rstrip("\n").split("\t")))
            if r["homologous"] != "1":
                continue
            if r["kind"] == "TTA":
                tta_hom[r["genome"]] += 1
            elif r["kind"] == "STOP":
                ctg, b, e, st, *_ = genes[f"{r['genome']}|{r['gene']}"]
                pos = e - 2 if st == 1 else b
                ev[r["genome"]][(ctg, pos)] = (r["codon"], "5p", float(r["score"]))
    lal = Align.PairwiseAligner()
    lal.substitution_matrix = substitution_matrices.load("BLOSUM62")
    lal.open_gap_score, lal.extend_gap_score = -11, -1
    lal.mode = "local"
    rng = random.Random(1)
    rows = []
    opp = {}
    for lab in labs:
        contigs = dict(read_fasta(f"{gendir}/{acc[lab]}.fna.gz"))
        cu = Counter()
        for g, (ctg, b, e, st, okb, oke, nt) in genes.items():
            if not g.startswith(lab + "|"):
                continue
            if okb and oke and len(nt) % 3 == 0:
                for i in range(21, len(nt) // 3 - 21):
                    cu[nt[3 * i:3 * i + 3]] += 1
            if not okb or g not in best:
                continue
            bits, s, pid, qs, qe, ss, se, ql, sl = best[g]
            if qs > 10 or ss - qs < 30:
                continue
            cs = contigs[ctg]
            n = 3 * 700
            if st == 1:
                up = cs[max(0, b - 1 - n):b - 1]
            else:
                up = rc(cs[e:e + n])
            up = up[len(up) % 3:]  # keep whole codons ending at the start
            cods = [up[i:i + 3] for i in range(0, len(up), 3)]
            j = None
            for kk in range(1, len(cods) + 1):
                if cods[-kk] in STOPS:
                    j = kk
                    break
            if j is None:
                continue
            P0 = ss - qs - j  # 0-based partner index aligned with X
            if P0 < 30:
                continue
            W = min(100, P0)
            seg = cods[max(0, len(cods) - j - (W + 10)):len(cods) - j]
            if len(seg) < 30:
                continue
            cont = tr("".join(seg))
            pp = tr(genes[s][6][:-3])
            ref = pp[P0 - W:P0]
            sc = lal.score(cont, ref)
            al = lal.align(cont, ref)[0]
            e_c, e_r = al.coordinates[0][-1], al.coordinates[1][-1]
            sh = list(cont)
            rng.shuffle(sh)
            null = lal.score("".join(sh), ref)
            hom = sc >= thr and len(cont) - e_c <= 10 and len(ref) - e_r <= 10
            # genomic position of X (leftmost base)
            if st == 1:
                pos = b - 3 * j
            else:
                pos = e + 3 * j - 2
            codon = cods[-j]
            rows.append((lab, role[lab], g.split("|", 1)[1], codon, j, P0, pid, sc, null, int(hom)))
            if hom and (ctg, pos) not in ev[lab]:
                ev[lab][(ctg, pos)] = (codon, "3p", sc)
            elif hom:
                ev[lab][(ctg, pos)] = (codon, "both", sc)
        n_tta = sum(v * one_step(c, {"TTA"}) for c, v in cu.items() if c not in STOPS and c != "TTA")
        n_stop = sum(v * one_step(c, STOPS) for c, v in cu.items() if c not in STOPS)
        n_stop2 = sum(v * one_step(c, {"TAA", "TAG"}) for c, v in cu.items() if c not in STOPS)
        opp[lab] = (n_tta, n_stop, n_stop2)
    with open(out_fn, "w") as f:
        f.write("genome\trole\tgene\tcodon\tcodons_upstream_of_start\tpartner_idx\tpident\tscore\tnull\thomologous\n")
        for r in rows:
            f.write("\t".join(map(str, r)) + "\n")
    nulls = sorted(r[8] for r in rows)
    print(f"3'-side candidates {len(rows)}; their shuffled-score 99th pct {nulls[int(0.99 * len(nulls))]} "
          f"(threshold used {thr}); shuffled >= thr: {sum(1 for x in nulls if x >= thr)}")
    print("genome\trole\tTTA_hom\tSTOP_hom_union\tby_side\tby_codon\tE_all\tE_noTGA\tobs/E_all\tobs/E_noTGA")
    tot = defaultdict(float)
    for lab in labs:
        d = ev[lab]
        sides = Counter(v[1] for v in d.values())
        cods = Counter(v[0] for v in d.values())
        nall = len(d)
        nno = nall - cods["TGA"]
        Ea = nall * opp[lab][0] / opp[lab][1]
        En = nno * opp[lab][0] / opp[lab][2]
        o = tta_hom[lab]
        if role[lab] == "clade":
            tot["o"] += o
            tot["Ea"] += Ea
            tot["En"] += En
            tot["n"] += nall
            tot["nno"] += nno
        print(f"{lab}\t{role[lab]}\t{o}\t{nall}\t{dict(sides)}\t{dict(cods)}\t{Ea:.2f}\t{En:.2f}\t"
              f"{(o / Ea) if Ea else float('nan'):.2f}\t{(o / En) if En else float('nan'):.2f}")
    print(f"clade total: TTA_hom {tot['o']:.0f}; STOP_hom {tot['n']:.0f} (no TGA {tot['nno']:.0f}); "
          f"E_all {tot['Ea']:.2f}; E_noTGA {tot['En']:.2f}")


if __name__ == "__main__":
    main(*sys.argv[1:8])
