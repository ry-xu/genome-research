"""Are interior TTAs in clade genes more frequent than the background of premature canonical stops?

usage: stop_background.py <panel.tsv> <genes_dir> <genomes_dir> <ortho_dir> <out.tsv>

If TTA is a stop (or an unassigned codon), an interior TTA in a conserved gene is a nonsense mutation, a pseudogene
or an assembly error, and should occur at about the rate of premature canonical stops (TAA/TAG/TGA) times the
mutational opportunity ratio. If TTA is still read as Leu at useful efficiency, interior TTAs should be in excess.

Same test for both codon kinds ("X" = the stop-like codon), each clade gene vs its best sister RBH partner:
  TTA   complete gene, TTA at residue r, >= 30 DIAMOND-aligned residues before it (qs <= r - 30)
  STOP  complete gene (pe == 0) ending at a canonical stop, DIAMOND alignment reaches within 2 residues of the
        gene's C-terminus and spans >= 30 residues
  CTA   control sense codon, same rule as TTA (sensitivity of the test for a real sense codon); up to 400 sites per
        genome, every k-th site
For each X: partner residues after the one aligned with X (TTA/CTA: via a global BLOSUM62 alignment; STOP: from the
DIAMOND end coordinates); require >= 30 of them. The clade contig after X, in frame, is translated (stops as '*') for
W + 10 codons, W = min(100, partner remainder). Local BLOSUM62 alignment (open -11, extend -1) of that continuation
to the first W partner residues. "Homologous continuation" = score >= the 99th percentile of the same score with the
continuation shuffled (one shuffle per site, pooled over all sites and genomes, seed 1) and the local alignment
starting within 10 residues of both starts (colinear).
Opportunity: sum over interior sense codons of the number of single-nucleotide changes that create TTA vs that
create a canonical stop (uniform mutation weights).
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


def downstream(contig, b, e, strand, off, n):
    if strand == 1:
        s = b - 1 + off
        return contig[s:s + n]
    end0 = e - 1 - off
    lo = max(0, end0 - n + 1)
    return rc(contig[lo:end0 + 1])


def tr(s):
    return "".join(CODE.get(s[i:i + 3], "X") for i in range(0, len(s) - 2, 3))


def one_step(c, targets):
    n = 0
    for i in range(3):
        for b in "ACGT":
            if b != c[i] and c[:i] + b + c[i + 1:] in targets:
                n += 1
    return n


def main(panel_fn, gdir, gendir, odir, out_fn):
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
    gal = Align.PairwiseAligner()
    gal.substitution_matrix = substitution_matrices.load("BLOSUM62")
    gal.open_gap_score, gal.extend_gap_score = -10, -0.5
    gal.mode = "global"
    gal.end_gap_score = 0
    lal = Align.PairwiseAligner()
    lal.substitution_matrix = substitution_matrices.load("BLOSUM62")
    lal.open_gap_score, lal.extend_gap_score = -11, -1
    lal.mode = "local"
    rng = random.Random(1)
    sites = []  # (lab, kind, codon, gene, res, ql, pid, score, null, start_c, start_r, W)
    opp = {}
    for lab in labs:
        contigs = dict(read_fasta(f"{gendir}/{acc[lab]}.fna.gz"))
        cu = Counter()
        cta_all = []
        tta_all = []
        stop_all = []
        for g, (ctg, b, e, st, okb, oke, nt) in genes.items():
            if not g.startswith(lab + "|"):
                continue
            ncod = len(nt) // 3
            if okb and oke and len(nt) % 3 == 0:
                for i in range(21, ncod - 21):
                    cu[nt[3 * i:3 * i + 3]] += 1
            if g not in best or len(nt) % 3:
                continue
            bits, s, pid, qs, qe, ss, se, ql, sl = best[g]
            if okb and oke:
                for r in range(31, ql + 1):
                    c = nt[3 * (r - 1):3 * r]
                    if c in ("TTA", "CTA") and qs <= r - 30:
                        (tta_all if c == "TTA" else cta_all).append((g, r))
            if oke and nt[-3:] in STOPS and qe >= ql - 2 and qe - qs + 1 >= 30:
                stop_all.append((g, ql + 1))
        step = max(1, len(cta_all) // 400)
        cands = [("TTA", x) for x in tta_all] + [("CTA", x) for x in cta_all[::step]] + [("STOP", x) for x in stop_all]
        if role[lab] == "sister":  # TTA is ordinary Leu there; subsample it like CTA
            st2 = max(1, len(tta_all) // 400)
            cands = [("TTA", x) for x in tta_all[::st2]] + [c for c in cands if c[0] != "TTA"]
        prot_cache = {}

        def prot(gg):
            if gg not in prot_cache:
                prot_cache[gg] = tr(genes[gg][6][:-3])
            return prot_cache[gg]

        for kind, (g, r) in cands:
            ctg, b, e, st, okb, oke, nt = genes[g]
            bits, s, pid, qs, qe, ss, se, ql, sl = best[g]
            pp = prot(s)
            if kind == "STOP":
                P = se + (ql - qe) + 1  # 1-based partner residue aligned with the stop codon
                rem = pp[P:]
            else:
                p = prot(g)
                a = gal.align(p, pp)[0]
                c1, c2 = a.coordinates
                m = {}
                for kk in range(len(c1) - 1):
                    if c1[kk + 1] - c1[kk] and c2[kk + 1] - c2[kk]:
                        for j in range(c1[kk + 1] - c1[kk]):
                            m[c1[kk] + j] = c2[kk] + j
                t = None
                for i in range(r - 1, max(-1, r - 12), -1):
                    if i in m:
                        t = m[i] + (r - 1 - i)
                        break
                if t is None:
                    continue
                rem = pp[t + 1:]
            if len(rem) < 30:
                continue
            W = min(100, len(rem))
            cont_nt = downstream(contigs[ctg], b, e, st, 3 * r, 3 * (W + 10))
            if len(cont_nt) < 90:
                continue
            cont = tr(cont_nt)
            ref = rem[:W]
            sc = lal.score(cont, ref)
            al = lal.align(cont, ref)[0]
            s0c, s0r = al.coordinates[0][0], al.coordinates[1][0]
            sh = list(cont)
            rng.shuffle(sh)
            null = lal.score("".join(sh), ref)
            codon = nt[3 * (r - 1):3 * r] if kind != "STOP" else nt[-3:]
            sites.append((lab, kind, codon, g.split("|", 1)[1], r, ql, pid, sc, null, s0c, s0r, W))
        n_tta = sum(v * one_step(c, {"TTA"}) for c, v in cu.items() if c not in STOPS and c != "TTA")
        n_stop = sum(v * one_step(c, STOPS) for c, v in cu.items() if c not in STOPS)
        opp[lab] = (n_tta, n_stop, sum(v for c, v in cu.items()))
        print(f"{lab}: {len(cands)} candidates", file=sys.stderr)
    nulls = sorted(x[8] for x in sites)
    thr = nulls[int(0.99 * len(nulls))]
    with open(out_fn, "w") as f:
        f.write("genome\trole\tkind\tcodon\tgene\tres\tql\tpident\tscore\tnull\tstart_c\tstart_r\tW\thomologous\n")
        for x in sites:
            hom = int(x[7] >= thr and x[9] <= 10 and x[10] <= 10)
            f.write("\t".join(map(str, x[:1] + (role[x[0]],) + x[1:])) + f"\t{hom}\n")
    print(f"threshold (99th pct of shuffled scores, n={len(nulls)}): {thr}")
    print("genome\trole\tTTA_sites\tTTA_hom\tCTA_sites\tCTA_hom\tCTA_sens\tSTOP_sites\tSTOP_hom\tSTOP_hom_by_codon\t"
          "opp_TTA/opp_STOP\tE_TTA_hom_if_stoplike\tobs/E")
    agg = defaultdict(lambda: defaultdict(lambda: [0, 0]))
    bycod = defaultdict(Counter)
    for x in sites:
        hom = x[7] >= thr and x[9] <= 10 and x[10] <= 10
        agg[x[0]][x[1]][0] += 1
        agg[x[0]][x[1]][1] += hom
        if x[1] == "STOP" and hom:
            bycod[x[0]][x[2]] += 1
    for lab in labs:
        d = agg[lab]
        sens = d["CTA"][1] / d["CTA"][0] if d["CTA"][0] else float("nan")
        ratio = opp[lab][0] / opp[lab][1]
        E = d["STOP"][1] * ratio
        print(f"{lab}\t{role[lab]}\t{d['TTA'][0]}\t{d['TTA'][1]}\t{d['CTA'][0]}\t{d['CTA'][1]}\t{sens:.2f}\t"
              f"{d['STOP'][0]}\t{d['STOP'][1]}\t{dict(bycod[lab])}\t{ratio:.3f}\t{E:.2f}\t"
              f"{(d['TTA'][1] / E) if E else float('nan'):.2f}")


if __name__ == "__main__":
    main(*sys.argv[1:6])
