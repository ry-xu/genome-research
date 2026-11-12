"""Stop-codon context: bases right after TTA-at-e*+1 vs after the clade's canonical stops vs after internal Leu codons.

usage: context.py <panel.tsv> <genes_dir> <genomes_dir> <ortho_dir>

Clade genomes only. Sites (gene orientation), from runon_events.tsv with R = sister:
  TTA_end   run-on events with TTA at e*+1 (one per gene)
  STOP_end  "same" events: the canonical stop at e*+1, i.e. ending where the sister ortholog ends
  STOP_all  the canonical stop of every complete gene
  LEU_int   CTG / CTC / TTG codons at residues 21 .. len-21 of complete genes (sense context)
  ANY_int   all codons at residues 21 .. len-21 (next-codon first-base composition)
Per class: base frequencies at +4, +5, +6 (first three bases after the codon) and -1 (last base of the codon before).
G-test (likelihood ratio, 3 df) of TTA_end against STOP_end and against LEU_int at each position, with per-genome
pooling (counts summed over the 8 clade genomes). Also: share of TTA_end followed by A or G (a TAA / TAG would then
sit one base downstream, in frame +1).
"""
import gzip
import math
import sys
from collections import Counter, defaultdict

COMP = str.maketrans("ACGTN", "TGCAN")
STOPS = {"TAA", "TAG", "TGA"}


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


def gtest(a, b):
    """G statistic and chi2(3) upper-tail p for two count vectors over ACGT."""
    na, nb = sum(a), sum(b)
    g = 0.0
    for x, y in zip(a, b):
        tot = x + y
        for o, n in ((x, na), (y, nb)):
            e = tot * n / (na + nb)
            if o > 0:
                g += 2 * o * math.log(o / e)
    # chi2 3 df survival: exact formula
    x = g / 2
    # chi2(3) survival at G: erfc(sqrt(G/2)) + sqrt(2G/pi) exp(-G/2)
    p = math.erfc(math.sqrt(x)) + 2 * math.sqrt(x / math.pi) * math.exp(-x)
    return g, p


def main(panel_fn, gdir, gendir, odir):
    panel = [l.rstrip("\n").split("\t") for l in open(panel_fn)]
    clade = [r[0] for r in panel if r[2] == "clade"]
    acc = {r[0]: r[1] for r in panel}
    ev = defaultdict(dict)
    with open(f"{odir}/runon_events.tsv") as f:
        h = next(f).rstrip("\n").split("\t")
        for line in f:
            r = dict(zip(h, line.rstrip("\n").split("\t")))
            if r["R"] == "sister" and r["B_role"] == "clade":
                if r["class"] == "runon" and r["w_p1"] == "TTA":
                    ev[r["B"]][r["gene"]] = ("TTA_end", int(r["e_star"]))
                elif r["class"] == "same":
                    ev[r["B"]][r["gene"]] = ("STOP_end", int(r["e_star"]))
    cnt = defaultdict(lambda: defaultdict(Counter))  # cls -> pos -> Counter
    per_genome_tta_next = Counter()
    for lab in clade:
        contigs = dict(read_fasta(f"{gendir}/{acc[lab]}.fna.gz"))
        with open(f"{gdir}/{lab}.genes.tsv") as f:
            next(f)
            for line in f:
                gid, ctg, b, e, st, pb, pe, nt = line.rstrip("\n").split("\t")
                if pb == "1" or pe == "1" or len(nt) % 3:
                    continue
                b, e, st = int(b), int(e), int(st)
                ext = downstream(contigs[ctg], b, e, st, 0, len(nt) + 6)
                if len(ext) < len(nt) + 3:
                    continue
                L = len(nt)
                # all canonical stops
                for k, pos in enumerate((4, 5, 6)):
                    if L + k < len(ext):
                        cnt["STOP_all"][pos][ext[L + k]] += 1
                cnt["STOP_all"][-1][nt[L - 4]] += 1
                ncod = L // 3
                for i in range(21, ncod - 21):
                    c = nt[3 * i:3 * i + 3]
                    nxt = nt[3 * i + 3:3 * i + 6]
                    cnt["ANY_int"][4][nxt[0]] += 1
                    if c in ("CTG", "CTC", "TTG"):
                        for k, pos in enumerate((4, 5, 6)):
                            cnt["LEU_int"][pos][nxt[k]] += 1
                        cnt["LEU_int"][-1][nt[3 * i - 1]] += 1
                if gid in ev[lab]:
                    cl, es = ev[lab][gid]
                    s = 3 * es  # codon at e*+1 starts here
                    cod = ext[s:s + 3]
                    assert (cod == "TTA") if cl == "TTA_end" else (cod in STOPS), (lab, gid, cod)
                    for k, pos in enumerate((4, 5, 6)):
                        cnt[cl][pos][ext[s + 3 + k]] += 1
                    cnt[cl][-1][ext[s - 1]] += 1
    order = "ACGT"
    print("class\tpos\tn\tA\tC\tG\tT")
    for cl in ("TTA_end", "STOP_end", "STOP_all", "LEU_int", "ANY_int"):
        for pos in (-1, 4, 5, 6):
            c = cnt[cl][pos]
            n = sum(c.values())
            if n:
                print(f"{cl}\t{pos}\t{n}\t" + "\t".join(f"{c[x] / n:.3f}" for x in order))
    print()
    print("test\tpos\tG\tp")
    for other in ("STOP_end", "STOP_all", "LEU_int"):
        for pos in (-1, 4, 5, 6):
            a = [cnt["TTA_end"][pos][x] for x in order]
            b = [cnt[other][pos][x] for x in order]
            if sum(b):
                g, p = gtest(a, b)
                print(f"TTA_end vs {other}\t{pos}\t{g:.2f}\t{p:.3g}")
    for other in ("LEU_int", "ANY_int"):
        a = [cnt["STOP_end"][4][x] for x in order]
        b = [cnt[other][4][x] for x in order]
        g, p = gtest(a, b)
        print(f"STOP_end vs {other}\t4\t{g:.2f}\t{p:.3g}")


if __name__ == "__main__":
    main(*sys.argv[1:5])
