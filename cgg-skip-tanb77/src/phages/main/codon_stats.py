"""Per-genome disputed-codon statistics (code-walls RULES.md).

usage: python -I codon_stats.py IDS.txt|- FASTA|- OUT.tsv [--codon CGG|TGA] [--table 11|4]
       [--mode meta|single|longest] [--whole NAME] [--plant f --seed s] [--minlen N] [--append]

For every sequence whose ID is in IDS.txt (or all with "-"): CDS calls (pyrodigal meta / single trained on the
sequence(s) / longest stop-to-stop ORF on the given strand), CDS >= 300 nt, not partial; in-frame codon counts.
O = disputed codon count; E = n(family) x third-position share (RULES); for CGG a homology-free skip signature:
an in-frame CGG within the last 16 codons of a CDS whose +1 frame (from the CGG's 2nd base) is open >= 50 codons.
--whole NAME aggregates all records of the FASTA into one row (genomes with several contigs).
--plant f replaces a share f of in-frame CGG by another CGN drawn from the genome's own CGN use (P3 recall).
"""
import sys
import random
import argparse
import pyrodigal
try:
    import pyaragorn
except ImportError:
    pyaragorn = None

STOPS = {11: {"TAA", "TAG", "TGA"}, 4: {"TAA", "TAG"}}
FOURFOLD_NONARG = ["GG", "GC", "CC", "AC", "GT", "CT", "TC"]
NAR = ["AA", "GA", "CA"]  # Lys, Glu, Gln two-fold purine-ending
COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def read_fasta(fh, keep):
    name, seq = None, []
    for line in fh:
        if line.startswith(">"):
            if name is not None and (keep is None or name in keep):
                yield name, "".join(seq).upper()
            name, seq = line[1:].split()[0], []
        else:
            seq.append(line.strip())
    if name is not None and (keep is None or name in keep):
        yield name, "".join(seq).upper()


def revcomp(s):
    return s[::-1].translate(COMP)


def open_len(seq, start, stops):
    n = 0
    for j in range(start, len(seq) - 2, 3):
        if seq[j:j + 3] in stops:
            return n
        n += 1
    return n


class G:
    def __init__(self, b, e, st):
        self.begin, self.end, self.strand = b, e, st
        self.partial_begin = self.partial_end = False


def longest_orf(seq, stops):
    best = None
    for f in range(3):
        st = f
        for j in range(f, len(seq) - 2, 3):
            if seq[j:j + 3] in stops:
                if best is None or j + 3 - st > best[1] - best[0]:
                    best = (st, j + 3)
                st = j + 3
    return [G(best[0] + 1, best[1], 1)] if best else []


class Acc:
    def __init__(self):
        self.counts, self.skip, self.ncds, self.totlen, self.gcn = {}, 0, 0, 0, 0
        self.trna = {"ccg": 0, "tcg": 0, "all": 0}


def add_seq(acc, seq, genes, a, stops, rng):
    acc.totlen += len(seq)
    acc.gcn += seq.count("G") + seq.count("C")
    rc = None
    for g in genes:
        if g.partial_begin or g.partial_end:
            continue
        s = seq[g.begin - 1:g.end]
        if g.strand == -1:
            s = revcomp(s)
        if len(s) < 300 or len(s) % 3:
            continue
        acc.ncds += 1
        cod = [s[i:i + 3] for i in range(0, len(s) - 3, 3)]  # terminal stop dropped
        if a.plant > 0 and a.codon == "CGG":
            alts = [c for c in cod if c[:2] == "CG" and c != "CGG"] or ["CGT", "CGA", "CGC"]
            cod = [(rng.choice(alts) if (c == "CGG" and rng.random() < a.plant) else c) for c in cod]
        for i, c in enumerate(cod):
            acc.counts[c] = acc.counts.get(c, 0) + 1
            if a.codon == "CGG" and c == "CGG" and len(cod) - i <= 16:
                if g.strand == 1:
                    ok = open_len(seq, g.begin - 1 + 3 * i + 1, stops) >= 50
                else:
                    if rc is None:
                        rc = revcomp(seq)
                    ok = open_len(rc, (len(seq) - g.end) + 3 * i + 1, stops) >= 50
                acc.skip += ok


def row(name, acc, a):
    c = acc.counts
    ncod = sum(c.values())
    O = c.get(a.codon, 0)
    if a.codon == "CGG":
        nfam = sum(v for k, v in c.items() if k[:2] == "CG")
        tot = sum(v for k, v in c.items() if k[:2] in FOURFOLD_NONARG)
        third = sum(v for k, v in c.items() if k[:2] in FOURFOLD_NONARG and k[2] == "G") / tot if tot else 0
    elif a.codon == "TGA":
        nfam = c.get("TGA", 0) + c.get("TGG", 0)
        tot = sum(v for k, v in c.items() if k[:2] in NAR and k[2] in "AG")
        third = sum(v for k, v in c.items() if k[:2] in NAR and k[2] == "A") / tot if tot else 0
    else:
        raise SystemExit("codon not supported: " + a.codon)
    E = nfam * third
    r = O / E if E > 0 else float("nan")
    sib = ""
    if a.codon == "CGG":  # same-family siblings: CGx observed vs n(CGN) x share of x at four-fold non-Arg thirds
        for x in "ATC":
            sx = sum(v for k, v in c.items() if k[:2] in FOURFOLD_NONARG and k[2] == x) / tot if tot else 0
            sib += f"\t{c.get('CG' + x, 0)}\t{nfam * sx:.3f}"
    else:
        sib = "\tNA" * 6
    gc = acc.gcn / max(1, acc.totlen)
    return (f"{name}\t{acc.totlen}\t{gc:.4f}\t{acc.ncds}\t{ncod}\t{O}\t{nfam}\t{third:.4f}\t{E:.3f}\t{r:.4f}\t"
            f"{acc.skip}\t{acc.trna['ccg']}\t{acc.trna['tcg']}\t{acc.trna['all']}{sib}\n")


HEADER = "id\tlen\tgc\tn_cds\tn_codons\tO\tn_fam\tthird_share\tE\tr\tskip_sites\ttrna_ccg\ttrna_tcg\ttrna_all\tO_CGA\tE_CGA\tO_CGT\tE_CGT\tO_CGC\tE_CGC\n"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ids"); ap.add_argument("fasta"); ap.add_argument("out")
    ap.add_argument("--codon", default="CGG")
    ap.add_argument("--table", type=int, default=11)
    ap.add_argument("--plant", type=float, default=0.0)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--minlen", type=int, default=5000)
    ap.add_argument("--mode", default="meta", choices=["meta", "single", "longest"])
    ap.add_argument("--whole", default=None)
    ap.add_argument("--append", action="store_true")
    a = ap.parse_args()
    keep = None if a.ids == "-" else set(l.split()[0] for l in open(a.ids) if l.strip())
    fh = sys.stdin if a.fasta == "-" else open(a.fasta)
    rng = random.Random(a.seed)
    stops = STOPS[a.table]
    meta = pyrodigal.GeneFinder(meta=True)
    rf = pyaragorn.RNAFinder(translation_table=a.table, tmrna=False) if pyaragorn else None

    def trna(acc, seq):
        if rf is None:
            return
        for g in rf.find_rna(seq.encode()):
            if getattr(g, "anticodon", None):
                acc.trna["all"] += 1
                ac = g.anticodon.lower()
                if ac in acc.trna:
                    acc.trna[ac] += 1
    out = open(a.out, "a" if a.append else "w")
    if not a.append:
        out.write(HEADER)

    def genes_for(seq, trained=None):
        if a.mode == "longest":
            return longest_orf(seq, stops)
        if trained is not None:
            return trained.find_genes(seq.encode())
        if a.mode == "single" and len(seq) >= 20000:
            g1 = pyrodigal.GeneFinder(meta=False)
            g1.train(seq.encode(), translation_table=a.table)
            return g1.find_genes(seq.encode())
        return meta.find_genes(seq.encode())

    if a.whole:
        recs = [(n, x) for n, x in read_fasta(fh, keep) if len(x) >= a.minlen]
        trained = None
        if a.mode == "single":
            trained = pyrodigal.GeneFinder(meta=False)
            trained.train(*[x.encode() for _, x in recs], translation_table=a.table)
        acc = Acc()
        for n, x in recs:
            add_seq(acc, x, genes_for(x, trained), a, stops, rng)
            trna(acc, x)
        out.write(row(a.whole, acc, a))
    else:
        for n, x in read_fasta(fh, keep):
            if len(x) < a.minlen:
                continue
            acc = Acc()
            add_seq(acc, x, genes_for(x), a, stops, rng)
            trna(acc, x)
            out.write(row(n, acc, a))
    out.close()


if __name__ == "__main__":
    main()
