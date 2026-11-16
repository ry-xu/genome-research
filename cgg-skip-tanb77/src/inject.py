"""Plant frame breaks in a real genome to measure recall, or homopolymer indels to model sequencing error.

usage: python -I src/inject.py <genome.fna.gz> <out.fna.gz> <truth.tsv> <mode> <frac> <seed>
  mode plus1  : in a fraction <frac> of genes (pyrodigal, >= 300 nt), insert a T before an interior codon starting
                AA/AG/GA, so the old frame reads a stop at the site and the gene continues +1 (Euplotes type)
  mode minus1 : in a fraction <frac> of genes, delete one nt at an interior codon (gene continues -1)
  mode plus1at:XYZ   : as plus1, but only where the codon before the site is XYZ (a recurring Euplotes-type motif)
  mode quad:XYZ:N    : insert base N right after an interior codon XYZ (a quadruplet codon XYZN; no stop needed)
  mode minus1at:XYZ  : delete the first nt of the codon right after an interior codon XYZ (a recurring -1 site)
  mode hp     : <frac> indels per kb of genome, each a 1-nt deletion or insertion inside a homopolymer run >= 4
                (the nanopore/454 error class), on either strand, genes or not
Truth rows: contig, gene strand, forward 1-based position of the event in the OUTPUT genome, type (+1 / -1 as the
scanner reports them: +1 = genome has an extra nt), gene id (or '.').
Interior = codon 25% to 75% along the gene, so both flanks carry >= 25% of the protein.
"""
import gzip
import random
import re
import sys

import pyrodigal

COMP = str.maketrans("ACGTN", "TGCAN")


def rc(s):
    return s.translate(COMP)[::-1]


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf).upper()
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf).upper()
    return seqs


def main():
    src, out, truth, mode, frac, seed = sys.argv[1:7]
    frac = float(frac)
    rng = random.Random(int(seed))
    g = read_fasta(src)
    edits = {c: [] for c in g}   # contig -> [(fwd 0-based pos, 'ins'/'del', base, strand, type, gene)]
    kind, *arg = mode.split(":")
    if kind in ("plus1", "minus1", "plus1at", "quad", "minus1at"):
        finder = pyrodigal.GeneFinder()
        finder.train(*[s.encode() for s in g.values() if len(s) > 20000][:50])
        genes = []
        for c, s in g.items():
            for i, gn in enumerate(finder.find_genes(s.encode())):
                if gn.end - gn.begin + 1 >= 300 and not gn.partial_begin and not gn.partial_end:
                    genes.append((c, gn.begin, gn.end, "+" if gn.strand == 1 else "-", f"{c}_{i + 1}"))
        pick = rng.sample(genes, max(1, round(frac * len(genes))))
        for c, b, e, strand, gid in pick:
            gs = g[c][b - 1:e] if strand == "+" else rc(g[c][b - 1:e])
            ncod = len(gs) // 3
            cods = list(range(ncod // 4, 3 * ncod // 4))
            rng.shuffle(cods)
            for k in cods:
                cod = gs[3 * k:3 * k + 3]
                prev = gs[3 * k - 3:3 * k]
                if kind in ("plus1", "plus1at") and cod[:2] not in ("AA", "AG", "GA"):
                    continue
                if kind in ("plus1at", "quad", "minus1at") and prev != arg[0]:
                    continue
                ins = {"plus1": "T", "plus1at": "T", "quad": arg[-1] if arg else ""}.get(kind)
                gp = 3 * k                      # gene-strand offset of the codon start
                if strand == "+":
                    fwd = b - 1 + gp            # insert before this forward position
                    if ins:
                        edits[c].append((fwd, "ins", ins, strand, "+1", gid))
                    else:
                        edits[c].append((fwd, "del", "", strand, "-1", gid))
                else:
                    # gene-strand offset gp maps to forward index e-1-gp; inserting before it on the gene strand
                    # means inserting after it on the forward strand
                    fwd = e - 1 - gp
                    if ins:
                        edits[c].append((fwd + 1, "ins", rc(ins), strand, "+1", gid))
                    else:
                        edits[c].append((fwd, "del", "", strand, "-1", gid))
                break
    elif kind == "hp":
        for c, s in g.items():
            runs = [m.start() for m in re.finditer(r"A{4,}|C{4,}|G{4,}|T{4,}", s)]
            n = min(len(runs), round(frac * len(s) / 1000))
            for st in rng.sample(runs, n):
                if rng.random() < 0.5:
                    edits[c].append((st, "del", "", ".", "-1", "."))
                else:
                    edits[c].append((st, "ins", s[st], ".", "+1", "."))
    else:
        sys.exit(f"bad mode {mode}")
    with gzip.open(out, "wt", compresslevel=3) as fo, open(truth, "w") as ft:
        ft.write("contig\tstrand\tpos\ttype\tgene\n")
        for c, s in g.items():
            ev = sorted(edits[c])
            parts, last, shift = [], 0, 0
            for fwd, kind, base, strand, typ, gid in ev:
                parts.append(s[last:fwd])
                if kind == "ins":
                    parts.append(base)
                    last = fwd
                    newpos = fwd + shift + 1    # 1-based position of the inserted base
                    shift += 1
                else:
                    last = fwd + 1
                    newpos = fwd + shift + 1    # 1-based position of the base after the deletion
                    shift -= 1
                ft.write(f"{c}\t{strand}\t{newpos}\t{typ}\t{gid}\n")
            parts.append(s[last:])
            fo.write(f">{c}\n{''.join(parts)}\n")


if __name__ == "__main__":
    main()
