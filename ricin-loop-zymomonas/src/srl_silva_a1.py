"""Call the sarcin-ricin tetraloop (E. coli 2659-2662) in a SILVA LSU export, rule A1 of ../RULES_srl_silva.md:
an exact AGTAC NNNN GGACC match starting at >= 70% of the sequence length; 0 matches -> unread (empty loop),
more than one -> "multi". usage: srl_silva_a1.py SILVA_..._tax_silva.fasta.gz OUT.tsv"""
import gzip, re, sys

RX = re.compile("AGTAC([ACGT]{4})GGACC")


def records(path):
    head, buf = None, []
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                if head:
                    yield head, "".join(buf)
                head, buf = line[1:].rstrip(), []
            else:
                buf.append(line.strip())
    if head:
        yield head, "".join(buf)


def main(path, out):
    with open(out, "w") as fo:
        fo.write("acc\tlen\tpos\tloop\ttaxonomy\n")
        for head, seq in records(path):
            acc, _, tax = head.partition(" ")
            seq = seq.upper().replace("U", "T")
            hits = [m for m in RX.finditer(seq) if m.start() >= 0.7 * len(seq)]
            if len(hits) == 1:
                pos, loop = hits[0].start(), hits[0].group(1)
            else:
                pos, loop = "", ("multi" if hits else "")
            fo.write(f"{acc}\t{len(seq)}\t{pos}\t{loop}\t{tax}\n")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
