"""pyaragorn tRNA calls with anticodon offset (codes_census/org_trnas.py plus ac_off/ac_len columns).
usage: python trnas.py out.tsv genomes.fa [table=4]. acceptor() copied unchanged from org_trnas.py."""
import sys
import pyaragorn

PAIRS = {("G", "C"), ("C", "G"), ("A", "T"), ("T", "A"), ("G", "T"), ("T", "G")}
TAIL_RANK = {"CCA": 0, "": 1, "CC": 2, "C": 3}


def acceptor(seq):
    s = seq.upper()
    best = None
    for o in range(0, 3):
        stem5 = s[o:o + 7]
        for t in range(0, 6):
            e = len(s) - t
            if e - o < 50:
                break
            stem3, tail = s[e - 8:e - 1], s[e:]
            wc = sum((stem5[i], stem3[6 - i]) in PAIRS for i in range(7))
            tloop = s[e - 20:e - 17] == "TTC"
            key = (-(wc + 3 * tloop), o, TAIL_RANK.get(tail, 4), t)
            if best is None or key < best[0]:
                best = (key, stem5, stem3, s[e - 1], tail == "CCA", wc)
    if best is None:
        return s[:7], "", "?:?", "?", False, 0
    _, stem5, stem3, n73, cca, wc = best
    return stem5, stem3, f"{stem5[2]}:{stem3[4]}", n73, cca, wc


def fasta(path):
    name, buf = None, []
    for line in open(path):
        if line.startswith('>'):
            if name:
                yield name, ''.join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name:
        yield name, ''.join(buf)


def main(out, fa, table=4):
    finder = pyaragorn.RNAFinder(int(table), tmrna=False, linear=False)
    with open(out, 'w') as fh:
        fh.write('acc\tbegin\tend\tstrand\taa\tanticodon\tac_off\tac_len\tstem5\tstem3\tpair3_70\tN73\tstem_pairs\tenergy\tseq\n')
        for acc, s in fasta(fa):
            for g in finder.find_rna(s):
                if g.type != 'tRNA':
                    continue
                seq = g.sequence() if callable(g.sequence) else g.sequence
                st5, st3, p, n73, cca, wc = acceptor(seq)
                fh.write(f'{acc}\t{g.begin}\t{g.end}\t{g.strand}\t{g.amino_acid}\t{g.anticodon.upper()}\t{g.anticodon_offset}\t'
                         f'{g.anticodon_length}\t{st5}\t{st3}\t{p}\t{n73}\t{wc}\t{g.energy:.1f}\t{seq}\n')


if __name__ == '__main__':
    main(*sys.argv[1:])
