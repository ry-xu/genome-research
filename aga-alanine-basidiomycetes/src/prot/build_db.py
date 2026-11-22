"""Build the protein database for the AGA = Ala proteomics test.
usage: python build_db.py mito.fa ACC[,ACC...] nuclear.fasta out.fasta sites.tsv
Mitochondrial ORFs: all six frames of each circular genome, stop to stop under table 4 (TAA/TAG stop), >= 50 codons.
Each ORF is written three times: |R (table 4, AGA = Arg), |A (AGA = Ala) and |C (control: AGA = Arg, CGN = Ala).
An engine that favours Ala versions for search reasons alone would favour |C too; a real AGA = Ala favours only |A.
sites.tsv lists, per ORF, the protein index of each AGA and CGN codon."""
import sys

B = 'TCAG'
AA4 = 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
T4 = {a + b + c: AA4[16 * i + 4 * j + k] for i, a in enumerate(B) for j, b in enumerate(B) for k, c in enumerate(B)}
RC = str.maketrans('ACGTacgt', 'TGCAtgca')


def read_fa(path, want):
    seqs, cur = {}, None
    for line in open(path):
        if line.startswith('>'):
            acc = line[1:].split()[0]
            cur = acc if acc in want else None
            if cur:
                seqs[cur] = []
        elif cur:
            seqs[cur].append(line.strip().upper())
    return {k: ''.join(v) for k, v in seqs.items()}


def orfs(seq, minc=50):
    """Stop-to-stop ORFs on a circular sequence; yields (strand, frame, start_nt, codons)."""
    n = len(seq)
    for strand, s in (('+', seq), ('-', seq.translate(RC)[::-1])):
        ss = s + s  # circular: walk twice the length, keep ORFs that start in the first copy
        for f in range(3):
            cod = [ss[i:i + 3] for i in range(f, 2 * n - 2, 3)]
            stops = [i for i, c in enumerate(cod) if T4.get(c) == '*']
            for a, b in zip(stops, stops[1:]):
                st = f + 3 * (a + 1)
                if st >= n or b - a - 1 < minc:
                    continue
                yield strand, f, st, cod[a + 1:b]


def main(mito_fa, accs, nuc_fa, out_fa, sites_tsv):
    seqs = read_fa(mito_fa, set(accs.split(',')))
    with open(out_fa, 'w') as out, open(sites_tsv, 'w') as st:
        st.write('orf\tcodon\tpos\n')
        k = 0
        for acc, seq in seqs.items():
            for strand, f, start, cod in orfs(seq):
                if any(T4.get(c, 'X') == 'X' for c in cod):
                    continue
                name = f'{acc}_{strand}{f}_{start}'
                r = ''.join(T4[c] for c in cod)
                a = ''.join('A' if c == 'AGA' else T4[c] for c in cod)
                c4 = ''.join('A' if c.startswith('CG') else T4[c] for c in cod)
                for tag, p in (('R', r), ('A', a), ('C', c4)):
                    out.write(f'>mito|{name}|{tag}\n{p}\n')
                for i, c in enumerate(cod):
                    if c == 'AGA' or c.startswith('CG'):
                        st.write(f'{name}\t{c}\t{i}\n')
                k += 1
        print(f'{k} mito ORFs from {len(seqs)} genomes', file=sys.stderr)
        n = 0
        for line in open(nuc_fa):
            if line.startswith('>'):
                n += 1
                out.write(f'>nuc|{line[1:].split()[0]}\n')
            else:
                out.write(line)
        print(f'{n} nuclear proteins', file=sys.stderr)


if __name__ == '__main__':
    main(*sys.argv[1:])
