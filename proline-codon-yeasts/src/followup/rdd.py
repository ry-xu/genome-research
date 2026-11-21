"""RNA-DNA differences (RDD) in coding sequence, RULES.md arm R.
usage:
  python rdd.py genes_gb  ref.gb  out_genes.tsv            # CDS from GenBank (joins kept, codon position along CDS)
  python rdd.py genes_orf ref.fa  orfs.tsv.gz acc out.tsv   # ORFs >= 300 nt from gbmitoorf, placed by exact match
  python rdd.py call ref.fa genes.tsv rna.bam [dna.bam]     # caller
Caller: positions with depth >= 20 (BQ >= 20, MAPQ >= 20 already applied); RDD = the top non-reference base
is >= 20% of reads. Positions inside CDS on both strands are dropped. Types are given on the transcript (gene) strand.
DNA filter: an RDD is removed if the DNA run has the same alt at >= 5% (depth >= 10); if DNA depth < 10 it is kept
and counted as 'dna_unchecked'."""
import sys, gzip, csv, subprocess
from collections import Counter, defaultdict

COMP = str.maketrans('ACGTN', 'TGCAN')


def fasta1(path):
    return ''.join(l.strip() for l in open(path) if not l.startswith('>')).upper()


def genes_gb(gb, out):
    from Bio import SeqIO
    rec = SeqIO.read(gb, 'genbank')
    with open(out, 'w') as fh:
        fh.write('gene\tpos1\tstrand\tcodon_pos\n')
        for f in rec.features:
            if f.type != 'CDS':
                continue
            name = f.qualifiers.get('gene', f.qualifiers.get('locus_tag', ['?']))[0]
            idx = 0
            for part in f.location.parts:
                rng = range(int(part.start), int(part.end))
                if part.strand == -1:
                    rng = reversed(rng)
                for p in rng:
                    fh.write(f'{name}\t{p + 1}\t{part.strand}\t{idx % 3 + 1}\n')
                    idx += 1


def genes_orf(fa, orfs, acc, out):
    g = fasta1(fa)
    circ = g + g[:5000]
    rc = circ.translate(COMP)[::-1]
    L = len(g)
    with open(out, 'w') as fh:
        fh.write('gene\tpos1\tstrand\tcodon_pos\n')
        for line in gzip.open(orfs, 'rt'):
            x = line.rstrip('\n').split('\t')
            if x[0] != acc or len(x[-1]) < 300:
                continue
            s = x[-1].upper()
            i = circ.find(s)
            if i >= 0:
                for k in range(len(s)):
                    fh.write(f'orf{x[1]}\t{(i + k) % L + 1}\t1\t{k % 3 + 1}\n')
                continue
            j = rc.find(s)
            if j >= 0:  # rc index j -> forward position len(circ)-1-j going down
                for k in range(len(s)):
                    p = (len(circ) - 1 - (j + k)) % L
                    fh.write(f'orf{x[1]}\t{p + 1}\t-1\t{k % 3 + 1}\n')
            else:
                print('not placed', x[1], file=sys.stderr)


def pileup(fa, bam):
    """{pos1: Counter(base)} from samtools mpileup (BQ>=20)."""
    cmd = ['samtools', 'mpileup', '-B', '-Q', '20', '-q', '20', '-d', '1000000', '-f', fa, bam]
    out = {}
    for line in subprocess.run(cmd, capture_output=True, text=True).stdout.splitlines():
        c, p, ref, d, bases = line.split('\t')[:5]
        cnt, i = Counter(), 0
        while i < len(bases):
            b = bases[i]
            if b == '^':
                i += 2; continue
            if b in '+-':
                j = i + 1
                while j < len(bases) and bases[j].isdigit():
                    j += 1
                i = j + int(bases[i + 1:j]); continue
            if b in '.,':
                cnt[ref.upper()] += 1
            elif b.upper() in 'ACGT':
                cnt[b.upper()] += 1
            i += 1
        out[int(p)] = cnt
    return out


def call(fa, genes, rna, dna=None):
    g = fasta1(fa)
    gene = defaultdict(list)
    for r in csv.DictReader(open(genes), delimiter='\t'):
        gene[int(r['pos1'])].append((int(r['strand']), int(r['codon_pos']), r['gene']))
    pr = pileup(fa, rna)
    pd = pileup(fa, dna) if dna else {}
    tested = Counter()
    types, bycp, unchecked, removed = Counter(), defaultdict(Counter), 0, 0
    sites = []
    for p, cnt in pr.items():
        gs = gene.get(p)
        if not gs or len({s for s, _, _ in gs}) > 1:
            continue
        strand, cp, name = gs[0]
        ref = g[p - 1]
        d = sum(cnt.values())
        if d < 20:
            continue
        tested[cp] += 1
        if (ref if strand == 1 else ref.translate(COMP)) == 'T' and cp == 2:
            tested['U2'] += 1
        alt, n = max(((b, c) for b, c in cnt.items() if b != ref), key=lambda x: x[1], default=('N', 0))
        if n / d < 0.2:
            continue
        if dna:
            dc = pd.get(p, Counter())
            dd = sum(dc.values())
            if dd >= 10 and dc[alt] / dd >= 0.05:
                removed += 1
                continue
            if dd < 10:
                unchecked += 1
        r_, a_ = (ref, alt) if strand == 1 else (ref.translate(COMP), alt.translate(COMP))
        t = f"{r_.replace('T', 'U')}>{a_.replace('T', 'U')}"
        types[t] += 1
        bycp[t][cp] += 1
        sites.append((p, name, strand, cp, t, n, d))
    print(f'# rna {rna} dna {dna}')
    print(f"coding positions tested (depth>=20): {sum(v for k, v in tested.items() if k != 'U2')} by codon pos "
          f"{ {k: v for k, v in tested.items() if k != 'U2'} }; U at codon pos 2: {tested['U2']}")
    print(f'RDD in CDS: {sum(types.values())}; removed by DNA {removed}; dna_unchecked {unchecked}')
    for t, n in types.most_common():
        print(f'  {t}\t{n}\tcodon pos {dict(sorted(bycp[t].items()))}')
    tot = sum(types.values())
    if tot:
        print(f'C>U share {types["C>U"] / tot:.3f}; U>C count {types["U>C"]}')
    with open(rna.replace('.bam', '.rdd.tsv'), 'w') as fh:
        fh.write('pos\tgene\tstrand\tcodon_pos\ttype\talt_n\tdepth\n')
        for s in sorted(sites):
            fh.write('\t'.join(map(str, s)) + '\n')


if __name__ == '__main__':
    {'genes_gb': genes_gb, 'genes_orf': genes_orf, 'call': call}[sys.argv[1]](*sys.argv[2:])
