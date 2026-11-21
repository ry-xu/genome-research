"""(d) hyper-editing rescue (Porath et al. 2014 style). Writes a converted reference and converted reads.
usage: convert.py ref.fa reads.fq X Y outprefix
  ref: + strand X->Y. reads: set S = X->Y (sense to + strand), set A = comp(X)->comp(Y) (antisense); the original
  read sequence is kept in the read name so the original mismatches can be counted after mapping."""
import sys, gzip
fa, fq, X, Y, out = sys.argv[1:6]
C = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
g = ''.join(l.strip() for l in open(fa) if not l.startswith('>')).upper()
open(out + '.ref.fa', 'w').write('>conv\n' + g.replace(X, Y) + '\n')
op = gzip.open if fq.endswith('.gz') else open
with op(fq, 'rt') as fh, open(out + '.S.fq', 'w') as s, open(out + '.A.fq', 'w') as a:
    while True:
        h = fh.readline()
        if not h: break
        seq = fh.readline().strip(); fh.readline(); qual = fh.readline().strip()
        name = h[1:].split()[0]
        s.write(f'@{name}|{seq}\n{seq.replace(X, Y)}\n+\n{qual}\n')
        a.write(f'@{name}|{seq}\n{seq.replace(C[X], C[Y])}\n+\n{qual}\n')
