"""usage: python3 -I getseq.py FASTA CONTIG START END STRAND NAME  (1-based inclusive) -> fasta on stdout"""
import sys
fa, c, s, e, st, nm = sys.argv[1:7]
s, e = int(s), int(e)
seq, on = [], False
for l in open(fa):
    if l.startswith(">"):
        on = l[1:].split()[0].split("|")[-1] == c
        continue
    if on: seq.append(l.strip())
x = "".join(seq).upper()[s-1:e]
if st == "-":
    x = x[::-1].translate(str.maketrans("ACGTN", "TGCAN"))
print(">" + nm)
for i in range(0, len(x), 80): print(x[i:i+80])
