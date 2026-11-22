"""Write the read pairs (and singletons) whose names are in NAMES from an interleaved fastq-dump --split-spot -Z stream.

usage: fastq-dump --split-spot -Z RUN.sra | python3 -I pick.py NAMES OUT12.fq OUTSINGLE.fq
(rewritten 15 Nov after the scratchpad copy was lost; same behaviour as the Rule 2 pick.py)
"""
import sys

names = set(open(sys.argv[1]).read().split())
o12, os_ = open(sys.argv[2], "w"), open(sys.argv[3], "w")


def flush(b):
    if b and b[0][0][1:].split()[0] in names:
        (o12 if len(b) == 2 else os_).write("".join("".join(r) for r in b))


buf = []
it = iter(sys.stdin)
for h in it:
    rec = (h, next(it), next(it), next(it))
    n = h[1:].split()[0]
    if buf and buf[0][0][1:].split()[0] != n:
        flush(buf)
        buf = []
    buf.append(rec)
flush(buf)
