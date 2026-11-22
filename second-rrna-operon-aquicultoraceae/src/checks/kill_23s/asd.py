"""For each 16S hit (cmsearch tblout), report whether the 3' end region (last 60 nt of the hit + 40 nt downstream) contains
the anti-SD core CCUCC and the 3'-terminal GAUCACC motif."""
import sys
seqs, n = {}, None
for l in open("ctg16.fa"):
    if l.startswith(">"): n = l[1:].strip(); seqs[n] = []
    else: seqs[n].append(l.strip())
seqs = {k: "".join(v) for k, v in seqs.items()}
rc = lambda s: s[::-1].translate(str.maketrans("ACGTN", "TGCAN"))
for l in open("cms_ctg16.tbl"):
    if l.startswith("#"): continue
    f = l.split(); c, mf, mt, a, b, st, sc = f[0], int(f[5]), int(f[6]), int(f[7]), int(f[8]), f[9], f[14]
    s = seqs[c]
    if st == "+": w = s[b - 60:b + 40]
    else: w = rc(s[b - 41:b + 59])
    print(c, a, b, st, f"model {mf}-{mt}", "score", sc, "CCTCC" in w, "GATCACC" in w, "GTAACAAGG" in (s[a-1:b] if st=='+' else rc(s[b-1:a])))
