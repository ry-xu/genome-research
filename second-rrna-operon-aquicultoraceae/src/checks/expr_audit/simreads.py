"""Simulate reads from CG loci in the format the pipeline used (fastq-dump -Z without split: R1+R2 concatenated, 302 nt).
usage: simreads.py REF.fa OUT.fq
Sources: the 16 CG loci (div and norm). Fragments of length L in {40,60,90,150,250,400} every 15 nt along the locus,
sense orientation (RNA). R1 = revcomp(fragment) (dUTP stranded), R2 = fragment; each mate is 151 nt, padded with
adapter + random bases when the fragment is shorter than 151. Substitution error 0 or 1%. Read name encodes source."""
import sys, random
CG = ["PFKS", "PFIW", "MNXI", "JAACTV", "PFTP", "PFNG", "PEXG", "PFFU"]
rnd = random.Random(7)
comp = str.maketrans("ACGTN", "TGCAN")
rc = lambda s: s.translate(comp)[::-1]
A1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"; A2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
seqs = {}; n = None
for l in open(sys.argv[1]):
    if l.startswith(">"):
        n = l[1:].split()[0]
        if len(seqs) >= 51: break
        seqs[n] = []
    else: seqs[n].append(l.strip().upper())
def mut(s, e):
    return "".join(rnd.choice("ACGT".replace(c, "")) if rnd.random() < e else c for c in s)
def mate(fr, ad):
    s = fr[:151]
    if len(s) < 151: s = s + ad + "".join(rnd.choice("ACGT") for _ in range(151))
    return s[:151]
out = open(sys.argv[2], "w"); k = 0
for name, parts in seqs.items():
    if not any(t in name.split("|")[1] for t in CG): continue
    s = "".join(parts)
    for L in (40, 60, 90, 150, 250, 400):
        for st in range(0, len(s) - L + 1, 15):
            fr = s[st:st + L]
            for e in (0.0, 0.01):
                f2 = mut(fr, e)
                r = mate(rc(f2), A1) + mate(f2, A2)
                k += 1
                out.write(f"@sim{k}|{name}|{st}|{L}|{e}\n{r}\n+\n{'I' * len(r)}\n")
print("reads", k)
