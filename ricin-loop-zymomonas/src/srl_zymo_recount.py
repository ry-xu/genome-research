"""SRL tetraloop per 23S copy: regex AGTAC....GGACC on both strands (E. coli 2654-2667 AGUACGAGAGGACC; loop = 2659-2662).
Also an anchor upstream to confirm it is 23S: report 30 nt context."""
import gzip, re, sys, glob, collections
comp = str.maketrans("ACGTN", "TGCAN")
pat = re.compile(r"(?=(CTGCTCCTAGTAC(....)GGACCG))")  # E. coli 2647-2667 region (GCUCCUAGUACGAGAGGACCG)
loose = re.compile(r"(?=(AGTAC(....)GGACCG))")
for f in sorted(glob.glob(sys.argv[1] + "/*.fna.gz")):
    seq = {}; name = None
    for l in gzip.open(f, "rt"):
        if l[0] == ">": name = l[1:].split()[0]; seq[name] = []
        else: seq[name].append(l.strip().upper())
    c = collections.Counter(); cl = collections.Counter()
    for n, s in seq.items():
        s = "".join(s)
        for st in (s, s.translate(comp)[::-1]):
            for m in pat.finditer(st): c[m.group(2)] += 1
            for m in loose.finditer(st): cl[m.group(2)] += 1
    print(f.split("/")[-1][:34], "strict:", dict(c), "loose:", dict(cl))
