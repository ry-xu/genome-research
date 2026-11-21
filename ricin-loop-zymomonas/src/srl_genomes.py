"""Count SRL tetraloop contexts (AGTAC NNNN GGACC, both strands) in NCBI Datasets genome zips."""
import sys, zipfile, re, collections

RX = re.compile("AGTAC([ACGT]{4})GGACC")
COMP = str.maketrans("ACGT", "TGCA")


def seqs(zp):
    with zipfile.ZipFile(zp) as z:
        for n in z.namelist():
            if n.endswith((".fna", ".fasta")):
                cur = []
                for line in z.read(n).decode().splitlines():
                    if line.startswith(">"):
                        if cur:
                            yield "".join(cur)
                        cur = []
                    else:
                        cur.append(line.strip().upper())
                if cur:
                    yield "".join(cur)


for zp in sys.argv[1:]:
    c = collections.Counter()
    for s in seqs(zp):
        for strand in (s, s.translate(COMP)[::-1]):
            c.update(m.group(1) for m in RX.finditer(strand))
    acc = zp.split("/")[-1].removesuffix(".zip")
    print(acc, sum(c.values()), " ".join(f"{k}:{v}" for k, v in sorted(c.items())), sep="\t")
