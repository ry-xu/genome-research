"""Shared helpers: genetic code tables, FASTA IO."""
import itertools

BASES = "TCAG"
CODONS = ["".join(c) for c in itertools.product(BASES, repeat=3)]  # TTT, TTC, TTA, TTG, TCT, ...
_STD = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
STD = dict(zip(CODONS, _STD))
T4 = dict(STD)
T4["TGA"] = "W"          # NCBI table 4
T3 = dict(T4)            # NCBI table 3 (yeast mito)
T3["ATA"] = "M"
for n in BASES:
    T3["CT" + n] = "T"

# codons whose meaning differs among candidate yeast-mito codes (tables 3 vs 4; CGA/CGC absent in t3)
CONTESTED = {"CTT", "CTC", "CTA", "CTG", "ATA", "CGA", "CGC"}


def read_fasta(path):
    name, seq, out = None, [], []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    out.append((name, "".join(seq)))
                name, seq = line[1:], []
            else:
                seq.append(line.strip())
    if name is not None:
        out.append((name, "".join(seq)))
    return out


def revcomp(s):
    return s.translate(str.maketrans("ACGTacgtNn", "TGCAtgcaNn"))[::-1]


def codons_of(nt):
    nt = nt.upper().replace("U", "T")
    return [nt[i:i + 3] for i in range(0, len(nt) - len(nt) % 3, 3)]


def translate_codons(cods, table=T4, mask=None, stop="*"):
    out = []
    for c in cods:
        if mask and c in mask:
            out.append("X")
        elif c not in table:
            out.append("X")
        else:
            aa = table[c]
            out.append(stop if aa == "*" else aa)
    return "".join(out)


def write_fasta(path, recs, width=80):
    with open(path, "w") as fh:
        for n, s in recs:
            fh.write(">" + n + "\n")
            for i in range(0, len(s), width):
                fh.write(s[i:i + width] + "\n")
