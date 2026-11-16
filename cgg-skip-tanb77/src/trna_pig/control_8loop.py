"""Positive control: does ARAGORN (pyaragorn, codon-hunt's settings) report a tRNA-Gly whose anticodon loop has one
extra base (5'-UCCG-3', able to read CGGA)? Real tRNA-Gly(TCC) of MGYG000302328 (tRNAscan 11.trna1), unchanged and with a
G inserted after anticodon base 36 (sequence from the .ss, upper-cased), each embedded in 3 kb of random sequence (GC 45%, seed 1), both strands.
Writes control_8loop.fna for tRNAscan-SE too."""
import random

import pyaragorn

wt = "GCAGGTGTCGTATAATGGCTATTACCTCAGCCTTCCAAGCTGATGACGAGGGTTCGATTCCCTCTACCTGCTCCA"
i = wt.index("CTTCCAA") + 5  # after C36
mut = wt[:i] + "G" + wt[i:]
random.seed(1)
flank = lambda: "".join(random.choices("ACGT", weights=[55, 45, 45, 55], k=1500))
rc = lambda s: s[::-1].translate(str.maketrans("ACGT", "TGCA"))
f = pyaragorn.RNAFinder(11, tmrna=False, linear=True)
with open("scratch/trnascan_pig/control_8loop.fna", "w") as out:
    for name, t in (("wt", wt), ("ins_G", mut), ("ins_G_rc", rc(mut))):
        s = flank() + t + flank()
        out.write(f">{name}\n{s}\n")
        hits = [(g.begin, g.end, g.amino_acids, g.anticodon, g.anticodon_length, round(g.energy, 1))
                for g in f.find_rna(s)]
        print(name, t[i - 8:i + 8], hits)
