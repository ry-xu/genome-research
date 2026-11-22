"""One-sided Fisher: Leu-consensus share at Choristocarpus ATA (all genes) vs at Met-codon ATG in the 5 NEG controls."""
from math import comb
def fisher_greater(a, b, c, d):  # P(X >= a), rows: (a,b) test, (c,d) control
    n1, n2, k = a + b, c + d, a + c; N = n1 + n2
    return sum(comb(n1, x) * comb(n2, k - x) for x in range(a, min(n1, k) + 1)) / comb(N, k)
ata_L, ata_n = 18 + 12, 108 + 38          # cons.out POS ATA resp + ribo
atg_L, atg_n = 47 + 15, 627 + 96          # cons.out NEG ATG resp + ribo
cat_L, cat_n = 2 + 1, 35 + 5              # POS ATG
print(f"POS ATA L {ata_L}/{ata_n} = {ata_L/ata_n:.3f}; NEG ATG L {atg_L}/{atg_n} = {atg_L/atg_n:.3f}; POS ATG L {cat_L}/{cat_n} = {cat_L/cat_n:.3f}")
print("one-sided Fisher ATA vs NEG ATG:", f"{fisher_greater(ata_L, ata_n-ata_L, atg_L, atg_n-atg_L):.2g}")
print("one-sided Fisher ATA vs POS ATG:", f"{fisher_greater(ata_L, ata_n-ata_L, cat_L, cat_n-cat_L):.2g}")
