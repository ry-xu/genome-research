import sys
from Bio import SeqIO
from Bio.Seq import reverse_complement
a = str(SeqIO.read(sys.argv[1], 'fasta').seq).upper()
b = str(SeqIO.read(sys.argv[2], 'genbank').seq).upper()
print(len(a), len(b), 'identical' if a == b else 'differ')
if a != b:
    # rotation / reverse complement check
    for name, s in (('fwd', b), ('rc', reverse_complement(b))):
        i = (s + s).find(a[:200])
        print(name, 'rotation offset of OR start in BK:', i)
        if i >= 0:
            r = (s + s)[i:i + len(s)]
            print('  identical after rotation:', r == a, 'mismatches', sum(x != y for x, y in zip(r, a)))
