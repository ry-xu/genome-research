import sys
from Bio import SeqIO
recs = {r.id.split('.')[0]: str(r.seq).upper() for r in SeqIO.parse(sys.argv[1], 'genbank')}
a, b = recs['BK074532'], recs['OR703884']
print(len(a), len(b))
from Bio.Seq import reverse_complement
for name, s in (('fwd', b), ('rc', reverse_complement(b))):
    i = (a + a).find(s)
    print(name, 'rotation offset', i)
# codon usage CTC/CTG in all six-frame ORFs >= 100 codons (table 4 stops TAA/TAG) for the four
import re
for acc in ['BK074262','BK074532','BK074538','BK074542']:
    s = recs[acc]
    tot = {'CTC':0,'CTG':0,'CTT':0,'CTA':0}
    for seq in (s, reverse_complement(s)):
        sx = seq + seq[:3000]
        for fr in range(3):
            cods = [sx[i:i+3] for i in range(fr, len(sx)-2, 3)]
            orf = []
            for c in cods:
                if c in ('TAA','TAG'):
                    if len(orf) >= 100:
                        for x in orf:
                            if x in tot: tot[x] += 1
                    orf = []
                else: orf.append(c)
    print(acc, 'ORF>=100 codons, 6 frames:', tot)
