"""Reading frame through each Methanosuratincola insert. Phase from c0: CDS ACO0C9_06155 starts at 0-based 1089196 (+),
insert at 1089656, offset 460 = 1 mod 3, so the codon spanning the insert start is the last upstream base T + insert[:2].
Prints the translation of T+insert+downstream 30 nt in that frame, and of the joined form. usage: msur5_frame.py"""
code = {a + b + c: aa for (a, b, c), aa in zip([(x, y, z) for x in 'TCAG' for y in 'TCAG' for z in 'TCAG'],
        'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG')}
tr = lambda s: ''.join(code.get(s[i:i + 3], 'X') for i in range(0, len(s) - 2, 3))
assert (1089656 - 1089196) % 3 == 1
DN = 'TTAGCCTGGGCTGCAACTTCGCCTGCATCC'
ins = {'GCA_049570475.1 (c0)': 'ACGGATCTAGAGTCCGGCCGAA', 'GCA_038801655.1': 'ACGGATCTAAAGTCCGGCCGAA',
       'GCA_030522375.1 (y0)': 'ACGGGTGAGAGAATCCGGCCGAA', 'GCA_038801705.1': 'ACGGATGAGATCCGGCCGAA',
       'GCA_029907755.1': 'GCGGATGAAAGAGTCCGGCCGAA', 'GCA_030055545.1': 'ACGGATGAGAAAGTCCGGCCGAA',
       'GCA_037478935.1': 'GCGGATAAAAAAAGTCCGGCCGAA', 'GCA_937930465.1': 'GCGGATTCTAAAGTCCGGCCGAA'}
print('joined (any genome):', tr('T' + DN))
for k, v in ins.items():
    p = tr('T' + v + DN)
    print(f'{k}: {len(v)} nt, mod3 {len(v) % 3}, unspliced frame {p}  first stop at aa {p.find("*") + 1 if "*" in p else "-"}')
