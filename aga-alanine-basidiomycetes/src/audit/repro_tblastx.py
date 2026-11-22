"""Independent placement test, no profile HMMs, no gene tables: tblastx (code 4 both sides) of each query mitogenome
against 3-4 Arg-reading relatives of the same class; for every query codon covered by an HSP (best HSP per relative),
record the relative's residue. Tally relatives' residues at query AGA vs CGN vs GCN codons; decoy: every codon.
usage: python -I repro_tblastx.py genome_dir out_prefix"""
import sys, os, subprocess, csv
from collections import Counter, defaultdict
from Bio.Seq import Seq
from Bio.Data import CodonTable
T4 = CodonTable.unambiguous_dna_by_id[4].forward_table
gdir, outp = sys.argv[1], sys.argv[2]
TREM = ['NC_004336.1', 'BK073231.1', 'BK073612.1', 'BK074202.1']      # Cryptococcus, Kwoniella, Naganishia, Solicoccozyma
EXO = ['CP096881.1', 'BK073411.1', 'BK072594.1', 'NC_023248.1']       # Exobasidium, Meira, Ceraceosorus, Jaminaea
PUC = ['CM018325.1', 'NC_036340.1', 'BK073353.1', 'BK074216.1']       # Symmetrospora, Rhodotorula, Leucosporidium, Sporobolomyces
Q = [  # query, label, role, relatives
    ('NC_054273.1', 'Trichosporon inkin', 'in 32', TREM),
    ('MZ439919.1', 'Apiotrichum gamsii', 'in 32', TREM),
    ('BK074411.1', 'Tremella mesenterica', 'in 32', TREM),
    ('NC_009880.1', 'Tilletia indica', 'in 32', EXO),
    ('BK072866.1', 'Cystobasidium benthicum', 'in 32', PUC),
    ('BK072867.1', 'Cystobasidium minutum', 'in 32', PUC),
    ('BK074202.1', 'Solicoccozyma (Tremellomycetes ctrl)', 'control', TREM[:3]),
    ('BK072956.1', 'Exobasidium (Exobasidiomycetes ctrl)', 'control', EXO[1:]),
    ('CM018325.1', 'Symmetrospora (Cystobasidiomycetes ctrl)', 'control', PUC[1:]),
    ('BK074082.1', 'Rhodotorula diobovata (Pucciniomycotina ctrl)', 'control', PUC[1:]),
    ('BK074509.1', 'Vanrija (G3:U70 carrier, not in 32)', 'extra', TREM),
    ('NC_036422.1', 'Tremella fuciformis (G3:U70 carrier, not in 32)', 'extra', TREM),
    ('PP995850.1', 'Sympodiomycopsis (side)', 'extra', EXO),
]
def seq(a): return ''.join(l.strip() for l in open(f'{gdir}/{a}.fa') if l[0] != '>').upper()
CGN = {'CGT', 'CGC', 'CGA', 'CGG'}; GCN = {'GCT', 'GCC', 'GCA', 'GCG'}
summ = open(outp + '.summary.tsv', 'w')
summ.write('query\tlabel\trole\trelatives\tclass\tcodons_aligned\trel_residues(top)\tunanimous_pos\tunanimous_top\n')
dec = open(outp + '.decoy.tsv', 'w')
dec.write('query\tcodon\tT4\tunanimous_pos\tstd_share\ttop\n')
for q, label, role, rels in Q:
    s = seq(q)
    cov = defaultdict(dict)  # (strand, pos) -> rel -> (bits, residue)
    bad = 0; tot = 0
    for r in rels:
        cmd = ['tblastx', '-query', f'{gdir}/{q}.fa', '-subject', f'{gdir}/{r}.fa', '-query_gencode', '4', '-db_gencode', '4',
               '-evalue', '1e-10', '-seg', 'no', '-outfmt', '6 qframe sframe qstart qend sstart send evalue bitscore pident qseq sseq']
        res = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
        for line in res.splitlines():
            qf, sf, qs, qe, ss, se, ev, bits, pid, qa, sa = line.split('\t')
            qf, qs, bits = int(qf), int(qs), float(bits)
            k = 0
            for x, y in zip(qa, sa):
                if x == '-':
                    continue
                if qf > 0:
                    p = qs + 3 * k; cod = s[p - 1:p + 2]; key = ('+', p)
                else:
                    p = qs - 3 * k; cod = str(Seq(s[p - 3:p]).reverse_complement()); key = ('-', p)
                k += 1
                tot += 1
                if T4.get(cod, '*') != x and x != 'X':
                    bad += 1
                if y in '-X*':
                    continue
                if bits > cov[key].get(r, (0, ''))[0]:
                    cov[key][r] = (bits, y, cod)
    print(q, label, 'HSP residues', tot, 'codon/translation mismatches', bad, file=sys.stderr)
    bycod = defaultdict(Counter); unan = defaultdict(Counter); npos = Counter()
    for key, d in cov.items():
        cod = next(iter(d.values()))[2]
        npos[cod] += 1
        for b, y, c in d.values():
            bycod[cod][y] += 1
        ys = {y for b, y, c in d.values()}
        if len(d) >= 2 and len(ys) == 1:
            unan[cod][ys.pop()] += 1
    def agg(codset):
        c1, c2, n = Counter(), Counter(), 0
        for c in codset:
            c1.update(bycod[c]); c2.update(unan[c]); n += npos[c]
        return n, c1, c2
    for name, cs in (('AGA', {'AGA'}), ('AGG', {'AGG'}), ('CGN', CGN), ('GCN', GCN)):
        n, c1, c2 = agg(cs)
        summ.write(f"{q}\t{label}\t{role}\t{','.join(rels)}\t{name}\t{n}\t{dict(c1.most_common(5))}\t{sum(c2.values())}\t{dict(c2.most_common(5))}\n")
    for c in sorted(npos):
        u = unan[c]; n = sum(u.values()); aa = T4.get(c, '*')
        if n:
            dec.write(f"{q}\t{c}\t{aa}\t{n}\t{u.get(aa, 0) / n:.2f}\t{dict(u.most_common(3))}\n")
summ.close(); dec.close()
