"""Independent placement test (no profile HMMs, no project gene tables): NCBI-annotated core proteins of one
Arg-reading reference per class are tblastn'd (code 4) against the query genome and against 2-3 more Arg-reading
relatives; every reference position then has a query codon and the residue of each relative (reference included).
Tallies relatives' residues at the query's AGA, AGG, CGN and GCN codons; 'unanimous' = >= 3 residues (reference + relatives) aligned and all
the same residue. Decoy: the same tally for every query codon.
usage: python -I repro_star.py genome_dir protein_dir out_prefix"""
import sys, subprocess, re
from collections import Counter, defaultdict
from Bio.Seq import Seq
from Bio.Data import CodonTable
T4 = CodonTable.unambiguous_dna_by_id[4].forward_table
gdir, pdir, outp = sys.argv[1:4]
CORE = re.compile(r'\[gene=(cox[123]|co[123i]+|coi+|cob|cytb|atp[689]|nad[1-6]l?|nd[1-6]l?|nad4L)\]', re.I)
REF = {'TREM': 'NC_004336.1', 'EXO': 'NC_023248.1', 'PUC': 'NC_036340.1'}
TREM = ['BK073231.1', 'BK073612.1', 'BK074202.1']      # Kwoniella, Naganishia, Solicoccozyma (+ ref Cryptococcus)
EXO = ['CP096881.1', 'BK073411.1', 'BK072594.1']       # Exobasidium, Meira, Ceraceosorus (+ ref Jaminaea)
PUC = ['CM018325.1', 'BK073353.1', 'BK074216.1']       # Symmetrospora, Leucosporidium, Sporobolomyces (+ ref Rhodotorula)
Q = [
    ('NC_054273.1', 'Trichosporon inkin', 'in 32', 'TREM', TREM),
    ('MZ439919.1', 'Apiotrichum gamsii', 'in 32', 'TREM', TREM),
    ('BK074411.1', 'Tremella mesenterica', 'in 32', 'TREM', TREM),
    ('NC_009880.1', 'Tilletia indica', 'in 32', 'EXO', EXO),
    ('BK072866.1', 'Cystobasidium benthicum', 'in 32', 'PUC', PUC),
    ('BK072867.1', 'Cystobasidium minutum', 'in 32', 'PUC', PUC),
    ('BK074202.1', 'Solicoccozyma aeria-group (Tremellomycetes control)', 'control', 'TREM', TREM[:2]),
    ('BK072956.1', 'Exobasidium (Exobasidiomycetes control)', 'control', 'EXO', EXO[1:]),
    ('CM018325.1', 'Symmetrospora (Cystobasidiomycetes control)', 'control', 'PUC', PUC[1:]),
    ('BK074082.1', 'Rhodotorula diobovata (Pucciniomycotina control)', 'control', 'PUC', PUC[1:]),
    ('BK074509.1', 'Vanrija (G3:U70 carrier, not in 32)', 'extra', 'TREM', TREM),
    ('NC_036422.1', 'Tremella fuciformis (G3:U70 carrier, not in 32)', 'extra', 'TREM', TREM),
    ('PP995850.1', 'Sympodiomycopsis (side)', 'extra', 'EXO', EXO),
]
def seq(a): return ''.join(l.strip() for l in open(f'{gdir}/{a}.fa') if l[0] != '>').upper()
def core_faa(ref):
    out = f'{pdir}/{ref}.core.faa'; keep = False; n = 0
    with open(out, 'w') as o:
        for l in open(f'{pdir}/{ref}.faa'):
            if l[0] == '>':
                keep = bool(CORE.search(l)); n += keep
                if keep: o.write(f'>{CORE.search(l).group(1)}_{n}\n')
            elif keep: o.write(l)
    return out, n
cache = {}
def place(faa, g):
    """-> {(prot, refpos): (codon, residue)} from best tblastn HSP per position"""
    if (faa, g) in cache: return cache[(faa, g)]
    s = seq(g)
    cmd = ['tblastn', '-query', faa, '-subject', f'{gdir}/{g}.fa', '-db_gencode', '4', '-seg', 'no', '-evalue', '1e-5',
           '-outfmt', '6 qseqid sframe qstart qend sstart send bitscore qseq sseq']
    best = {}; bad = 0
    for line in subprocess.run(cmd, capture_output=True, text=True, check=True).stdout.splitlines():
        qid, sf, qs, qe, ss, se, bits, qa, sa = line.split('\t')
        sf, qs, ss, bits = int(sf), int(qs), int(ss), float(bits)
        qp, k = qs, 0
        for x, y in zip(qa, sa):
            if y != '-':
                if sf > 0: p = ss + 3 * k; cod = s[p - 1:p + 2]
                else: p = ss - 3 * k; cod = str(Seq(s[p - 3:p]).reverse_complement())
                k += 1
                if x != '-':
                    if T4.get(cod, '*') != y and y not in 'X*': bad += 1
                    if bits > best.get((qid, qp), (0,))[0]: best[(qid, qp)] = (bits, cod, y)
            if x != '-': qp += 1
    cache[(faa, g)] = {k: v[1:] for k, v in best.items()}
    if bad: print('WARNING codon/translation mismatches', g, bad, file=sys.stderr)
    return cache[(faa, g)]
CL = {'AGA': {'AGA'}, 'AGG': {'AGG'}, 'CGN': {'CGT', 'CGC', 'CGA', 'CGG'}, 'GCN': {'GCT', 'GCC', 'GCA', 'GCG'}}
summ = open(outp + '.summary.tsv', 'w'); dec = open(outp + '.decoy.tsv', 'w'); sites = open(outp + '.aga_sites.tsv', 'w')
summ.write('query\tlabel\trole\tref+relatives\tclass\tquery_codons_placed\tunanimous\tunan_top\tall_relative_residues_top\n')
dec.write('query\tcodon\tT4\tunanimous\tstd_share\ttop3\n')
sites.write('query\tprot\trefpos\tquery_codon\tref_residue\trelative_residues\n')
for q, label, role, cls, rels in Q:
    faa, n = core_faa(REF[cls])
    ref = {}
    name = None
    for l in open(faa):
        if l[0] == '>': name = l[1:].strip(); ref[name] = ''
        else: ref[name] += l.strip()
    qp = place(faa, q); rp = [place(faa, r) for r in rels]
    unan = defaultdict(Counter); allr = defaultdict(Counter); npl = Counter()
    for (prot, pos), (cod, aa) in qp.items():
        res = [ref[prot][pos - 1]] + [r[(prot, pos)][1] for r in rp if (prot, pos) in r]
        res = [x for x in res if x not in 'X*-']
        npl[cod] += 1; allr[cod].update(res)
        if len(res) >= 3 and len(set(res)) == 1: unan[cod][res[0]] += 1
        if cod == 'AGA': sites.write(f'{q}\t{prot}\t{pos}\t{cod}\t{ref[prot][pos - 1]}\t{"".join(res[1:])}\n')
    for c, cs in CL.items():
        u = sum((unan[x] for x in cs), Counter()); a = sum((allr[x] for x in cs), Counter())
        summ.write(f"{q}\t{label}\t{role}\t{REF[cls]}+{','.join(rels)}\t{c}\t{sum(npl[x] for x in cs)}\t{sum(u.values())}\t"
                   f"{dict(u.most_common(4))}\t{dict(a.most_common(4))}\n")
    for c in sorted(npl):
        u = unan[c]; t = sum(u.values()); aa = T4.get(c, '*')
        if t: dec.write(f"{q}\t{c}\t{aa}\t{t}\t{u.get(aa, 0) / t:.2f}\t{dict(u.most_common(3))}\n")
summ.close(); dec.close(); sites.close()
