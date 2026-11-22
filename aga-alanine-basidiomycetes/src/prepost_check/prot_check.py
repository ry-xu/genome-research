"""Independent tally of set A Sage PSMs; database composition; expected false 'A-only' PSMs; AGA opportunity in seen proteins."""
import csv, re, sys
from collections import Counter, defaultdict
RES = sys.argv[1]; DB = 'lab/work/new-stops-euk/results/aga_ala/prot/db.fasta'
seqs = {}; name = None
for line in open(DB):
    if line.startswith('>'): name = line[1:].strip(); seqs[name] = []
    else: seqs[name].append(line.strip())
seqs = {k: ''.join(v) for k, v in seqs.items()}
mito = {k: v for k, v in seqs.items() if k.startswith('mito|')}
nuc = {k: v for k, v in seqs.items() if k.startswith('nuc|')}
print('db: mito entries', len(mito), 'residues', sum(map(len, mito.values())), '| nuclear entries', len(nuc), 'residues', sum(map(len, nuc.values())))
orfs = Counter(k.split('|')[1].rsplit('_', 2)[0] for k in mito)
print('mito ORFs per genome (x3 versions):', {k: v // 3 for k, v in orfs.items()})
# PSM tally
tot = Counter(); cls = Counter(); classes = defaultdict(Counter); peps = defaultdict(set)
for r in csv.DictReader(open(RES), delimiter='\t'):
    q = float(r['spectrum_q'])
    if q > 0.01: continue
    dec = r['label'] == '-1'
    prots = r['proteins'].split(';')
    m = [p for p in prots if p.replace('rev_', '').startswith('mito|')]
    tot['decoy' if dec else 'target'] += 1
    if not m: continue
    v = ''.join(sorted({p.rsplit('|', 1)[1] for p in m}, key='RAC'.index)) + ('+nuc' if len(m) < len(prots) else '')
    classes[v]['decoy' if dec else 'target'] += 1
    if not dec: peps[v].add(re.sub(r'\[[^\]]*\]-?', '', r['peptide']))
print('all PSMs at spectrum q<=0.01:', dict(tot))
for v, c in classes.items(): print('  mito class', v, dict(c), sorted(peps[v]))
# Of all decoy PSMs at any q, what share land on mito entries? (proxy for where false matches fall)
dm = Counter()
for r in csv.DictReader(open(RES), delimiter='\t'):
    if r['label'] != '-1': continue
    prots = r['proteins'].split(';')
    m = [p for p in prots if p.replace('rev_', '').startswith('mito|')]
    v = ''.join(sorted({p.rsplit('|', 1)[1] for p in m}, key='RAC'.index)) if m else 'nuc'
    dm[v] += 1
n = sum(dm.values()); print('decoy PSMs (all q) by class:', {k: f'{c} ({c/n:.4f})' for k, c in dm.most_common()})
# expected false target PSMs at q<=0.01 in class A ~ (decoys passing) * share of decoys in class A
print('expected false A-only PSMs at q<=0.01 ~', tot['decoy'], 'x', f"{dm.get('A',0)/n:.4f}", '=', round(tot['decoy'] * dm.get('A', 0) / n, 3))
# AGA opportunity: AGA positions in the cox2 ORF and in the peptide's ORF
sites = defaultdict(list)
for r in csv.DictReader(open('lab/work/new-stops-euk/results/aga_ala/prot/sites.tsv'), delimiter='\t'):
    if r['codon'] == 'AGA': sites[r['orf']].append(int(r['pos']))
for o in ('CM001777.1_+2_24401', 'CM001777.1_+1_24946'):
    s = mito[f'mito|{o}|R']
    print(o, 'len', len(s), 'AGA count', len(sites[o]), 'positions', sites[o][:20])
total_aga = sum(len(v) for k, v in sites.items() if k.startswith('CM001777.1'))
print('AGA codons in all CM001777.1 ORFs >=50 codons:', total_aga)
