"""Gene table for arm R from Clover's census inputs: ORFs of gbmitoorf (CDS = nt[:cds_len]) that hit a census
profile (data/mito_profiles.hmm, E < 1e-5, as reads_codons.positions), placed on the mitogenome by exact match.
usage: python census_genes.py ref.fa acc out.tsv   (writes gene, pos1, strand, codon_pos; prints genes kept)"""
import sys, gzip, os
sys.path[:0] = ['src/codes_census', 'src']  # run from the folder root
import pyhmmer
import census_org, scan_org
from codes import translate_ext

NS = os.environ.get('NEW_STOPS', 'new-stops-euk') + '/'  # the lab's census tree (profiles, gbmitoorf ORFs); not included
COMP = str.maketrans('ACGTN', 'TGCAN')
fa, acc, out = sys.argv[1:4]
census_org.load(NS + 'data/mito_profiles.hmm')
genes = next(g for a, g in scan_org.genomes(NS + 'data/gbmitoorf', {acc}))
seqs, nts = [], []
for g in genes:
    nt = ''.join(c if c in 'ACGT' else 'N' for c in g[7])
    prot, _ = translate_ext(nt, g[5], g[3])
    if len(prot) < 30:
        continue
    nts.append((g[0], nt[:g[5]]))
    seqs.append(pyhmmer.easel.TextSequence(name=str(len(seqs)).encode(), sequence=prot))
block = pyhmmer.easel.TextSequenceBlock(seqs).digitize(scan_org.ALPHA)
best = {}
for hits in pyhmmer.hmmsearch(scan_org.HMMS, block, cpus=1, E=1e-5):
    m = hits.query.name
    m = m.decode() if isinstance(m, bytes) else m
    for hit in hits.reported:
        i = int(hit.name.decode() if isinstance(hit.name, bytes) else hit.name)
        if hit.score > best.get(i, (-1e9,))[0]:
            best[i] = (hit.score, m)
g = ''.join(l.strip() for l in open(fa) if not l.startswith('>')).upper()
L = len(g)
circ = g + g[:5000]
rc = circ.translate(COMP)[::-1]
kept = 0
with open(out, 'w') as fh:
    fh.write('gene\tpos1\tstrand\tcodon_pos\n')
    for i, (sc, m) in sorted(best.items()):
        idx, s = nts[i]
        name = f'{m}_orf{idx}'
        j = circ.find(s)
        if j >= 0:
            pos = [(j + k) % L + 1 for k in range(len(s))]; st = 1
        else:
            j = rc.find(s)
            if j < 0:
                print('not placed', name, file=sys.stderr); continue
            pos = [(len(circ) - 1 - (j + k)) % L + 1 for k in range(len(s))]; st = -1
        kept += 1
        for k, p in enumerate(pos):
            fh.write(f'{name}\t{p}\t{st}\t{k % 3 + 1}\n')
print(acc, 'ORFs', len(genes), 'translated>=30aa', len(seqs), 'with profile hit', kept, sorted({m for _, m in best.values()}))
