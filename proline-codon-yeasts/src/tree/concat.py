"""Concatenated mito protein alignment for a set of gbmitoorf genomes: genes -> table-4 translation with CTN masked as X ->
best profile per gene (E<1e-5) -> hmmalign per profile, match columns only -> concatenated FASTA (missing genes as gaps)."""
import sys, os, collections
sys.path.insert(0, 'src'); sys.path.insert(0, 'src/codes_census')
import pyhmmer, scan_org
from codes import translate_ext
prefix, hmmf, accf, out = sys.argv[1:5]
accs = [l.split()[0] for l in open(accf) if l.strip()]
with pyhmmer.plan7.HMMFile(hmmf) as fh: hmms = [h for h in fh if (h.name.decode() if isinstance(h.name, bytes) else h.name).startswith('mito_Fun_')]
alpha = pyhmmer.easel.Alphabet.amino()
per = collections.defaultdict(dict)  # model -> acc -> seq
for acc, genes in scan_org.genomes(prefix, set(accs)):
    seqs = []
    for g in genes:
        nt = ''.join(c if c in 'ACGT' else 'N' for c in g[7])
        prot, _ = translate_ext(nt, g[5], g[3])
        prot = ''.join('X' if nt[3*i:3*i+2] == 'CT' else a for i, a in enumerate(prot))
        if len(prot) >= 30: seqs.append(pyhmmer.easel.TextSequence(name=str(len(seqs)).encode(), sequence=prot))
    blk = pyhmmer.easel.TextSequenceBlock(seqs).digitize(alpha)
    best = {}
    for hits in pyhmmer.hmmsearch(hmms, blk, cpus=2, E=1e-5):
        m = hits.query.name; m = m.decode() if isinstance(m, bytes) else m
        for h in hits.reported:
            if h.score > best.get(m, (-1e9, None))[0]: best[m] = (h.score, int(h.name))
    for m, (_, i) in best.items(): per[m][acc] = seqs[i]
cols = {}
for h in hmms:
    m = h.name.decode() if isinstance(h.name, bytes) else h.name
    if len(per[m]) < 0.8 * len(accs): continue
    names = list(per[m]); blk = pyhmmer.easel.TextSequenceBlock([per[m][a] for a in names]).digitize(alpha)
    msa = pyhmmer.hmmer.hmmalign(h, blk, trim=True)
    rows = {}
    for a, s in zip(names, msa.alignment):
        rows[a] = ''.join(c for c in s if not c.islower() and c != '.')
    L = len(next(iter(rows.values())))
    cols[m] = {a: rows.get(a, '-' * L) for a in accs}
    print(m, len(per[m]), L, file=sys.stderr)
with open(out, 'w') as fh:
    for a in accs: fh.write(f'>{a}\n' + ''.join(cols[m][a] for m in cols) + '\n')
