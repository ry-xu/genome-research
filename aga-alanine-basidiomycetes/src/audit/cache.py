"""Load the 3 mito census files with the project's own sitescan_org.load() and cache fungal sites + global column
consensus to audit/cache.npz (so later audit scripts don't reload 80M sites)."""
import sys, numpy as np
P = 'lab/work/new-stops-euk'
sys.path.insert(0, P + '/src/codes_census'); sys.path.insert(0, P + '/src')
from sitescan_org import load
paths = [P + f'/data/codes/{x}.census.npz' for x in ('mito', 'gbmito', 'gbmitoorf')]
D = load(paths)
tax = np.array(D['tax'])
fungi = np.array([';Fungi;' in t + ';' or t.startswith('Eukaryota;Fungi') for t in D['tax']])
m = fungi[D['sg']]
src = np.concatenate([np.full(len(np.load(p)['acc']), i) for i, p in enumerate(paths)])
from collections import Counter
fam = Counter(t.rsplit(';', 1)[0] for t in D['tax'])
w = np.array([1 / fam[t.rsplit(';', 1)[0]] for t in D['tax']])
np.savez(P + '/results/aga_ala/audit/cache.npz', acc=np.array(D['acc']), tax=tax, table=np.array(D['table']), src=src,
         fungi=fungi, w=w, sg=D['sg'][m], sc=D['sc'][m], col=D['col'][m], res=D['res'][m], strict=D['strict'],
         cons=D['cons'], freq=D['freq'], tot=D['tot'], models=np.array(D['models']), K=D['K'])
print(len(D['acc']), 'genomes', fungi.sum(), 'fungal', m.sum(), 'fungal sites', D['strict'].sum(), 'strict cols')
