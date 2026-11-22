"""Red algal genome counts per census file, to locate the post's '319 red algal mitogenomes' (audit)."""
import sys, numpy as np
from collections import Counter
for p in sys.argv[1:]:
    d = np.load(p)
    tax = d['taxonomy'].tolist(); acc = d['acc'].tolist()
    red = [i for i, t in enumerate(tax) if 'Rhodophyta' in t]
    sg = d['site_genome']
    n_with = len(set(np.unique(sg).tolist()) & set(red))
    print(p.split('/')[-1], 'genomes', len(acc), 'red', len(red), 'red with sites', n_with,
          'classes', dict(Counter(t.split(';')[2] if len(t.split(';')) > 2 else '?' for t in (tax[i] for i in red)).most_common(6)))
