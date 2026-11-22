import csv, numpy as np
from collections import Counter
exec(open('trna_parse_check.py').read().split('miss, cls')[0])
n = Counter()
for r in ucu:
    if r['pair3_70'] == 'G:T': continue
    k = 'ala' if r['acc'] in ala else 'arg' if r['acc'] in arg else 'other'
    rd = readings(r['seq'])
    good = [x for x in rd if x[2] == 'G:T' and (x[5] or x[6])]
    n[(k, bool(good))] += 1
    if good: print(k, r['acc'], r['seq'], good, r['stem5'], r['stem3'], r['pair3_70'])
print(dict(n))
