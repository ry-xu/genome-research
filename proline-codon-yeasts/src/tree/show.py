import sys
from Bio import Phylo
names = {l.split('\t')[0]: l.split('\t')[1] for l in open(sys.argv[2]) if not l.startswith('acc')}
sp = {}
import gzip
t = Phylo.read(sys.argv[1], 'newick')
t.root_with_outgroup('BK072907.1')
pro = {'BK074532.1', 'BK074538.1', 'BK074542.1', 'BK074262.1'}
for c in t.get_terminals():
    c.name = ('**PRO** ' if c.name in pro else '') + c.name + ' ' + names.get(c.name, '?') + ' ' + sp.get(c.name, '')
Phylo.draw_ascii(t, column_width=90)
