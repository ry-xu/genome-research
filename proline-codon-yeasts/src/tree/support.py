import sys
from Bio import Phylo
t = Phylo.read(sys.argv[1], 'newick')
t.root_with_outgroup('BK072907.1')
names = {l.split('\t')[0]: l.split('\t')[1] for l in open(sys.argv[2]) if not l.startswith('acc')}
wick3 = ['BK074532.1', 'BK074538.1', 'BK074542.1']
m = t.common_ancestor(wick3)
print('3 Wick PRO clade: support', m.confidence, 'size', len(m.get_terminals()), 'branch', m.branch_length)
for target in ['BK074262.1'] + wick3:
    path = t.get_path(target)
    print('==', target, names.get(target))
    for n in reversed(path[:-1][-6:]):
        tips = [c.name for c in n.get_terminals()]
        gen = sorted(set(names.get(x, '?').split()[0] for x in tips))
        print('  ancestor support', n.confidence, 'tips', len(tips), gen[:5])
# genus monophyly
for g in ['Sugiyamaella', 'Wickerhamiella', 'Starmerella']:
    tips = [c.name for c in t.get_terminals() if names.get(c.name, '').startswith(g)]
    a = t.common_ancestor(tips)
    print(g, len(tips), 'mrca tips', len(a.get_terminals()), 'support', a.confidence)
