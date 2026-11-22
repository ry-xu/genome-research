import sys, gzip, csv, re, random, collections
inp, out_fa, out_tsv = sys.argv[1], sys.argv[2], sys.argv[3]
def fam(pn, gn):
    p = pn.lower()
    g = gn.lower().split()
    rules = [
        (r'^cytochrome c oxidase subunit 1\b', 'cox1'), (r'^cytochrome c oxidase subunit 2\b', 'cox2'),
        (r'^cytochrome c oxidase subunit 3\b', 'cox3'), (r'^cytochrome b\b(?! ?5)(?!-c)', 'cob'),
        (r'^atp synthase subunit a\b', 'atp6'), (r'^atp synthase (protein|subunit) 8\b', 'atp8'),
        (r'^atp synthase subunit (c|9)\b', 'atp9'), (r'^atp synthase subunit alpha', 'atp1'),
        (r'^nadh-ubiquinone oxidoreductase chain 4l', 'nad4L'),
        (r'^nadh-ubiquinone oxidoreductase chain (\d+)', None),
        (r'^nadh dehydrogenase \[ubiquinone\] iron-sulfur protein', None),
        (r'(\d+)s ribosomal protein ([sl]\d+)', None), (r'ribosomal protein ([sl]\d+)', None),
    ]
    for pat, name in rules:
        m = re.search(pat, p)
        if m:
            if name: return name
            if 'chain' in pat: return 'nad' + m.group(1)
            if 'iron-sulfur' in pat: return None
            return 'r' + m.group(m.lastindex)
    return None
def group(lin):
    for k in ['Vertebrata', 'Arthropoda', 'Fungi', 'Viridiplantae']:
        if k + ' (' in lin: return k
    if 'Metazoa (' in lin: return 'Metazoa_other'
    return 'Euk_other'
cap = {'Vertebrata': 3, 'Arthropoda': 3, 'Metazoa_other': 6, 'Fungi': 10, 'Viridiplantae': 8, 'Euk_other': 12}
rows = collections.defaultdict(list)
with gzip.open(inp, 'rt') as f:
    r = csv.DictReader(f, delimiter='\t')
    for row in r:
        if 'Eukaryota' not in row['Taxonomic lineage']: continue
        fm = fam(row['Protein names'], row['Gene Names'])
        if fm is None: continue
        if 'fragment' in row['Protein names'].lower(): continue
        if int(row['Length']) < 40: continue
        rows[(fm, group(row['Taxonomic lineage']))].append(row)
rng = random.Random(7)
sel = []
for (fm, gr), lst in sorted(rows.items()):
    lst.sort(key=lambda x: x['Entry']); rng.shuffle(lst)
    seen = set(); n = 0
    for row in lst:
        genus = row['Organism'].split()[0]
        if genus in seen: continue
        seen.add(genus); sel.append((fm, gr, row)); n += 1
        if n >= cap[gr]: break
with open(out_fa, 'w') as fo, open(out_tsv, 'w') as ft:
    for fm, gr, row in sel:
        fo.write('>%s|%s|%s\n%s\n' % (row['Entry'], fm, gr, row['Sequence']))
        ft.write('%s\t%s\t%s\t%s\n' % (row['Entry'], fm, gr, row['Organism']))
c = collections.Counter((fm, gr) for fm, gr, _ in sel)
fams = sorted(set(f for f, g in c))
print(len(sel), 'selected')
for fm in fams:
    print(fm, {g: c[(fm, g)] for g in cap if c[(fm, g)]})
