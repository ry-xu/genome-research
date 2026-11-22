"""T3: tRNAs (ARAGORN) in the 16S-23S spacer and in 1500 nt after the 23S; 5S already from cmsearch (w_5s.tbl) plus a
lenient cmsearch of the 1500-nt downstream region. Output out/t3_content.tsv"""
import os, re, subprocess, sys, tempfile
sys.path.insert(0, 'src')
from common import operons, region
ARAG = '~/dev/tools/bin/aragorn1.2.41'; CMS = '~/dev/tools/bin/cmsearch'
CM5 = 'lab/work/frozen-core/data/rfam/RF00001.cm'
pat = re.compile(r'^\s*\d+\s+tRNA-(\S+)\s+(c?)\[(\d+),(\d+)\]\s+(\d+)\s+\((\w+)\)', re.M)
def trnas(s):
    with tempfile.NamedTemporaryFile('w', suffix='.fa', delete=False) as f: f.write(f'>x\n{s}\n')
    t = subprocess.run([ARAG, '-t', '-gc11', '-l', '-w', f.name], capture_output=True, text=True).stdout
    os.unlink(f.name)
    res = []
    for mm in pat.finditer(t):
        aa, comp, a, b, pos, ac = mm.groups()
        a, b = int(a), int(b)
        end3 = s[b:b + 3] if not comp else ''
        res.append(f"{aa}-{ac.upper()}{'(minus)' if comp else ''}@{a}-{b}:len{b-a+1}:3'{s[b-4:b]}|{end3}")
    return res
def fives(s):
    with tempfile.TemporaryDirectory() as td:
        open(f'{td}/x.fa', 'w').write(f'>x\n{s}\n')
        subprocess.run([CMS, '--cpu', '1', '-E', '10', '--tblout', f'{td}/t', CM5, f'{td}/x.fa'], capture_output=True)
        res = []
        for l in open(f'{td}/t'):
            if l[0] == '#': continue
            x = l.split()
            if x[9] == '+': res.append(f"5S@{x[7]}-{x[8]}:m{x[5]}-{x[6]}:sc{x[14]}:E{x[15]}")
        return res
o = open('out/t3_content.tsv', 'w')
o.write('genome\tcontig\tcls\tspacer_len\tspacer_tRNAs\tdown_len\tdown_tRNAs\tdown_5S\n')
for r in operons():
    if r['l23a'] is None and r['s16b'] is None: continue
    sp = region(r, r['s16b'] + 1, r['l23a'] - 1) if (r['s16b'] and r['l23a']) else None
    dn = ''
    if r['l23b']:
        for L in range(1500, 0, -50):
            dn = region(r, r['l23b'] + 1, r['l23b'] + L)
            if dn: break
    o.write('\t'.join(map(str, [r['genome'], r['contig'], r['cls'], len(sp) if sp is not None else 'NA', ';'.join(trnas(sp)) if sp else '',
                                len(dn) if r['l23b'] else 'NA', ';'.join(trnas(dn)) if dn else '', ';'.join(fives(dn)) if dn else ''])) + '\n')
