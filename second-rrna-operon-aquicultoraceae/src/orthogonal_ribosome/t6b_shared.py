"""T6-B (RULES_T6.md): is the A. secundus divergent-16S rewrite (bridges + decoding residues) shared by other carrier genera?
usage: uv run python src/t6b_shared.py   (writes out/t6/t6b.txt via stdout)"""
import collections, os, subprocess, sys, tempfile
sys.path.insert(0, os.path.dirname(__file__))
import t4_census as C
from t5_ref import refine, IDX
from t5_count import binom_sf

GEN = {'GCA_036504445.1': 'DASXVE01', 'GCA_902810715.1': 'CADDYS01', 'GCA_027392615.1': 'CADDYS01',
       'GCA_014859275.1': 'Kmv33', 'GCA_035643475.1': 'DASQQM01'}
CLS = {'GCA_027392615.1|DAHYRC010000125.1|2|1404|-': 'norm', 'GCA_027392615.1|DAHYRC010000072.1|9445|11057|-': 'div',
       'GCA_902810715.1|CADDYS010000032.1|17839|19398|+': 'div', 'GCA_035643475.1|DASQQM010000255.1|8358|10072|+': 'div',
       'GCA_036504445.1|DASXVE010000010.1|40680|42216|+': 'norm', 'GCA_036504445.1|DASXVE010000023.1|2051|3668|+': 'div',
       'GCA_014859275.1|JACUUO010000012.1|11422|12967|+': 'div'}  # blastn bitscore classing (T6-B)
CHG_B16 = [702, 783, 784, 1419, 1428, 1430, 1473, 1474, 1475, 1476, 1493, 1495, 1496, 1517]
R = sorted(set(CHG_B16 + [1054, 1492, 1400, 1498, 790, 966, 967, 1518]))
tmp = tempfile.mkdtemp(dir='out/t6/tmp')
seqs = {}
for k, c in CLS.items():
    acc, contig, lo, hi, st = k.split('|')
    g = C.read_fasta(C.ATLAS + IDX[acc])
    s, p0, p1, m = refine(g, (contig, int(lo), int(hi), st), tmp)
    seqs[f'{acc}_{c}'] = s[p0:p1]
    print(f'{acc} {GEN[acc]} {c}: full-CM model {m}, 16S {p1 - p0} nt')
for l in open('out/t5/sec16.fa').read().split('>')[1:]:
    h, *s = l.split('\n'); seqs[f'sec_{h}'] = ''.join(s)
ec = ''.join(l.split('\t')[1] for l in open('out/t6/ecoli16S_num.tsv'))
nums = [int(l.split('\t')[0]) for l in open('out/t6/ecoli16S_num.tsv')]
f = f'{tmp}/all.fa'; open(f, 'w').write(f'>ecoli\n{ec}\n' + ''.join(f'>{k}\n{v}\n' for k, v in seqs.items()))
subprocess.run(['cmalign', '--cpu', '1', '-o', f'{tmp}/all.sto', C.CM, f], check=True, stdout=subprocess.DEVNULL)
aln = collections.defaultdict(str)
for l in open(f'{tmp}/all.sto'):
    if l.strip() and not l.startswith(('#', '//')): n, s = l.split(); aln[n] += s.upper().replace('T', 'U')
col, k = {}, 0
for c, x in enumerate(aln['ecoli']):
    if x in 'ACGU': col[nums[k]] = c; k += 1
b = lambda n, i: aln[n][col[i]] if aln[n][col[i]] in 'ACGU' else None
sd, sn = 'sec_div', 'sec_norm'
G = [i for i in nums if i not in R and b(sd, i) and b(sn, i) and b(sd, i) != b(sn, i)]
assert all(b(sd, i) != b(sn, i) for i in R), 'R must be changed in A. secundus'
# gate
gate = [(n, i) for n in aln if n.endswith('_norm') and n != sn for i in R if b(n, i)]
g_ok = sum(b(n, i) == b(sn, i) for n, i in gate)
print(f'gate: other normal copies equal A. secundus norm at R in {g_ok}/{len(gate)} ({g_ok / max(1, len(gate)):.2f}; bar 0.90)')
print('\nper-position (E. coli no.: ecoli / sec norm / sec div | other divs in order', [n for n in aln if n.endswith('_div') and n != sd], '| other norms)')
divs = [n for n in aln if n.endswith('_div') and n != sd]; norms = [n for n in aln if n.endswith('_norm') and n != sn]
for i in R:
    print(f'  {i}: {b("ecoli", i)} / {b(sn, i)} / {b(sd, i)} | ' + ' '.join(b(n, i) or '-' for n in divs) + ' | ' + ' '.join(b(n, i) or '-' for n in norms))
res = {}
for n in divs:
    rr = [i for i in R if b(n, i)]; gg = [i for i in G if b(n, i)]
    sR = sum(b(n, i) == b(sd, i) for i in rr) / len(rr); sG = sum(b(n, i) == b(sd, i) for i in gg) / len(gg)
    x = sum(b(n, i) == b(sd, i) for i in rr); p = binom_sf(x, len(rr), sG)
    allc = [c for c in range(len(aln[n])) if aln[n][c] in 'ACGU' and aln[sd][c] in 'ACGU']
    ident = sum(aln[n][c] == aln[sd][c] for c in allc) / len(allc)
    canon = sum(b(n, i) == b('ecoli', i) for i in rr) / len(rr)
    acc = n.rsplit('_', 1)[0]; res[n] = (GEN[acc], sR, sG)
    print(f'{n} {GEN[acc]}: s_R {sR:.2f} ({x}/{len(rr)}), s_G {sG:.2f} ({len(gg)} cols), P(s_R>s_G) {p:.3g}; '
          f'E. coli state at R {canon:.2f}; div-div identity {ident:.3f}')
gen = collections.defaultdict(list)
for n, (g, sR, sG) in res.items(): gen[g].append((sR, sG))
shared = sum(all(sR >= 0.8 and sR >= sG for sR, sG in v) for v in gen.values())
spec = sum(all(sR <= 0.4 for sR, sG in v) for v in gen.values())
v = 'GATE FAIL' if g_ok / max(1, len(gate)) < 0.9 else 'SHARED REWRITE' if shared >= 3 else 'A. SECUNDUS-SPECIFIC' if spec >= 3 else 'MIXED'
print(f'genera passing shared {shared}/{len(gen)}, specific {spec}/{len(gen)}\nVERDICT (T6-B): {v}')
