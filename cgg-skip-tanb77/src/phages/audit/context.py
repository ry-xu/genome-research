# context of CGG sites: next base, position from gene end, +1 open length; T1 skip vs R skip vs T1 non-skip
import sys, csv, collections
rows = list(csv.DictReader(open(sys.argv[1]), delimiter='\t'))
LOW = {'24', '26', '28', '30'}
def show(name, sel):
    nb = collections.Counter(r['next3'][:1] for r in sel)
    pos = collections.Counter(min((int(r['n']) - int(r['k']) - 1) // 4, 3) for r in sel)
    op = sorted(int(r['open']) for r in sel)
    pre = collections.Counter()
    print(f"{name}: n={len(sel)} next base {dict(sorted(nb.items()))}; codons-from-end bins(0-3,4-7,8-11,12-15) {[pos[i] for i in range(4)]}; open median {op[len(op)//2] if op else None}; reached_end {sum(int(r['reached_end']) for r in sel)}")
show('T1 skip', [r for r in rows if r['set'] == 'T1' and r['skip'] == '1' and r['gc_bin'] in LOW])
show('T1 CGG non-skip', [r for r in rows if r['set'] == 'T1' and r['skip'] == '0' and r['gc_bin'] in LOW])
show('T1 CGG last16 non-skip', [r for r in rows if r['set'] == 'T1' and r['skip'] == '0' and r['last16'] == '1' and r['gc_bin'] in LOW])
show('R skip', [r for r in rows if r['set'] == 'R' and r['skip'] == '1' and r['gc_bin'] in LOW])
show('R CGG all', [r for r in rows if r['set'] == 'R' and r['gc_bin'] in LOW])
# T1 CGG not in last 16: how far is the +1 frame open?
sel = [r for r in rows if r['set'] == 'T1' and r['last16'] == '0' and r['gc_bin'] in LOW]
print('T1 internal CGG (not last16): n', len(sel), 'with +1 open>=50:', sum(int(r['open']) >= 50 for r in sel))
sel = [r for r in rows if r['set'] == 'R' and r['last16'] == '0' and r['gc_bin'] in LOW]
print('R internal CGG (not last16): n', len(sel), 'with +1 open>=50:', sum(int(r['open']) >= 50 for r in sel))
