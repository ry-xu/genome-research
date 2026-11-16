# Gene calling for the kill1 audit. modes: meta (pyrodigal meta), single (pyrodigal trained per phage, >=20kb),
# orf (plain longest-ORF caller: ATG/GTG/TTG..stop, >=300 nt incl. stop, first start after previous in-frame stop).
# Output TSV: id begin end strand partial_l partial_r ttable  (1-based inclusive coords, as pyrodigal)
import sys, pyrodigal
fa, ids_path, mode, out = sys.argv[1:5]
keep = set(l.split('\t')[0] for l in open(ids_path) if not l.startswith('id\t'))
def read_fa(p):
    name, buf = None, []
    for line in open(p):
        if line.startswith('>'):
            if name: yield name, ''.join(buf)
            name = line[1:].split()[0]; buf = []
        else: buf.append(line.strip())
    if name: yield name, ''.join(buf)
COMP = str.maketrans('ACGTacgtN', 'TGCAtgcaN')
STOPS = {'TAA', 'TAG', 'TGA'}; STARTS = {'ATG', 'GTG', 'TTG'}
def orfs(seq):
    L = len(seq); res = []
    for strand in (1, -1):
        s = seq if strand == 1 else seq.translate(COMP)[::-1]
        for fr in range(3):
            start = None; prev_stop_seen = False
            i = fr
            while i + 3 <= L:
                c = s[i:i+3]
                if c in STOPS:
                    if start is not None and prev_stop_seen and (i + 3 - start) >= 300:
                        res.append((start, i + 3, strand))
                    start = None; prev_stop_seen = True
                elif start is None and c in STARTS and prev_stop_seen:
                    start = i
                i += 3
            # ORFs without an upstream stop or downstream stop are partial: dropped
    out = []
    for a, b, strand in res:
        if strand == 1: out.append((a + 1, b, 1))
        else: out.append((L - b + 1, L - a, -1))
    return out
meta = pyrodigal.GeneFinder(meta=True) if mode == 'meta' else None
n = 0; skipped_short = 0
with open(out, 'w') as fo:
    fo.write('id\tbegin\tend\tstrand\tpartial_l\tpartial_r\tttable\n')
    for name, seq in read_fa(fa):
        if name not in keep or len(seq) < 5000: continue
        seq = seq.upper()
        if mode == 'meta':
            genes = meta.find_genes(seq)
            for g in genes:
                fo.write(f'{name}\t{g.begin}\t{g.end}\t{g.strand}\t{int(g.partial_begin)}\t{int(g.partial_end)}\t{g.translation_table}\n')
        elif mode == 'single':
            if len(seq) < 20000: skipped_short += 1; continue
            gf = pyrodigal.GeneFinder(meta=False)
            gf.train(seq, translation_table=11)
            for g in gf.find_genes(seq):
                fo.write(f'{name}\t{g.begin}\t{g.end}\t{g.strand}\t{int(g.partial_begin)}\t{int(g.partial_end)}\t{g.translation_table}\n')
        else:
            for a, b, strand in orfs(seq):
                fo.write(f'{name}\t{a}\t{b}\t{strand}\t0\t0\t11\n')
        n += 1
        if n % 500 == 0: print(n, flush=True)
print('done phages', n, 'skipped_short(<20kb, single mode)', skipped_short)
