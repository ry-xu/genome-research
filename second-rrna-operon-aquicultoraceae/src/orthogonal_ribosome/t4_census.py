"""T4 census of 16S anti-SD tail heterogeneity (RULES_T4.md, frozen). Per genome: cmsearch --hmmonly RF00177, copies
spanning >= 1,200 model positions, 3' regions, majority tail columns (aligned to E. coli rrsB 1530-1542), per-copy
class. Usage:
  t4_census.py scan LIST OUT.jsonl          (resumable: skips accessions already in OUT)
  t4_census.py one FASTA                    (a local FASTA, e.g. a planted copy; prints the record)
"""
import gzip, json, os, subprocess, sys, tempfile
from collections import Counter
from Bio import Align

ATLAS = 'data/external/codon-atlas/'
CM = 'data/external/rfam/RF00177.cm'
CMSEARCH = 'cmsearch'
ECO = open('data/external/ssu/ecoli_rrsB.fna').read().split('\n', 1)[1].replace('\n', '')
assert ECO.endswith('GGATCACCTCCTTA'), ECO[-20:]
ECO_REF = ECO[1439:1542]          # E. coli 1440-1542, tail = last 13 (1530-1542)
TAIL0 = len(ECO_REF) - 13
COMP = str.maketrans('ACGTN', 'TGCAN')

AL = Align.PairwiseAligner(mode='global', match_score=2, mismatch_score=-3, open_gap_score=-5, extend_gap_score=-2)
AL_ECO = Align.PairwiseAligner(mode='global', match_score=2, mismatch_score=-3, open_gap_score=-5, extend_gap_score=-2,
                               end_gap_score=0)  # free end gaps: the 3' region overhangs the E. coli reference


def rc(s):
    return s.translate(COMP)[::-1]


def has_core(t):
    return 'CCTC' in t or 'CTCC' in t


def read_fasta(path):
    op = gzip.open if path.endswith('.gz') else open
    seqs, n = {}, None
    with op(path, 'rt') as f:
        for l in f:
            if l[0] == '>':
                n = l[1:].split()[0]; seqs[n] = []
            else:
                seqs[n].append(l.strip().upper())
    return {k: ''.join(v) for k, v in seqs.items()}


def copies(fasta_plain):
    """[(contig, start, end, strand, mdl_from, mdl_to)] 1-based, merged."""
    with tempfile.NamedTemporaryFile(suffix='.tbl', delete=False) as t:
        tbl = t.name
    subprocess.run([CMSEARCH, '--cpu', '1', '--hmmonly', '--noali', '-E', '1e-10', '--tblout', tbl, CM, fasta_plain],
                   check=True, stdout=subprocess.DEVNULL)
    hits = []
    for l in open(tbl):
        if l[0] == '#': continue
        f = l.split()
        mf, mt, a, b, st = int(f[5]), int(f[6]), int(f[7]), int(f[8]), f[9]
        if mt - mf + 1 < 1200: continue
        lo, hi = min(a, b), max(a, b)
        hits.append([f[0], lo, hi, st, mf, mt, float(f[14])])
    os.unlink(tbl)
    hits.sort(key=lambda h: -h[6])
    keep = []
    for h in hits:
        if any(k[0] == h[0] and k[3] == h[3] and min(k[2], h[2]) - max(k[1], h[1]) + 1 > 0.5 * min(k[2] - k[1], h[2] - h[1]) for k in keep):
            continue
        keep.append(h)
    return sorted(keep, key=lambda h: (h[0], h[1]))


def region(seqs, h):
    c, lo, hi, st, mf, mt = h[:6]
    push = max(0, 1533 - mt)
    s = seqs[c]
    if st == '+':
        e = hi + push; return s[max(0, e - 90):e + 40]
    e = lo - push; return rc(s[max(0, e - 41):e + 89])


def cols_in_target(aln, target_cols):
    """residues of the query aligned to the given target columns (0-based); gaps skipped."""
    out = []
    t, q = aln.indices
    for k in range(t.size):
        if t[k] in target_cols and q[k] >= 0:
            out.append(aln.query[q[k]])
    return ''.join(out)


def tail_cols(maj):
    """indices in maj aligned to E. coli 1530-1542."""
    a = AL_ECO.align(ECO_REF, maj)[0]
    t, q = a.indices
    return [int(q[k]) for k in range(t.size) if t[k] >= TAIL0 and q[k] >= 0]


def only_hp_indel(a, b):
    """True if a and b differ only by one 1-nt indel inside a homopolymer of >= 4."""
    if abs(len(a) - len(b)) != 1: return False
    lo, hi = (a, b) if len(a) < len(b) else (b, a)
    for i in range(len(hi)):
        if hi[:i] + hi[i + 1:] == lo:
            c = hi[i]; j = i
            while j > 0 and hi[j - 1] == c: j -= 1
            k = i
            while k + 1 < len(hi) and hi[k + 1] == c: k += 1
            if k - j + 1 >= 4: return True
    return False


def hp_shift(a, b):
    """T4-A7: regions are fixed-length windows, so a 1-nt indel shifts one window end by 1 nt; test with one end trimmed."""
    return any(only_hp_indel(x, y) for x, y in ((a, b[:-1]), (a, b[1:]), (a[:-1], b), (a[1:], b)))


def cm_end(seqs, h, tmp, fl=300):
    """T4-A7: for a hit ending before model 1530, re-anchor its 3' end at the full-CM hit end (handles insertions near
    helix 45). Returns (status, hit); status ok | breaks (CM end < 1520, or < 0.8 nt per model position = truncated) | edge."""
    c, lo, hi, st, mf, mt, sc = h
    s = seqs[c]; a, b = lo - 1 - fl, hi + fl
    if a < 0 or b > len(s): return 'edge', h
    loc = s[a:b] if st == '+' else rc(s[a:b])
    end_hmm = (hi - a) if st == '+' else (b - lo + 1)
    f = os.path.join(tmp, 'loc.fna'); t = os.path.join(tmp, 'loc.tbl')
    open(f, 'w').write(f'>loc\n{loc}\n')
    subprocess.run([CMSEARCH, '--cpu', '1', '--noali', '-E', '1e-5', '--toponly', '--tblout', t, CM, f], check=True, stdout=subprocess.DEVNULL)
    best = None
    for l in open(t):
        if l[0] == '#': continue
        x = l.split(); cmt, s1 = int(x[6]), int(x[8])
        if s1 >= end_hmm - 50 and (best is None or float(x[14]) > best[2]): best = (cmt, s1, float(x[14]))
    if best is None or best[0] < 1520: return 'breaks', h
    cmt, s1, _ = best
    dm, dn = cmt - mt, s1 - end_hmm
    if dm > 0 and dn / dm < 0.8: return 'breaks', h
    if st == '+': return 'ok', [c, lo, a + s1, st, mf, cmt, sc]
    return 'ok', [c, b - s1 + 1, hi, st, mf, cmt, sc]


def body_identity(seqs, h, hmaj):
    def full(x):
        c, lo, hi, st = x[:4]; s = seqs[c][lo - 1:hi]; return s if st == '+' else rc(s)
    a = AL.align(full(hmaj), full(h))[0]
    t, q = a.indices
    m = sum(1 for k in range(t.size) if t[k] >= 0 and q[k] >= 0 and a.target[t[k]] == a.query[q[k]])
    cols = sum(1 for k in range(t.size) if t[k] >= 0 and q[k] >= 0)
    return round(m / cols, 4) if cols else None


def census(seqs, hits):
    regs = [region(seqs, h) for h in hits]
    cnt = Counter(regs)
    top = max(cnt.values())
    imaj = max((i for i, r in enumerate(regs) if cnt[r] == top), key=lambda i: hits[i][6])  # T4-A1: tie -> best bit score
    maj = regs[imaj]
    tc = tail_cols(maj)
    mtail = ''.join(maj[i] for i in tc)
    rec = dict(n_copies=len(hits), maj_tail=mtail, maj_held=len(tc), copies=[])
    if len(tc) < 7:
        rec['cls'] = 'maj_undefined'; return rec
    tcs = set(tc)
    for i, (h, r) in enumerate(zip(hits, regs)):
        if r == maj:
            t, cls, hp = mtail, 'same', False
        else:
            a = AL.align(maj, r)[0]
            # map: target = maj, query = copy
            tt, qq = a.indices
            res = [a.query[qq[k]] for k in range(tt.size) if tt[k] in tcs and qq[k] >= 0]
            t = ''.join(res)
            if len(res) < 7: cls = 'undefined'
            elif t == mtail: cls = 'same'
            elif has_core(mtail) and not has_core(t): cls = 'core-lost'
            else: cls = 'point'
            hp = cls != 'same' and hp_shift(maj, r)
        c = dict(contig=h[0], start=h[1], end=h[2], strand=h[3], mdl=[h[4], h[5]], tail=t, cls=cls, hp=hp)
        if cls != 'same':
            c['body_id'] = body_identity(seqs, h, hits[imaj])
        rec['copies'].append(c)
    order = ['core-lost', 'undefined', 'point', 'same']
    got = {c['cls'] for c in rec['copies']}
    rec['cls'] = next(o for o in order if o in got)
    rec['cls'] = 'uniform' if rec['cls'] == 'same' else rec['cls']
    return rec


def run_fasta(path, tmpdir):
    seqs = read_fasta(path)
    plain = os.path.join(tmpdir, 'g.fna')
    with open(plain, 'w') as o:
        for k, v in seqs.items(): o.write(f'>{k}\n{v}\n')
    hits = copies(plain)
    if len(hits) < 2:
        return dict(n_copies=len(hits), cls='single')
    # T4-A7: re-anchor short hits at the full-CM end; set aside copies that are truncated or clipped by a contig end
    keep, unres = [], []
    for h in hits:
        stt, h2 = cm_end(seqs, h, tmpdir) if h[5] < 1520 else ('ok', h)  # T4-A10
        if stt == 'ok' and len(region(seqs, h2)) < 130: stt = 'edge'
        if stt == 'ok': keep.append(h2)
        else: unres.append(dict(contig=h[0], start=h[1], end=h[2], strand=h[3], mdl=[h[4], h[5]], status=stt))
    if len(keep) < 2:
        return dict(n_copies=len(hits), cls='single', n_resolved=len(keep), unresolved=unres)
    rec = census(seqs, keep)
    rec['n_copies'] = len(hits); rec['n_resolved'] = len(keep); rec['unresolved'] = unres
    return rec


if __name__ == '__main__':
    mode = sys.argv[1]
    tmp = tempfile.mkdtemp(prefix='t4_', dir=os.environ.get('TMPDIR'))
    if mode == 'one':
        print(json.dumps(run_fasta(sys.argv[2], tmp)))
    else:
        idx = {l.split('\t')[0]: l.split('\t')[1] for l in open(ATLAS + 'index/genomes.tsv')}
        out = sys.argv[3]
        done = set()
        if os.path.exists(out):
            done = {json.loads(l)['acc'] for l in open(out)}
        with open(out, 'a') as o:
            for acc in open(sys.argv[2]).read().split():
                if acc in done: continue
                try:
                    r = run_fasta(ATLAS + idx[acc], tmp)
                except Exception as e:
                    r = dict(cls='error', err=repr(e)[:200])
                r['acc'] = acc
                o.write(json.dumps(r) + '\n'); o.flush()
