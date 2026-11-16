# Task 3 eval: does any single target HSP span the fusion junction (>=20 aa on both sides)? gaps near junction?
# diamond fields: qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen qseq_gapped sseq_gapped
import sys, csv, collections
qtsv, hits_p, label = sys.argv[1:4]
Q = {r['qid']: r for r in csv.DictReader(open(qtsv), delimiter='\t')}
best = {}  # qid -> dict
info = collections.defaultdict(lambda: {'left': False, 'span': False, 'span_clean': False, 'best_span': None, 'right_only': False})
for l in open(hits_p):
    f = l.rstrip('\n').split('\t')
    qid, sid = f[0], f[1]
    if sid.split('|')[0] == Q[qid]['id']: continue   # own phage
    qs, qe, ss = int(f[4]), int(f[5]), int(f[6]); ev = float(f[8]); bits = float(f[9])
    qg, sg = f[12], f[13]
    J = int(Q[qid]['junction'])
    d = info[qid]
    if qs <= J - 19: d['left'] = True
    if qs >= J - 5 and qe >= J + 20: d['right_only'] = True
    if not (qs <= J - 19 and qe >= J + 20): continue
    # walk alignment
    qp = qs - 1; sp = ss - 1; gaps_win = 0; idL = nL = idR = nR = 0; maxrun = run = 0
    for a, b in zip(qg, sg):
        if a != '-': qp += 1
        if b != '-': sp += 1
        inwin = J - 9 <= qp <= J + 10
        if a == '-' or b == '-':
            if inwin:
                gaps_win += 1; run += 1; maxrun = max(maxrun, run)
        else:
            run = 0
            if qp <= J: nL += 1; idL += a == b
            else: nR += 1; idR += a == b
    clean = maxrun <= 2
    d['span'] = True
    if clean: d['span_clean'] = True
    rec = (bits, sid, ev, f[2], qs, qe, J, gaps_win, maxrun, round(100*idL/max(nL,1)), round(100*idR/max(nR,1)))
    if d['best_span'] is None or bits > d['best_span'][0]: d['best_span'] = rec
cls = collections.defaultdict(lambda: collections.Counter())
for qid, r in Q.items():
    c = r['class']; d = info[qid]
    cls[c]['n'] += 1
    cls[c]['left_hit'] += d['left']
    cls[c]['span'] += d['span']
    cls[c]['span_clean'] += d['span_clean']
    cls[c]['right_hit'] += d['right_only'] or d['span']
print(f'== {label}')
for c in ('T1skip', 'T1ctrl', 'Rskip'):
    x = cls[c]
    print(f"{c}: n={x['n']} queries with HSP >=20 aa left of junction={x['left_hit']}; spanning HSP={x['span']} (of which gap run<=2 in +-10 aa window: {x['span_clean']}); any HSP reaching >=20 aa right of junction={x['right_hit']}; span/left-hit={x['span']}/{x['left_hit']}")
print('\nper-query best spanning HSP (bits, target, evalue, pident, qs, qe, J, gapcols_in_win, max_gaprun, idL%, idR%):')
for qid in sorted(Q):
    d = info[qid]
    if d['best_span']: print(' ', qid, Q[qid]['codon'], Q[qid]['next3'], d['best_span'])
