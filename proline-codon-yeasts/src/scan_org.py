"""Signals A and B for organelle genomes (blinded translation, full-length gene profiles).

usage: python scan_org.py prefix profiles.hmm out_prefix [workers] [acc_list]
reads prefix.records.tsv / prefix.genes.tsv.gz from org_extract.py.
writes out_prefix.npz: acc, N_in [g,64] (in-frame codons at aligned positions, posterior >= 0.9, not the last 3
       aligned positions of a domain), N_off [g,64] (the two off-frame triplets at the same positions)
       out_prefix.ends.tsv: one row per gene with a best alignment (any end point):
       acc, gene_idx, gene, model, hmm_to, M, last (codon index of last aligned residue), cut (codon index of the first
       canonical stop at/after the annotated end; -1 if none in the extension), cds_codons, pre (last 5 aligned
       codons), post (up to 40 codons after the last aligned one), body (aligned in-frame span)
"""
import sys, os, gzip
from collections import defaultdict
from multiprocessing import Pool
import numpy as np
import pyhmmer
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # runs under python -I
from codes import CI, translate_ext

HMMS = None
MLEN = None
ALPHA = pyhmmer.easel.Alphabet.amino()


def load(path):
    global HMMS, MLEN
    with pyhmmer.plan7.HMMFile(path) as f:
        HMMS = list(f)
    MLEN = {h.name.decode() if isinstance(h.name, bytes) else h.name: h.M for h in HMMS}


def scan(item):
    acc, genes = item
    seqs, meta = [], []
    for g in genes:
        idx, gene, prod, table, strand, cds_len, ext_len, nt = g
        if set(nt) - set('ACGT'):
            nt = ''.join(c if c in 'ACGT' else 'N' for c in nt)
        prot, cut = translate_ext(nt, cds_len, table)
        if len(prot) < 30:
            continue
        meta.append((idx, gene, nt, cut, cds_len // 3))
        seqs.append(pyhmmer.easel.TextSequence(name=str(len(seqs)).encode(), sequence=prot))
    N_in = np.zeros(64, np.int32)
    N_off = np.zeros(64, np.int32)
    ends = []
    if not seqs:
        return acc, N_in, N_off, ends
    block = pyhmmer.easel.TextSequenceBlock(seqs).digitize(ALPHA)
    best = {}  # seq index -> best domain over all models
    for hits in pyhmmer.hmmsearch(HMMS, block, cpus=1, E=1e-5):
        model = hits.query.name if hasattr(hits, 'query') else None
        model = model.decode() if isinstance(model, bytes) else model
        for hit in hits.reported:
            i = int(hit.name)
            for d in hit.domains.reported:
                if d.score > best.get(i, (-1e9,))[0]:
                    best[i] = (d.score, d.alignment, model)
    for i, (sc, al, model) in best.items():
        idx, gene, nt, cut, cdsc = meta[i]
        pp = al.posterior_probabilities
        t = al.target_from - 1
        pos = []  # target codon indices at match states with posterior >= 0.9
        for hc, tc, p in zip(al.hmm_sequence, al.target_sequence, pp):
            if tc == '-':
                continue
            if hc != '.' and (p == '*' or p == '9'):
                pos.append(t)
            t += 1
        last = al.target_to - 1
        for t in pos:
            if t > last - 3:
                continue
            c = nt[3 * t:3 * t + 3]
            if c in CI:
                N_in[CI[c]] += 1
            for o in (1, 2):
                c = nt[3 * t + o:3 * t + o + 3]
                if c in CI:
                    N_off[CI[c]] += 1
        M = MLEN[model]
        if True:  # every best alignment; scoring decides how close to the model end is close enough
            pre = nt[3 * max(last - 4, 0):3 * (last + 1)]
            post = nt[3 * (last + 1):3 * (last + 41)]
            body = nt[3 * (al.target_from - 1):3 * (last + 1)]
            ends.append((acc, idx, gene, model, al.hmm_to, M, last, -1 if cut is None else cut, cdsc, pre, post, body))
    return acc, N_in, N_off, ends


def genomes(prefix, keep=None):
    cur, buf = None, []
    with gzip.open(prefix + '.genes.tsv.gz', 'rt') as f:
        for line in f:
            acc, idx, gene, prod, table, strand, cds_len, ext_len, nt = line.rstrip('\n').split('\t')
            if acc != cur:
                if buf and (keep is None or cur in keep):
                    yield cur, buf
                cur, buf = acc, []
            buf.append((int(idx), gene, prod, int(table), strand, int(cds_len), int(ext_len), nt))
    if buf and (keep is None or cur in keep):
        yield cur, buf


def main(prefix, hmm, out, workers=6, acc_list=None):
    keep = set(open(acc_list).read().split()) if acc_list else None
    accs, A, B = [], [], []
    with open(out + '.ends.tsv', 'w') as fe, Pool(int(workers), initializer=load, initargs=(hmm,)) as pool:
        fe.write('acc\tgene_idx\tgene\tmodel\thmm_to\tM\tlast\tcut\tcds_codons\tpre\tpost\tbody\n')
        for n, (acc, a, b, ends) in enumerate(pool.imap_unordered(scan, genomes(prefix, keep), chunksize=4)):
            accs.append(acc)
            A.append(a)
            B.append(b)
            for e in ends:
                fe.write('\t'.join(map(str, e)) + '\n')
            if n % 1000 == 0:
                print(n, flush=True)
    np.savez_compressed(out + '.npz', acc=np.array(accs), N_in=np.array(A), N_off=np.array(B))
    print('done', len(accs))


if __name__ == '__main__':
    main(*sys.argv[1:])
