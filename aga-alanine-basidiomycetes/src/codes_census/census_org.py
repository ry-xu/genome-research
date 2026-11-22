"""Codon-meaning census of organelle genomes: which residue do strictly conserved profile columns want where each
codon sits? (mito-codes' conserved-site test, on this project's gene tables and lineage profiles.)

usage: python census_org.py prefix profiles.hmm out.npz [workers] [exclude_taxon]
reads prefix.records.tsv / prefix.genes.tsv.gz (org_extract.py / orf_genes.py). Genes are translated as in
scan_org.py (declared table with sense-to-stop erased; in-frame canonical stops -> X), each gene's best domain over
all profiles is kept, and every match-state residue with posterior >= 0.9 except the start codon gives one site (codon, model, column).
Sites are de-duplicated per genome. Genomes whose taxonomy contains exclude_taxon (e.g. Metazoa) are skipped.
out.npz: acc, table, taxonomy [g]; site_genome, site_codon, site_model, site_col [n]; models [m];
emax, earg [m, max M+1] (the column's max match emission and its residue index in AA).
"""
import sys, os
from multiprocessing import Pool
import numpy as np
import pyhmmer
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))  # runs under python -I
from codes import CI, translate_ext
import scan_org

AA = 'ACDEFGHIKLMNPQRSTVWY'
MIDX = None


def load(path):
    global MIDX
    scan_org.load(path)
    MIDX = {(h.name.decode() if isinstance(h.name, bytes) else h.name): i for i, h in enumerate(scan_org.HMMS)}


def census(item):
    acc, genes = item
    seqs, meta = [], []
    for g in genes:
        idx, gene, prod, table, strand, cds_len, ext_len, nt = g
        if set(nt) - set('ACGT'):
            nt = ''.join(c if c in 'ACGT' else 'N' for c in nt)
        prot, cut = translate_ext(nt, cds_len, table)
        if len(prot) < 30:
            continue
        meta.append(nt)
        seqs.append(pyhmmer.easel.TextSequence(name=str(len(seqs)).encode(), sequence=prot))
    sites = set()
    if not seqs:
        return acc, sites
    block = pyhmmer.easel.TextSequenceBlock(seqs).digitize(scan_org.ALPHA)
    best = {}
    for hits in pyhmmer.hmmsearch(scan_org.HMMS, block, cpus=1, E=1e-5):
        model = hits.query.name
        model = model.decode() if isinstance(model, bytes) else model
        for hit in hits.reported:
            i = int(hit.name)
            for d in hit.domains.reported:
                if d.score > best.get(i, (-1e9,))[0]:
                    best[i] = (d.score, d.alignment, model)
    for i, (sc, al, model) in best.items():
        nt, m = meta[i], MIDX[model]
        t, k = al.target_from - 1, al.hmm_from
        for hc, tc, p in zip(al.hmm_sequence, al.target_sequence, al.posterior_probabilities):
            if tc == '-':  # deletion: model column advances
                k += 1
                continue
            if hc == '.':  # insertion
                t += 1
                continue
            if p in '*9' and t > 0:  # not the start codon (ATA/ATT/GTG starts would vote Met)
                c = nt[3 * t:3 * t + 3]
                if c in CI:
                    sites.add((CI[c], m, k))
            t += 1
            k += 1
    return acc, sites


def main(prefix, hmm, out, workers=4, exclude=None):
    rec = {}
    for line in open(prefix + '.records.tsv'):
        f = line.rstrip('\n').split('\t')
        rec[f[0]] = (f[2], int(f[3]))
    keep = {a for a, (tx, t) in rec.items() if not (exclude and exclude in tx.split(';'))}
    load(hmm)
    models = [h.name.decode() if isinstance(h.name, bytes) else h.name for h in scan_org.HMMS]
    Mx = max(h.M for h in scan_org.HMMS) + 1
    emax = np.zeros((len(models), Mx), np.float32)
    earg = np.zeros((len(models), Mx), np.int8)
    for i, h in enumerate(scan_org.HMMS):
        em = np.asarray(h.match_emissions)[:, :20]  # row k = match state k (row 0 unused)
        emax[i, :h.M + 1] = em.max(1)
        earg[i, :h.M + 1] = em.argmax(1)
    accs, sg, sc, sm, sk = [], [], [], [], []
    with Pool(int(workers), initializer=load, initargs=(hmm,)) as pool:
        for n, (acc, sites) in enumerate(pool.imap_unordered(census, scan_org.genomes(prefix, keep), chunksize=4)):
            gi = len(accs)
            accs.append(acc)
            for c, m, k in sites:
                sg.append(gi); sc.append(c); sm.append(m); sk.append(k)
            if n % 1000 == 0:
                print(n, len(sg), flush=True)
    np.savez_compressed(out, acc=np.array(accs), table=np.array([rec[a][1] for a in accs]),
                        taxonomy=np.array([rec[a][0] for a in accs]), site_genome=np.array(sg, np.int32),
                        site_codon=np.array(sc, np.int8), site_model=np.array(sm, np.int16),
                        site_col=np.array(sk, np.int16), models=np.array(models), emax=emax, earg=earg)
    print('done', len(accs), 'genomes', len(sg), 'sites')


if __name__ == '__main__':
    main(*sys.argv[1:])
