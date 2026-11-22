"""Tally Sage PSMs on the three versions of each mitochondrial ORF (build_db.py: |R table 4, |A AGA = Ala,
|C CGN = Ala control). usage: python tally_aga.py results.sage.tsv db.fasta [q=0.01]
Sage lists every protein whose digest yields a peptide, so the set of versions a peptide maps to says what it needs:
  RAC  no AGA or CGN codon inside, tryptic in every version: neutral (shows the protein is seen at all)
  RC   needs AGA = Arg (contains it, or ends or starts at it)        -> Arg reading of AGA
  A    needs AGA = Ala                                                 -> Ala reading of AGA
  RA   needs CGN = Arg (positive control: the standard reading)
  C    needs CGN = Ala (negative control: should be decoy-level)
  R    needs both AGA = Arg and CGN = Arg
Counts target and decoy PSMs at spectrum q <= q, and lists every non-neutral target peptide."""
import sys, csv, re
from collections import defaultdict


def main(res, db, q=0.01):
    q = float(q)
    seqs = {}
    name = None
    for line in open(db):
        if line.startswith('>'):
            name = line[1:].strip()
        elif name.startswith('mito|'):
            seqs[name] = seqs.get(name, '') + line.strip()
    psm = defaultdict(lambda: [0, 0])
    peps = defaultdict(lambda: [0, set(), None, set()])
    for r in csv.DictReader(open(res), delimiter='\t'):
        if float(r['spectrum_q']) > q:
            continue
        prots = r['proteins'].split(';')
        dec = r['label'] == '-1'
        mito = [p for p in prots if p.split('|')[0].replace('rev_', '') == 'mito']
        if not mito:
            continue
        nuc = len(mito) < len(prots)
        vers = ''.join(sorted({p.rsplit('|', 1)[1] for p in mito}, key='RAC'.index))
        key = vers + ('+nuc' if nuc else '')
        psm[key][dec] += 1
        if not dec and vers != 'RAC':
            pep = re.sub(r'\[[^\]]*\]-?', '', r['peptide'])
            e = peps[(vers, pep)]
            e[0] += 1; e[1].add(r['filename']); e[3] |= {p.split('|')[1] for p in mito}
    print(f'PSMs at spectrum q <= {q} on mitochondrial ORFs, by versions matched (target / decoy):')
    for k in sorted(psm, key=lambda k: (-sum(psm[k]), k)):
        print(f'  {k:8s} {psm[k][0]:6d} {psm[k][1]:6d}')
    print('non-neutral target peptides: versions, PSMs, files, peptide, ORFs, residue at AGA/CGN positions')
    for (vers, pep), (n, files, _, orfs) in sorted(peps.items(), key=lambda x: (x[0][0], -x[1][0])):
        where = []
        for o in sorted(orfs):
            for v in vers[:1]:
                s = seqs.get(f'mito|{o}|{v}', '')
                i = s.find(pep)
                if i >= 0:
                    where.append(f'{o}:{i}')
        print(f'  {vers:3s} {n:4d} {len(files):3d}  {pep}  {",".join(where[:4])}')


if __name__ == '__main__':
    main(*sys.argv[1:])
