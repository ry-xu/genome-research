"""Parse relatives' fasta_cds_na, assign core gene names, translate with NCBI table 4 (own translation),
compare with GenBank translation (fasta_cds_aa), write one protein FASTA per gene.
usage: prep_relatives.py <dl_dir> <outdir>"""
import sys, os, re, glob
from Bio import SeqIO
from Bio.Seq import Seq
dl, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)
pats = [
 ('cox1', r'cytochrome (c )?oxidase subunit (1|I)$'), ('cox2', r'cytochrome (c )?oxidase subunit (2|II)$'),
 ('cox3', r'cytochrome (c )?oxidase subunit (3|III)$'), ('cob', r'(apo)?cytochrome b$'),
 ('atp6', r'ATP synthase.*(subunit (6|a)|subunit a)$'), ('atp8', r'ATP synthase.*subunit 8$'),
 ('atp9', r'ATP synthase.*subunit (9|c)$'),
 ('nad4L', r'NADH.*subunit 4L$'), ('nad1', r'NADH.*subunit 1$'), ('nad2', r'NADH.*subunit 2$'),
 ('nad3', r'NADH.*subunit 3$'), ('nad4', r'NADH.*subunit 4$'), ('nad5', r'NADH.*subunit 5$'),
 ('nad6', r'NADH.*subunit 6$'), ('rps3', r'ribosomal protein S3$')]
genes = {}
log = open(os.path.join(out, 'prep.log'), 'w')
for f in sorted(glob.glob(os.path.join(dl, '*.fasta_cds_na'))):
    acc = os.path.basename(f).replace('.fasta_cds_na', '')
    aa = {}
    for r in SeqIO.parse(f.replace('_na', '_aa'), 'fasta'):
        pid = re.search(r'\[protein_id=([^\]]+)\]', r.description)
        if pid: aa[pid.group(1)] = str(r.seq)
    for r in SeqIO.parse(f, 'fasta'):
        m = re.search(r'\[protein=([^\]]+)\]', r.description)
        prot = m.group(1).strip() if m else ''
        g = None
        for name, p in pats:
            if re.search(p, prot, re.I): g = name; break
        if g is None: continue
        if 'partial' in r.description: log.write(f'{acc} {g} partial, skipped\n'); continue
        s = str(r.seq).upper()
        if len(s) % 3: log.write(f'{acc} {g} length {len(s)} not multiple of 3\n')
        tr = str(Seq(s[:len(s)//3*3]).translate(table=4))
        if tr.endswith('*'): tr = tr[:-1]
        # initiator codon read as M
        tr = 'M' + tr[1:]
        pid = re.search(r'\[protein_id=([^\]]+)\]', r.description)
        gb = aa.get(pid.group(1)) if pid else None
        nstop = tr.count('*')
        diff = None if gb is None else sum(a != b for a, b in zip(tr, gb)) + abs(len(tr) - len(gb))
        log.write(f'{acc} {g} len={len(tr)} internal_stops={nstop} diff_vs_genbank={diff}\n')
        if nstop: continue
        genes.setdefault(g, []).append((acc, tr))
for g, lst in genes.items():
    with open(os.path.join(out, f'{g}.faa'), 'w') as fh:
        seen = {}
        for acc, tr in lst:
            seen[acc] = seen.get(acc, 0) + 1
            fh.write(f'>{acc}|{g}|{seen[acc]}\n{tr}\n')
    print(g, len(lst), sorted(set(a for a, _ in lst)).__len__())
