"""Ortholog-column tally for one query mitogenome.

usage: ortho_tally.py <genome.fna> <tblastn.tsv> <relaln_dir> <outdir>

tblastn.tsv: relatives' proteins (qseqid acc|gene|n) vs the genome, -db_gencode 4,
 outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore sframe'.

Steps: merge HSPs per (gene, strand, phase); resolve overlaps between genes by bitscore;
translate each merged segment with table 4, every codon read through (stops -> X);
add segments as fragments to the relatives' MAFFT alignment of that gene (mafft --addfragments);
for each query codon at an alignment column, record relatives' residues.
Writes <outdir>/<acc>.columns.tsv (one row per query codon in an aligned column).
"""
import sys, os, subprocess, collections
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data import CodonTable

genome_f, tsv, relaln, outdir = sys.argv[1:5]
os.makedirs(outdir, exist_ok=True)
rec = next(SeqIO.parse(genome_f, 'fasta'))
acc = rec.id
G = str(rec.seq).upper()
L = len(G)
t4 = CodonTable.unambiguous_dna_by_id[4].forward_table

# ---- read HSPs
hsps = []
for line in open(tsv):
    q, s, pid, ln, qs, qe, ss, se, ev, bs, sf = line.split('\t')
    gene = q.split('|')[1]
    ss, se, sf, bs = int(ss), int(se), int(sf), float(bs)
    if float(ev) > 1e-5:
        continue
    strand = 1 if sf > 0 else -1
    lo, hi = min(ss, se), max(ss, se)  # 1-based inclusive
    phase = (lo - 1) % 3 if strand == 1 else (L - hi) % 3
    hsps.append(dict(gene=gene, strand=strand, phase=phase, lo=lo, hi=hi, bs=bs))

# ---- merge per gene/strand/phase
segs = []
key = lambda h: (h['gene'], h['strand'], h['phase'])
groups = collections.defaultdict(list)
for h in hsps:
    groups[key(h)].append(h)
for k, lst in groups.items():
    lst.sort(key=lambda h: h['lo'])
    cur = None
    for h in lst:
        if cur and h['lo'] <= cur['hi'] + 1:
            cur['hi'] = max(cur['hi'], h['hi']); cur['bs'] = max(cur['bs'], h['bs']); cur['n'] += 1
        else:
            if cur: segs.append(cur)
            cur = dict(gene=k[0], strand=k[1], phase=k[2], lo=h['lo'], hi=h['hi'], bs=h['bs'], n=1)
    if cur: segs.append(cur)

# ---- resolve overlaps (any frame) between different genes / frames: keep higher bitscore
segs.sort(key=lambda s: -s['bs'])
acc_segs = []
for s in segs:
    ok = True
    for a in acc_segs:
        ov = min(s['hi'], a['hi']) - max(s['lo'], a['lo']) + 1
        if ov > 0.3 * (s['hi'] - s['lo'] + 1):
            ok = False; break
    if ok:
        acc_segs.append(s)

# ---- translate segments (codon-aligned)
def seg_codons(s):
    lo, hi = s['lo'], s['hi']
    if s['strand'] == 1:
        # trim to codon boundaries consistent with phase
        start = lo - 1
        end = hi - ((hi - lo + 1) % 3)
        nt = G[start:end]
        coords = [start + 1 + 3 * i for i in range(len(nt) // 3)]
    else:
        nt_full = str(Seq(G[lo - 1:hi]).reverse_complement())
        nt = nt_full[:len(nt_full) // 3 * 3]
        coords = [hi - 3 * i for i in range(len(nt) // 3)]  # 1-based coord of codon's first base (minus strand)
    cods = [nt[3 * i:3 * i + 3] for i in range(len(nt) // 3)]
    aa = ''.join(t4.get(c, 'X') for c in cods)  # stops (not in forward_table) -> X
    return cods, aa, coords

seg_by_gene = collections.defaultdict(list)
for i, s in enumerate(acc_segs):
    cods, aa, coords = seg_codons(s)
    s.update(cods=cods, aa=aa, coords=coords, id=f'Q{i}')
    seg_by_gene[s['gene']].append(s)

rows = []
for gene, lst in seg_by_gene.items():
    alnf = os.path.join(relaln, f'{gene}.aln')
    fragf = os.path.join(outdir, f'{acc}.{gene}.frag.faa')
    with open(fragf, 'w') as fh:
        for s in lst:
            fh.write(f">{s['id']}\n{s['aa']}\n")
    out = subprocess.run(['mafft', '--quiet', '--thread', '2', '--addfragments', fragf, alnf],
                         capture_output=True, text=True, check=True).stdout
    aln = {}
    name = None
    for line in out.splitlines():
        if line.startswith('>'):
            name = line[1:].split()[0]; aln[name] = []
        else:
            aln[name].append(line.strip())
    aln = {k: ''.join(v).upper() for k, v in aln.items()}
    rel = {k: v for k, v in aln.items() if not k.startswith('Q')}
    nrel = len(rel)
    rel_rows = list(rel.values())
    first = [next((i for i, ch in enumerate(r) if ch != "-"), -1) for r in rel_rows]
    for s in lst:
        a = aln[s['id']]
        ri = 0
        for col, ch in enumerate(a):
            if ch == '-':
                continue
            res = [r[col] for r in rel_rows if r[col] != '-']
            c = collections.Counter(res)
            top, ntop = (c.most_common(1)[0] if c else ('-', 0))
            nstart = sum(1 for fi in first if fi == col)  # relatives whose first residue is this column
            rows.append((acc, gene, s['strand'], s['coords'][ri], s['cods'][ri], s['aa'][ri],
                         nrel, len(res), top, ntop, ''.join(sorted(res)), col, nstart))
            ri += 1
        assert ri == len(s['aa']), (gene, s['id'], ri, len(s['aa']))

with open(os.path.join(outdir, f'{acc}.columns.tsv'), 'w') as fh:
    fh.write('acc\tgene\tstrand\tcoord\tcodon\tt4aa\tnrel\tnres\ttop\tntop\tresidues\tcol\tnstart\n')
    for r in rows:
        fh.write('\t'.join(map(str, r)) + '\n')
with open(os.path.join(outdir, f'{acc}.segments.tsv'), 'w') as fh:
    for s in sorted(acc_segs, key=lambda s: (s['gene'], s['lo'])):
        fh.write(f"{s['gene']}\t{s['strand']}\t{s['lo']}\t{s['hi']}\t{s['bs']}\t{s['n']}\t{len(s['aa'])}\n")
print(acc, 'segments', len(acc_segs), 'genes', len(seg_by_gene), 'codon rows', len(rows))
