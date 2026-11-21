"""Raw-read check of codons at conserved sites: do the reads show the assembled codon?

usage: python reads_codons.py prefix acc reads.fq CODON[,CODON...] out.tsv census.npz [census2.npz ...]
Aligns acc's genes to the profiles as census_org.py does, keeping each target codon's position and (model, column);
maps reads (minimap2 -ax sr) to the genes' nucleotide sequences and counts, per target codon, the reads covering all
three bases and what triplet they show (a deletion counts as a non-matching triplet; an insertion inside the codon
is not seen). Column consensus/strictness as in sitescan_org.py, computed from the census files given: pass the same
files as the census (mito, gbmito, gbmitoorf), or the site set differs.
Env: MIN_MAPQ (20), MAX_DIV (unset; e.g. 0.05 drops reads with NM > 5% of length), MM_PRESET (sr).
Writes one row per target codon: gene, codon index, codon, model:column, consensus residue, share, strict,
reads, reads showing the codon, other triplets."""
import sys, os, subprocess, tempfile
from collections import Counter, defaultdict
import pyhmmer
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import census_org, scan_org
from sitescan_org import load as load_census, AA
from codes import translate_ext


def positions(prefix, acc, targets, hmm):
    census_org.load(hmm)
    items = list(scan_org.genomes(prefix, {acc}))
    assert items, acc
    _, genes = items[0]
    seqs, nts = [], []
    for g in genes:
        nt = ''.join(c if c in 'ACGT' else 'N' for c in g[7])
        prot, _ = translate_ext(nt, g[5], g[3])
        if len(prot) < 30:
            continue
        nts.append((g[1] or f'g{g[0]}', nt))
        seqs.append(pyhmmer.easel.TextSequence(name=str(len(seqs)).encode(), sequence=prot))
    block = pyhmmer.easel.TextSequenceBlock(seqs).digitize(scan_org.ALPHA)
    best = {}
    for hits in pyhmmer.hmmsearch(scan_org.HMMS, block, cpus=1, E=1e-5):
        model = hits.query.name
        model = model.decode() if isinstance(model, bytes) else model
        for hit in hits.reported:
            for d in hit.domains.reported:
                i = int(hit.name)
                if d.score > best.get(i, (-1e9,))[0]:
                    best[i] = (d.score, d.alignment, model)
    out = []
    for i, (sc, al, model) in best.items():
        nt = nts[i][1]
        t, k = al.target_from - 1, al.hmm_from
        for hc, tc, p in zip(al.hmm_sequence, al.target_sequence, al.posterior_probabilities):
            if tc == '-':
                k += 1
                continue
            if hc == '.':
                t += 1
                continue
            c = nt[3 * t:3 * t + 3]
            if p in '*9' and t > 0 and c in targets:
                out.append((i, t, c, census_org.MIDX[model], k))
            t += 1
            k += 1
    return nts, out


def main(prefix, acc, reads, codons, out, *census):
    D = load_census(census)
    hmm = 'data/mito_profiles.hmm' if 'mito' in os.path.basename(census[0]) else 'data/plastid_profiles.hmm'
    nts, sites = positions(prefix, acc, set(codons.split(',')), hmm)
    tmp = tempfile.mkdtemp(dir=os.environ.get('TMPDIR'))
    ref = os.path.join(tmp, 'ref.fa')
    with open(ref, 'w') as fh:
        for j, (name, nt) in enumerate(nts):
            fh.write(f'>g{j}\n{nt}\n')
    want = defaultdict(dict)  # gene -> {start: key}
    for i, t, c, m, k in sites:
        want[i][3 * t] = (i, t)
    seen = defaultdict(Counter)
    p = subprocess.Popen(['minimap2', '-ax', os.environ.get('MM_PRESET', 'sr'), '-t', '2', '--secondary=no', ref, reads], stdout=subprocess.PIPE,
                         stderr=subprocess.DEVNULL, text=True)
    nmap = 0
    for line in p.stdout:
        if line[0] == '@':
            continue
        f = line.split('\t', 11)
        flag = int(f[1])
        if flag & 0x904 or int(f[4]) < int(os.environ.get('MIN_MAPQ', '20')):
            continue
        if 'MAX_DIV' in os.environ:  # optional: drop reads with NM above this fraction of read length (cross-mapped reads)
            nm = next((int(t[5:]) for t in (f[11].split('\t') if len(f) > 11 else []) if t.startswith('NM:i:')), 0)
            if nm > float(os.environ['MAX_DIV']) * len(f[9]):
                continue
        nmap += 1
        g = int(f[2][1:])
        if g not in want:
            continue
        pos, seq, cig = int(f[3]) - 1, f[9], f[5]
        rmap = {}
        q = 0
        num = ''
        for ch in cig:
            if ch.isdigit():
                num += ch
                continue
            n = int(num); num = ''
            if ch in 'M=X':
                for x in range(n):
                    rmap[pos + x] = seq[q + x]
                pos += n; q += n
            elif ch in 'IS':
                q += n
            elif ch in 'DN':
                for x in range(n):
                    rmap[pos + x] = '-'
                pos += n
        for s0, key in want[g].items():
            if s0 in rmap and s0 + 2 in rmap and s0 + 1 in rmap:
                tri = rmap[s0] + rmap[s0 + 1] + rmap[s0 + 2]
                seen[key][tri] += 1
    p.wait()
    K = D['K']
    with open(out, 'w') as fh:
        fh.write('gene\tcodon_index\tcodon\tmodel_col\tcons\tshare\tstrict\treads\tsame\tother\n')
        for i, t, c, m, k in sorted(sites):
            col = m * K + k
            sh = D['freq'][col].max() / D['tot'][col] if D['tot'][col] > 0 else 0
            cnt = seen[(i, t)]
            n = sum(cnt.values())
            fh.write(f"{nts[i][0]}\t{t}\t{c}\t{D['models'][m]}:{k}\t{AA[D['cons'][col]]}\t{sh:.2f}\t{int(D['strict'][col])}\t"
                     f"{n}\t{cnt[c]}\t{dict((x, y) for x, y in cnt.most_common(3) if x != c)}\n")
    print(f'{acc}: {nmap} reads mapped to genes; {len(sites)} target codons at aligned sites', file=sys.stderr)


if __name__ == '__main__':
    main(*sys.argv[1:])
