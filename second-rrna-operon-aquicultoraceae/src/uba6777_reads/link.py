"""Physical-link test for NODE_6 (rules in ux/NOTES.md, "Physical-link test" and its "amendment before rerun").

usage: python3 -I link.py SPADES_DIR ANCHOR.fa CEN.fa FAM.fa OUTDIR
Needs blastn/makeblastdb on PATH and cmsearch (Infernal) on PATH.
1. rRNA models (RF00177 16S, RF02541 23S) on scaffolds >= 1500 nt.
2. Positive control: scaffolds with a full-length 16S (model 1-10 .. >= 1400) and >= 500 nt of 23S on the same strand, with the
   23S start within 1.5 kb after the 16S end. Need >= 5.
3. Fingerprint of every 23S (RF02541 span >= 1500): dc-megablast, identity summed over HSPs per subject, subject alignment
   >= 1000 nt, against the census (CEN.fa) and the 3 family normal loci (FAM.fa). NODE_6-like = family is the best hit at
   0.826-0.846, and the best census hit is Myxococcota_A at 0.798-0.818.
4. Anchor (host flank + 16S 5' end) in scaffolds at >= 99% over >= 300 nt. The 23S downstream on the same scaffold is
   fingerprinted.
5. If no anchor scaffold carries a 23S downstream: walk the GFA (assembly_graph_with_scaffolds.gfa) from the segment
   holding the anchor's 3' end, in its orientation, following links up to 6 kb. Each reached segment >= 100 nt is assigned
   to the 23S scaffolds that contain it (>= 99% over >= 95% of the segment).
Writes OUTDIR/anchor_scaffolds.fa, OUTDIR/lsu_fingerprints.tsv and OUTDIR/node6_like.fa, and prints the reading.
"""
import collections
import gzip
import os
import subprocess
import sys
import tempfile

CM = "cmsearch"
RF = "data/external/rfam/"
SSU_ALL = "data/external/gtdb/ssu_all_r232.fna.gz"
RCT = str.maketrans("ACGTNacgtn", "TGCANtgcan")
FAM_BAND, CEN_BAND = (0.826, 0.846), (0.798, 0.818)


def rc(s):
    return s.translate(RCT)[::-1]


def fasta(p):
    d, n = {}, None
    for l in open(p):
        if l[0] == ">":
            n = l[1:].split()[0]
            d[n] = []
        else:
            d[n].append(l.strip())
    return {k: "".join(v) for k, v in d.items()}


def write(p, recs):
    with open(p, "w") as o:
        for k, v in recs:
            o.write(f">{k}\n{v}\n")


def cms(model, fa, tmp):
    tb = os.path.join(tmp, os.path.basename(model) + ".tbl")
    subprocess.run([CM, "--cpu", "2", "--tblout", tb, model, fa], stdout=subprocess.DEVNULL, check=True)
    hits = collections.defaultdict(list)
    for l in open(tb):
        if l[0] == "#":
            continue
        f = l.split()
        if f[16] == "!":
            hits[f[0]].append((int(f[5]), int(f[6]), int(f[7]), int(f[8]), f[9]))  # mdl from, mdl to, seq from, seq to, strand
    return hits


def blast(q, s, task=None, db=False, extra=()):
    cmd = ["blastn", "-query", q, "-db" if db else "-subject", s,
           "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send nident"]
    if task:
        cmd += ["-task", task]
    out = subprocess.run(cmd + list(extra), capture_output=True, text=True, check=True).stdout
    return [l.split("\t") for l in out.splitlines()]


def summed(rows):
    """Per (query, subject): identity summed over HSPs; keep subjects with summed subject alignment >= 1000 nt."""
    acc = collections.defaultdict(lambda: [0, 0, 0])
    for r in rows:
        a = acc[(r[0], r[1])]
        a[0] += int(r[8])
        a[1] += int(r[3])
        a[2] += abs(int(r[7]) - int(r[6])) + 1
    best = {}
    for (q, s), (nid, ln, sal) in acc.items():
        if sal >= 1000 and (q not in best or nid / ln > best[q][1]):
            best[q] = (s, nid / ln, sal)
    return best


def taxa(accs):
    want = {a[3:] if a[:3] in ("GB_", "RS_") else a for a in accs}
    out = {}
    with gzip.open(SSU_ALL, "rt") as f:
        for l in f:
            if l[0] != ">":
                continue
            a = l[1:].split()[0].split("~")[0]
            k = a[3:] if a[:3] in ("GB_", "RS_") else a
            if k in want and k not in out:
                out[k] = l.split(None, 1)[1].split(" [")[0]
    return out


def main(sp, anchor, cen, fam, outdir):
    os.makedirs(outdir, exist_ok=True)
    tmp = tempfile.mkdtemp(dir=os.path.dirname(os.path.abspath(anchor)))
    scf_p = os.path.join(sp, "scaffolds.fasta")
    scf = fasta(scf_p)
    big = os.path.join(tmp, "big.fa")
    write(big, ((k, v) for k, v in scf.items() if len(v) >= 1500))
    ssu = cms(RF + "RF00177.cm", big, tmp)
    lsu = cms(RF + "RF02541.cm", big, tmp)

    linked = 0
    for k in ssu:
        for a in ssu[k]:
            if a[0] > 10 or a[1] < 1400:
                continue
            for b in lsu.get(k, []):
                if b[4] != a[4] or b[1] - b[0] < 500:
                    continue
                gap = (b[2] - a[3]) if a[4] == "+" else (a[3] - b[2])
                if 0 < gap <= 1500:
                    linked += 1
                    break
    print(f"positive control: {linked} scaffolds with full 16S + >=500 nt 23S within 1.5 kb (need >= 5)")

    # 23S regions (RF02541 span >= 1500), oriented 5'->3'
    regions = {}
    for k, hs in lsu.items():
        for i, (m0, m1, s0, s1, st) in enumerate(hs):
            if m1 - m0 + 1 < 1500:
                continue
            seq = scf[k][s0 - 1:s1] if st == "+" else rc(scf[k][s1 - 1:s0])
            regions[f"{k}#{i}"] = (k, m0, m1, s0, s1, st, seq)
    rq = os.path.join(tmp, "lsu.fa")
    write(rq, ((n, r[6]) for n, r in regions.items()))
    db = os.path.join(tmp, "cendb")
    subprocess.run(["makeblastdb", "-in", cen, "-dbtype", "nucl", "-out", db], stdout=subprocess.DEVNULL, check=True)
    bc = summed(blast(rq, db, "dc-megablast", db=True, extra=("-max_target_seqs", "50", "-num_threads", "2")))
    bf = summed(blast(rq, fam, "dc-megablast"))
    tx = taxa({v[0].split("|", 1)[1] for v in bc.values()})
    fp = {}
    with open(os.path.join(outdir, "lsu_fingerprints.tsv"), "w") as o:
        o.write("region\tlen\tmodel\tfam_id\tcen_best\tcen_id\tcen_taxon\tnode6_like\n")
        for n, r in regions.items():
            f_id = bf.get(n, ("", 0.0, 0))[1]
            c = bc.get(n, ("", 0.0, 0))
            cacc = c[0].split("|", 1)[1] if c[0] else ""
            t = tx.get(cacc[3:] if cacc[:3] in ("GB_", "RS_") else cacc, "")
            like = (f_id > c[1] and FAM_BAND[0] <= f_id <= FAM_BAND[1] and CEN_BAND[0] <= c[1] <= CEN_BAND[1]
                    and "p__Myxococcota_A" in t)
            fp[n] = (f_id, cacc, c[1], t, like)
            o.write(f"{n}\t{len(r[6])}\t{r[1]}-{r[2]}\t{f_id:.3f}\t{cacc}\t{c[1]:.3f}\t{t}\t{like}\n")
    like = [n for n in fp if fp[n][4]]
    write(os.path.join(outdir, "node6_like.fa"), ((n, regions[n][6]) for n in like))
    print(f"23S regions >= 1500 model nt: {len(regions)}; NODE_6-like: {like}")
    near = sorted(fp.items(), key=lambda x: -x[1][0])[:5]
    print("closest to the family loci:", [(n, f"{v[0]:.3f}", f"{v[2]:.3f}", v[3].split(";")[1:3]) for n, v in near])

    hits = [h for h in blast(anchor, scf_p) if float(h[2]) >= 99 and int(h[3]) >= 300]
    print("anchor hits (>=99%, >=300 nt):", [(h[1], h[2], h[3], h[4], h[5], h[6], h[7]) for h in hits])
    write(os.path.join(outdir, "anchor_scaffolds.fa"), ((h[1], scf[h[1]]) for h in hits))
    reached = set()
    for h in hits:
        k, st, en = h[1], int(h[6]), int(h[7])
        plus = st < en
        down = [n for n, r in regions.items() if r[0] == k and (r[5] == "+") == plus
                and ((r[3] > en) if plus else (r[3] < en))]
        print(f"  {k} len {len(scf[k])}: anchor {st}-{en}; downstream 23S regions: "
              f"{[(n, f'{fp[n][0]:.3f}', fp[n][1], f'{fp[n][2]:.3f}', fp[n][4]) for n in down]}")
        reached.update(down)

    if linked < 5:
        print("READING: UNINFORMATIVE (positive control failed)")
        return
    if not hits:
        print("READING: UNINFORMATIVE (anchor not assembled)")
        return
    if not reached:
        reached = walk(sp, anchor, regions, tmp)
        print("GFA walk reached 23S regions:", [(n, f"{fp[n][0]:.3f}", fp[n][1], f"{fp[n][2]:.3f}", fp[n][4]) for n in sorted(reached)])
        if not reached:
            print("READING: UNINFORMATIVE (no path from the anchor into a 23S)")
            return
    # distinct 23S = regions < 99% to each other; regions from one scaffold family collapse
    distinct = cluster(sorted(reached), regions, tmp)
    good = [c for c in distinct if any(fp[n][4] for n in c)]
    if len(distinct) > 1:
        print(f"READING: AMBIGUOUS ({len(distinct)} distinct 23S reached; NODE_6-like among them: {bool(good)})")
    elif good:
        print("READING: LINK")
    else:
        print("READING: OTHER (the 23S reached is not NODE_6-like)")


def cluster(names, regions, tmp):
    """Group 23S regions that are >= 99% identical over >= 1000 nt."""
    if len(names) < 2:
        return [names]
    p = os.path.join(tmp, "reach.fa")
    write(p, ((n, regions[n][6]) for n in names))
    par = {n: n for n in names}

    def root(x):
        while par[x] != x:
            x = par[x]
        return x
    for r in blast(p, p):
        if r[0] != r[1] and float(r[2]) >= 99 and int(r[3]) >= 1000:
            par[root(r[0])] = root(r[1])
    g = collections.defaultdict(list)
    for n in names:
        g[root(n)].append(n)
    return list(g.values())


def walk(sp, anchor, regions, tmp):
    gfa = os.path.join(sp, "assembly_graph_with_scaffolds.gfa")
    seg, links = {}, collections.defaultdict(list)
    for l in open(gfa):
        f = l.rstrip("\n").split("\t")
        if f[0] == "S":
            seg[f[1]] = f[2]
        elif f[0] == "L":
            links[(f[1], f[2])].append((f[3], f[4]))
            links[(f[3], "-" if f[4] == "+" else "+")].append((f[1], "-" if f[2] == "+" else "+"))
    sp_fa = os.path.join(tmp, "segs.fa")
    write(sp_fa, ((k, v) for k, v in seg.items() if len(v) >= 100))
    ah = [h for h in blast(anchor, sp_fa) if float(h[2]) >= 99]
    a_len = len(next(iter(fasta(anchor).values())))
    start = [(h[1], "+" if int(h[6]) < int(h[7]) else "-") for h in ah if int(h[5]) >= a_len - 20]
    print("GFA segments holding the anchor's 3' end:", start)
    seen, todo, out = set(), [(s, 0) for s in start], []
    while todo:
        (s, o), d = todo.pop()
        if (s, o) in seen or d > 6000:
            continue
        seen.add((s, o))
        out.append(s)
        for nxt in links.get((s, o), []):
            todo.append((nxt, d + len(seg[s])))
    rp = os.path.join(tmp, "reach_segs.fa")
    write(rp, ((s, seg[s]) for s in set(out) if len(seg[s]) >= 100))
    lp = os.path.join(tmp, "lsu_all.fa")
    write(lp, ((n, r[6]) for n, r in regions.items()))
    got = set()
    for r in blast(rp, lp):
        if float(r[2]) >= 99 and int(r[3]) >= 0.95 * len(seg[r[0]]):
            got.add(r[1])
    print(f"GFA walk: {len(seen)} oriented segments within 6 kb")
    return got


if __name__ == "__main__":
    main(*sys.argv[1:])
