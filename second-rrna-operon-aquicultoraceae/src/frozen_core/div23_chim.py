"""Chimera and contig-edge check for flagged divergent 23S copies (after two co-assembly/assembly-edge kills:
GCA_025542975, RS_GCF_042651645.1).

usage: python3 -I src/div23_chim.py CLASS_OR_CHECK.tsv SCANDIR LIST > OUT.tsv
  CLASS_OR_CHECK.tsv  div23_class.py or div23_check.py output (needs columns acc, copy)
  SCANDIR             div23_scan output dir (copies.*.fna; headers acc|contig|from|to|strand|model span)
  LIST                acc<TAB>genome path relative to data/external/codon-atlas/ (the scan's list)
For each flagged copy: blastn (task blastn) against every other 23S copy of the same genome; the copy with the most
identical bases is the reference. blastn joins a 98% stretch and an 88% stretch into one long HSP (the Streptomyces
case read 85.6% over 1658 nt), so identity is measured in 300-column windows (step 50) along each HSP >= 300 nt.
chimera = a run of consecutive windows all >= 97% spanning >= 600 columns AND some window <= 90% against the same
reference (a normal piece joined to a foreign piece). A single 98% window is not enough: the 88% half of the
Streptomyces copy reached 99% in one window (conserved domain V), so a real ~90% copy could too. Real divergent copies
(Crystal Geyser Aquicultor_A, DASXVE01) peak at 80-90% windows; DASXVE010000021 also sits at a contig edge, so edge
alone is a flag, not a kill. Calibration (300 random same-family different-genus census pairs, seed 7, 14 Nov): a
>=600-column run with a <=90% window occurred in 0/134 pairs below 0.90 overall identity, 1/55 at 0.90-0.92, 26/75 above.
So chimera = "yes" only when the copy's identity to the reference over HSPs >= 300 nt is < 0.90, "maybe" at 0.90-0.95,
"no" otherwise; a "yes" is a reason to look at the locus, not an automatic kill. segments print qstart-qend:pident[min-max window]/longest high run in columns. edge = the locus ends within 100 nt of a contig end.
Prints the input row's acc and copy, then: ref copy, segments (qstart-qend:pident), chimera, edge, contig length.
"""
import collections
import glob
import gzip
import os
import subprocess
import sys
import tempfile

ATLAS = "data/external/codon-atlas/"


def read_copies(scandir):
    seqs, name = {}, None
    for p in glob.glob(f"{scandir}/copies.*.fna"):
        for line in open(p):
            if line[0] == ">":
                name = line[1:].split()[0]
                seqs[name] = []
            else:
                seqs[name].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def win(qa, sa, w=300, step=50):
    """Percent identity of each w-column window of a gapped alignment."""
    same = [a == b for a, b in zip(qa, sa)]
    return [100 * sum(same[i:i + w]) / w for i in range(0, len(same) - w + 1, step)]


def run(ws, hi=97, w=300, step=50):
    """Columns spanned by the longest run of consecutive windows all >= hi."""
    best = cur = 0
    for x in ws:
        cur = cur + 1 if x >= hi else 0
        best = max(best, cur)
    return w + (best - 1) * step if best else 0


def contig_len(path, contig):
    n, on = 0, False
    for line in gzip.open(path, "rt"):
        if line[0] == ">":
            if on:
                break
            on = line[1:].split()[0] == contig
        elif on:
            n += len(line.strip())
    return n


def main(tab, scandir, listp):
    paths = dict(l.rstrip("\n").split("\t")[:2] for l in open(listp))
    seqs = read_copies(scandir)
    by_acc = collections.defaultdict(list)
    for k in seqs:
        by_acc[k.split("|")[0]].append(k)
    rows = [l.rstrip("\n").split("\t") for l in open(tab)]
    hdr = rows[0]
    ia, ic = hdr.index("acc"), hdr.index("copy")
    print("acc\tcopy\tref_copy\tsegments\tchimera\tedge\tcontig_len")
    with tempfile.TemporaryDirectory() as td:
        for r in rows[1:]:
            acc, cp = r[ia], r[ic]
            others = [k for k in by_acc[acc] if k != cp]
            q, s = os.path.join(td, "q.fa"), os.path.join(td, "s.fa")
            open(q, "w").write(f">q\n{seqs[cp]}\n")
            seg_txt, chim, ref = "-", "na", "-"
            if others:
                open(s, "w").write("".join(f">{i}\n{seqs[k]}\n" for i, k in enumerate(others)))
                out = subprocess.run(["blastn", "-task", "blastn", "-query", q, "-subject", s, "-outfmt",
                                      "6 sseqid pident length nident qstart qend qseq sseq", "-evalue", "1e-20"],
                                     capture_output=True, text=True, check=True).stdout
                hsps = collections.defaultdict(list)
                for line in out.splitlines():
                    sid, pid, ln, ni, qs, qe, qa, sa = line.split("\t")
                    hsps[int(sid)].append((int(qs), int(qe), float(pid), int(ln), int(ni), win(qa, sa)))
                if hsps:
                    best = max(hsps, key=lambda i: sum(h[4] for h in hsps[i]))
                    ref = others[best]
                    segs = sorted(h for h in hsps[best] if h[3] >= 300)
                    seg_txt = ",".join(f"{a}-{b}:{p:.1f}[{min(w, default=0):.0f}-{max(w, default=0):.0f}]/{run(w)}"
                                       for a, b, p, _, _, w in segs) or "-"
                    hi = any(run(h[5]) >= 600 for h in segs)
                    lo = any(x <= 90 for h in segs for x in h[5])
                    ident = sum(h[4] for h in segs) / max(1, sum(h[3] for h in segs))
                    chim = "no" if not (hi and lo) else "yes" if ident < 0.90 else "maybe" if ident < 0.95 else "no"
            f = cp.split("|")
            clen = contig_len(ATLAS + paths[acc], f[1]) if acc in paths else 0
            lo_, hi_ = sorted((int(f[2]), int(f[3])))
            edge = "yes" if clen and (lo_ <= 100 or clen - hi_ <= 100) else ("no" if clen else "na")
            print(f"{acc}\t{cp}\t{ref}\t{seg_txt}\t{chim}\t{edge}\t{clen}")


if __name__ == "__main__":
    main(*sys.argv[1:4])
