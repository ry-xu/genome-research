"""V-C: release factors encoded by viruses (PREREG.md, "Phages and viruses").

Same candidate filter and search as the GTDB scan (src/gtdb_scan.py): proteins of 280-480 aa containing GGQ, DIAMOND
vs tag-bacteria's markers.dmnd (RF1, RF2, aRF1 + decoys), e-value 1e-20, best hit; RF1/RF2/aRF1 best hits are
written in rf_hits.tsv format (acc, gene, class, pident, len, seq) for src/rf_motifs.py.

usage: python src/viral_rf.py faa <proteins.faa[.gz] | -> <source> <out_prefix>     proteins (MGV, RefSeq viral)
       python src/viral_rf.py fna <genomes.fna[.gz] | -> <source> <out_prefix>     genomes, Pyrodigal meta calls
  -> <out_prefix>_cand.faa (candidates), <out_prefix>_rf_hits.tsv
  acc = genome or contig id (MGV: protein id minus its _N suffix; RefSeq: the [organism] of the header;
  gv = GVMAGs V2: the genome before "|")
"""
import gzip
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MARKERS = ROOT.parent / "tag-bacteria" / "data" / "markers.dmnd"
RFKEEP = ("RF1", "RF2", "aRF1")
LO, HI = 280, 480


def fasta(fn):
    f = sys.stdin if fn == "-" else (gzip.open(fn, "rt") if fn.endswith(".gz") else open(fn))
    head, buf = None, []
    for line in f:
        if line.startswith(">"):
            if head is not None:
                yield head, "".join(buf)
            head, buf = line[1:].rstrip("\n"), []
        else:
            buf.append(line.strip())
    if head is not None:
        yield head, "".join(buf)


def acc_of(head, source):
    pid = head.split()[0]
    if source == "refseq":
        org = head.rsplit("[", 1)[-1].rstrip("]") if "[" in head else ""
        return org.replace(" ", "_").replace("\t", "_"), pid
    if source == "mgv":
        return pid.rsplit("_", 1)[0], pid
    if source == "gv":  # GVMAGs V2: <genome>|<contig>_<n>
        return pid.split("|")[0], pid
    return pid, pid


def proteins_from_genomes(fn, threads=3):
    """Pyrodigal meta calls, `threads` genomes at a time (find_genes releases the GIL; output order kept)."""
    from concurrent.futures import ThreadPoolExecutor
    import pyrodigal
    f = pyrodigal.GeneFinder(meta=True)

    def one(rec):
        acc = rec[0].split()[0]
        return [(acc, f"{acc}_{k}", g.translate().rstrip("*")) for k, g in
                enumerate(f.find_genes(rec[1].upper().encode()), 1)]
    import itertools
    it = fasta(fn)
    with ThreadPoolExecutor(threads) as ex:
        while chunk := list(itertools.islice(it, 300)):  # bounded: map() would submit the whole file at once
            for prots in ex.map(one, chunk):
                yield from prots


def main(mode, fn, source, prefix):
    cand = Path(f"{prefix}_cand.faa")
    n = nc = 0
    accs = {}
    with open(cand, "w") as out:
        it = proteins_from_genomes(fn) if mode == "fna" else \
            ((*acc_of(h, source), s.rstrip("*")) for h, s in fasta(fn))
        for acc, gid, p in it:
            n += 1
            if LO <= len(p) <= HI and "GGQ" in p:
                nc += 1
                accs[gid] = acc
                out.write(f">{gid}\n{p}\n")
    print(f"{n} proteins, {nc} candidates", flush=True)
    h = Path(f"{prefix}_diamond.tsv")
    subprocess.run(["diamond", "blastp", "-q", str(cand), "-d", str(MARKERS), "-o", str(h), "--threads", "2",
                    "--quiet", "--evalue", "1e-20", "--max-target-seqs", "1", "--outfmt", "6", "qseqid", "sseqid",
                    "pident", "length"], check=True)
    seqs = {gid: s for gid, s in ((hd.split()[0], s) for hd, s in fasta(str(cand)))}
    seen, kept = set(), 0
    with open(f"{prefix}_rf_hits.tsv", "w") as out:
        for line in open(h):
            qid, sid, pid, ln = line.rstrip("\n").split("\t")
            if qid in seen:
                continue
            seen.add(qid)
            kind = sid.split("|")[0]
            if kind in RFKEEP:
                kept += 1
                out.write(f"{accs[qid]}\t{qid}\t{kind}\t{pid}\t{ln}\t{seqs[qid]}\n")
    print(f"{len(seen)} candidates with a marker hit, {kept} RF1/RF2/aRF1", flush=True)


if __name__ == "__main__":
    main(*sys.argv[1:5])
