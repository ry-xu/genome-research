"""Place flagged genomes in GTDB R232 without GTDB-Tk: the ribosomal-protein markers of every species representative
(bac120: 36, ar53: 26 incl. rpo*/ksgA, from GTDB's per-marker MSA files, gaps stripped) become DIAMOND subjects named <rep>|<marker>.
A genome's proteins (pyrodigal, table 4 so TGA genomes are not chopped) are searched against them. Per genome and rep:
the best bitscore per marker, summed over markers. The top rep gives the GTDB taxonomy; the number of markers and their
median identity say how close it is. Rough guide: > 80 % same genus, > 60 % same family/order, < 50 % phylum only.
usage: uv run --no-project --with pyrodigal python -I src/place.py <acc>...   (genomes from data/genomes or data/sym_genomes)
writes results/placement.tsv (appends; one row per genome)"""
import gzip
import statistics
import subprocess
import sys
import tarfile
from collections import defaultdict
from pathlib import Path

import pyrodigal

HERE = Path(__file__).resolve().parent.parent
DL = HERE / "data" / "gtdb_markers_dl"
INFO = HERE / "data" / "gtdb_marker_info_dl"
DB = HERE / "data" / "gtdb_place"
META = HERE.parent.parent / "shared" / "codon-atlas" / "ref"
THREADS = "3"


def read_fasta(path):
    name, seq = None, []
    with (gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)) as f:
        for line in f:
            if line.startswith(">"):
                if name:
                    yield name, "".join(seq)
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line.strip())
    if name:
        yield name, "".join(seq)


def read_fasta_text(text):
    name, seq = None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if name:
                yield name, "".join(seq)
            name, seq = line[1:].split()[0], []
        else:
            seq.append(line.strip())
    if name:
        yield name, "".join(seq)


def ribosomal(dom):
    with open(INFO / f"{dom}_msa_marker_info_r232.tsv") as f:
        next(f)
        return {r.split("\t")[0] for r in f if "ibosomal" in r or "_bact\t" in r or "rp" in r.split("\t")[1][:2]}


def build_db():
    DB.mkdir(parents=True, exist_ok=True)
    if (DB / "ribo.dmnd").exists():
        return
    n = defaultdict(int)
    with open(DB / "ribo.faa", "w") as out:
        for dom in ("bac120", "ar53"):
            want = ribosomal(dom)
            with tarfile.open(DL / f"{dom}_msa_marker_genes_reps_r232.tar.gz", "r|gz") as tar:
                for m in tar:
                    marker = next((w for w in want if m.name.endswith(f"_reps_{w.split('_', 1)[1]}.faa")), None)
                    if not m.isfile() or marker is None:
                        continue
                    for name, s in read_fasta_text(tar.extractfile(m).read().decode()):
                        s = s.replace("-", "").replace("*", "").upper()
                        if len(s) >= 40:
                            out.write(f">{name}|{marker}\n{s}\n")
                            n[marker] += 1
    print(f"{len(n)} markers, {sum(n.values())} sequences", file=sys.stderr)
    subprocess.run(["diamond", "makedb", "--in", str(DB / "ribo.faa"), "-d", str(DB / "ribo"), "--threads", THREADS],
                   check=True)
    (DB / "ribo.faa").unlink()


def taxonomy():
    tax = {}
    for fn in ("bac120_metadata.tsv.gz", "ar53_metadata.tsv.gz"):
        with gzip.open(META / fn, "rt") as f:
            head = f.readline().rstrip("\n").split("\t")
            ia, it, ir = head.index("accession"), head.index("gtdb_taxonomy"), head.index("gtdb_representative")
            for line in f:
                r = line.rstrip("\n").split("\t")
                if r[ir] == "t":
                    tax[r[ia]] = r[it]
    return tax


def proteins(acc, out):
    path = next(p for p in (HERE / "data" / "genomes" / f"{acc}.fna.gz", HERE / "data" / "sym_genomes" / f"{acc}.fna.gz")
                if p.exists())
    contigs = [s.upper().encode() for _, s in read_fasta(path)]
    f = pyrodigal.GeneFinder(meta=sum(map(len, contigs)) < 100_000)
    if not f.meta:
        f.train(*contigs, translation_table=4)
    for i, c in enumerate(contigs):
        for j, g in enumerate(f.find_genes(c)):
            p = g.translate(translation_table=4).rstrip("*")
            if len(p) >= 40:
                out.write(f">{acc}#{i}_{j}\n{p}\n")


def main():
    accs = sys.argv[1:]
    build_db()
    q = DB / "query.faa"
    with open(q, "w") as out:
        for a in accs:
            proteins(a, out)
    hits = DB / "query.m8"
    subprocess.run(["diamond", "blastp", "-q", str(q), "-d", str(DB / "ribo"), "-o", str(hits), "--threads", THREADS,
                    "-k", "50", "-e", "1e-10", "-b", "1", "-c", "1", "--quiet",
                    "--outfmt", "6", "qseqid", "sseqid", "pident", "bitscore"], check=True)
    best = defaultdict(dict)  # (acc, rep) -> marker -> (bits, pident)
    for line in open(hits):
        qid, sid, pid, bits = line.rstrip("\n").split("\t")
        rep, marker = sid.split("|")
        k = (qid.split("#")[0], rep)
        if float(bits) > best[k].get(marker, (0, 0))[0]:
            best[k][marker] = (float(bits), float(pid))
    per = defaultdict(dict)
    for (acc, rep), d in best.items():
        per[acc][rep] = d
    tax = taxonomy()
    new = not (HERE / "results" / "placement.tsv").exists()
    with open(HERE / "results" / "placement.tsv", "a") as out:
        if new:
            out.write("acc\tnearest_rep\tn_markers\tmedian_pident\tgtdb_taxonomy\trunner_up_other_genus\n")
        for a in accs:
            reps = per.get(a, {})
            ranked = sorted(reps, key=lambda r: -sum(b for b, _ in reps[r].values()))
            if not ranked or len(reps[ranked[0]]) < 5:
                out.write(f"{a}\t\t{len(reps[ranked[0]]) if ranked else 0}\t\tunplaced\t\n")
                continue
            r0 = ranked[0]
            genus = tax.get(r0, "").rsplit(";", 1)[0]
            alt = next((tax.get(r, r) for r in ranked[1:] if tax.get(r, "").rsplit(";", 1)[0] != genus), "")
            out.write(f"{a}\t{r0}\t{len(reps[r0])}\t{statistics.median(p for _, p in reps[r0].values()):.1f}\t"
                      f"{tax.get(r0, '?')}\t{alt}\n")
            out.flush()


if __name__ == "__main__":
    main()
