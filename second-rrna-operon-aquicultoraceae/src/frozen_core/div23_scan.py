"""Find every 23S copy in cached genomes, to look for highly divergent second copies (paralogs).

usage: python -I src/div23_scan.py LIST HMM OUTDIR LANE NLANES
  LIST   acc<TAB>path (path relative to data/external/codon-atlas)
  HMM    Rfam LSU HMM (RF02541 bacteria, RF02540 archaea)
Takes every NLANES-th genome starting at LANE. Per genome: nhmmer (1 cpu) on the decompressed genome; every hit with
score >= 200 goes to OUTDIR/hits.LANE.tsv (acc contig from to strand hmm_from hmm_to score); hits spanning >= 1200
model positions are written (rRNA strand) to OUTDIR/copies.LANE.fna as >acc|contig|from|to|strand|hmmF-hmmT.
OUTDIR/done.LANE lists finished accessions, so a rerun resumes.
"""
import gzip
import os
import subprocess
import sys

ROOT = "data/external/codon-atlas/"
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
IUPAC = str.maketrans("RYKMSWBDHVrykmswbdhv", "N" * 20)  # nhmmer cannot guess the alphabet of IUPAC-heavy MAGs


def read_fasta_gz(path):
    out, name, buf = {}, None, []
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                if name:
                    out[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name:
        out[name] = "".join(buf)
    return out


def main(lst, hmm, outdir, lane, nlanes):
    done_path = f"{outdir}/done.{lane}"
    done = set(open(done_path).read().split()) if os.path.exists(done_path) else set()
    tmp = f"{outdir}/tmp.{lane}.fna"
    tbl = f"{outdir}/tmp.{lane}.tbl"
    rows = [l.rstrip("\n").split("\t") for l in open(lst)]
    with open(f"{outdir}/hits.{lane}.tsv", "a") as fh, open(f"{outdir}/copies.{lane}.fna", "a") as fc, \
            open(done_path, "a") as fd:
        for i, (acc, path) in enumerate(rows):
            if i % nlanes != lane or acc in done:
                continue
            try:
                g = read_fasta_gz(ROOT + path)
            except (OSError, EOFError) as e:
                print(f"{acc}\tread error {e}", file=sys.stderr)
                continue
            with open(tmp, "w") as f:
                for k, v in g.items():
                    f.write(f">{k}\n{v.translate(IUPAC)}\n")
            r = subprocess.run(["nhmmer", "--cpu", "1", "--noali", "--tblout", tbl, hmm, tmp],
                               stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
            if r.returncode:
                print(f"{acc}\tnhmmer error {r.stderr.strip()}", file=sys.stderr)
                continue
            for line in open(tbl):
                if line.startswith("#"):
                    continue
                x = line.split()
                ctg, hf, ht, af, at, strand, score = x[0], int(x[4]), int(x[5]), int(x[6]), int(x[7]), x[11], float(x[13])
                if score < 200:
                    continue
                lo, hi = min(af, at), max(af, at)
                fh.write(f"{acc}\t{ctg}\t{lo}\t{hi}\t{strand}\t{hf}\t{ht}\t{score}\n")
                if ht - hf + 1 >= 1200:
                    s = g[ctg][lo - 1:hi]
                    if strand == "-":
                        s = s.translate(COMP)[::-1]
                    fc.write(f">{acc}|{ctg}|{lo}|{hi}|{strand}|{hf}-{ht}\n{s}\n")
            fd.write(acc + "\n")
            fh.flush(); fc.flush(); fd.flush()
    for p in (tmp, tbl):
        if os.path.exists(p):
            os.remove(p)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5]))
