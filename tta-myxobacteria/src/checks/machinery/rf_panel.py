"""RF1/RF2 copies in the 16 blind-panel genomes from tblastn on raw contigs (no gene calls; post hoc, machinery).

Queries: E. coli RF1 and RF2 from the census markers.faa. Each genomic locus (contig, strand, overlapping span) goes to
the query with the higher bitscore. The motif is read off the tblastn alignment at E. coli RF1 188-190 (PxT) or
RF2 205-207 (SPF). Asks whether JAXLYQ01 (SPF in GTDB's gene calls, 59% complete) carries a second, SPY copy.

usage: uv run --no-project python -I src/checks/machinery/rf_panel.py
"""
import gzip
import csv
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
MARKERS = ROOT.parent / "tag-bacteria" / "data" / "reference" / "markers.faa.gz"
REFS = {"RF1|RF1_ECOLI": ("RF1", 188), "RF2|RF2_ECOLI": ("RF2", 205)}
FMT = "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qseq sseq"


def read_faa(fn):
    seqs, name = {}, None
    for line in (gzip.open(fn, "rt") if str(fn).endswith(".gz") else open(fn)):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        elif name:
            seqs[name].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def motif(h, start):
    q, s, qp = h["qseq"], h["sseq"], int(h["qstart"])
    out = {}
    for a, b in zip(q, s):
        if a != "-":
            out[qp] = b
            qp += 1
    return "".join(out.get(i, "-") for i in range(start, start + 3))


def main():
    faa = read_faa(MARKERS)
    key = list(csv.DictReader(open(ROOT / "results" / "blind_panel_key.tsv"), delimiter="\t"))
    with tempfile.TemporaryDirectory() as td:
        qf = Path(td) / "q.faa"
        qf.write_text("".join(f">{k}\n{faa[k]}\n" for k in REFS))
        rows = []
        for k in sorted(key, key=lambda r: (r["role"], r["genus"], r["acc"])):
            r = subprocess.run(["tblastn", "-query", str(qf), "-subject", str(ROOT / "blind_panel" / f"{k['id']}.fna"),
                                "-evalue", "1e-10", "-outfmt", FMT, "-seg", "no"], capture_output=True, text=True,
                               check=True)
            hits = [dict(zip(FMT.split()[1:], l.split("\t"))) for l in r.stdout.splitlines()]
            loci = []
            for h in sorted(hits, key=lambda h: -float(h["bitscore"])):
                st, en = sorted((int(h["sstart"]), int(h["send"])))
                strand = int(h["send"]) > int(h["sstart"])
                if any(l["sseqid"] == h["sseqid"] and l["strand"] == strand and st < l["en"] and en > l["st"]
                       for l in loci):
                    continue
                cls, mpos = REFS[h["qseqid"]]
                loci.append(dict(h, st=st, en=en, strand=strand, cls=cls,
                                 motif=motif(h, mpos) if int(h["qstart"]) <= mpos and int(h["qend"]) >= mpos + 2
                                 else "n/a"))
            desc = "; ".join(f"{l['cls']} {l['motif']} ({l['pident']}%, q{l['qstart']}-{l['qend']}, bits "
                             f"{float(l['bitscore']):.0f})" for l in sorted(loci, key=lambda l: l["cls"]))
            rows.append(f"{k['id']}\t{k['role']}\t{k['genus'][3:]}\t{k['acc']}\tGTDB-call RF2 {k['rf2']}\t{desc}")
            print(rows[-1], flush=True)
    (ROOT / "scratch" / "rf_panel" / "rf_panel.tsv").write_text("\n".join(rows) + "\n")


if __name__ == "__main__":
    main()
