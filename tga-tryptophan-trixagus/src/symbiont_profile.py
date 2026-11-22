"""Does a genome read TGA as Trp? Three signals that do not use read-through or Swiss-Prot:
  1. codon usage: genes called with table 4 (TGA = Trp); internal TGA per 1000 codons and the stop codons ending genes.
     A TGA = Trp genome has many internal TGA and few TGA ends; a standard genome has ~0 internal TGA.
  2. coding density / mean gene length under table 11 vs table 4 (a recoded genome's genes are chopped by table 11).
  3. tRNA-Trp anticodons (tRNAscan-SE -B): UCA reads UGA; and release factors (RF1 / RF2 hits vs shared rf.dmnd).
usage: uv run --no-project --with pyrodigal,requests python -I src/symbiont_profile.py <out.tsv> <acc[:label]>...
genomes are cached in data/sym_genomes/ (downloaded data, kept apart from code)"""
import gzip
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path

import pyrodigal

HERE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(HERE.parent / "stop-codon-census" / "src"))
import readthrough as rt  # noqa: E402

rt.GEN = HERE / "data" / "sym_genomes"
RFDB = HERE.parent.parent / "shared" / "codon-atlas" / "ref" / "rf.dmnd"
TRNASCAN = Path.home() / "dev" / "tools" / "bin" / "tRNAscan-SE"
COLS = ["acc", "label", "bp", "contigs", "gc", "genes4", "cds_frac11", "cds_frac4", "mean_aa11", "mean_aa4",
        "TGA_per_k", "TGG_per_k", "end_TAA", "end_TAG", "end_TGA", "trp_anticodons", "RF1", "RF2"]


def genes(contigs, table):
    total = sum(map(len, contigs))
    f = pyrodigal.GeneFinder(meta=total < 100_000)
    if not f.meta:
        f.train(*contigs, translation_table=table)
    return [g for c in contigs for g in f.find_genes(c)]


def one(acc, label):
    seqs = dict(rt.read_fasta(rt.fetch(acc)))
    contigs = [s.upper().encode() for s in seqs.values()]
    total = sum(map(len, contigs))
    gc = sum(c.count(b"G") + c.count(b"C") for c in contigs) / total
    g11, g4 = genes(contigs, 11), genes(contigs, 4)
    cod, ends = Counter(), Counter()
    for g in g4:
        nt = g.sequence()
        cods = [nt[i:i + 3] for i in range(0, len(nt) - 2, 3)]
        complete_end = not (g.partial_end if g.strand == 1 else g.partial_begin)
        if complete_end:
            ends[cods[-1]] += 1
            cods = cods[:-1]
        cod.update(cods[1:])
    n = sum(cod.values())
    with tempfile.TemporaryDirectory() as td:
        fa = Path(td) / "g.fna"
        fa.write_text("".join(f">{k.split()[0]}\n{v}\n" for k, v in seqs.items()))
        out = Path(td) / "t.txt"
        subprocess.run([str(TRNASCAN), "-B", "-q", "--brief", "-o", str(out), str(fa)], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        trp = sorted(l.split("\t")[5].strip() for l in open(out) if l.split("\t")[4].strip() == "Trp") if out.exists() else []
        q = Path(td) / "q.faa"
        q.write_text("".join(f">p{i}\n{g.translate(translation_table=4).rstrip('*')}\n" for i, g in enumerate(g4)))
        h = Path(td) / "h.tsv"
        subprocess.run(["diamond", "blastp", "-q", str(q), "-d", str(RFDB), "-o", str(h), "--threads", "1", "--quiet",
                        "--evalue", "1e-20", "--max-target-seqs", "1", "--outfmt", "6", "qseqid", "sseqid"], check=True)
        rf = Counter("RF1" if "|RF1_" in l.split("\t")[1] else "RF2" for l in open(h))
    return dict(acc=acc, label=label, bp=total, contigs=len(contigs), gc=round(gc, 3), genes4=len(g4),
                cds_frac11=round(sum(len(g.sequence()) for g in g11) / total, 3),
                cds_frac4=round(sum(len(g.sequence()) for g in g4) / total, 3),
                mean_aa11=round(sum(len(g.sequence()) for g in g11) / 3 / max(1, len(g11))),
                mean_aa4=round(sum(len(g.sequence()) for g in g4) / 3 / max(1, len(g4))),
                TGA_per_k=round(1000 * cod["TGA"] / n, 2), TGG_per_k=round(1000 * cod["TGG"] / n, 2),
                end_TAA=ends["TAA"], end_TAG=ends["TAG"], end_TGA=ends["TGA"], trp_anticodons=",".join(trp) or "-",
                RF1=rf["RF1"], RF2=rf["RF2"])


if __name__ == "__main__":
    out = Path(sys.argv[1])
    new = not out.exists()
    with open(out, "a") as f:
        if new:
            f.write("\t".join(COLS) + "\n")
        for a in sys.argv[2:]:
            acc, _, label = a.partition(":")
            try:
                r = one(acc, label)
            except Exception as e:  # noqa: BLE001
                print(acc, "ERROR", e, file=sys.stderr)
                continue
            f.write("\t".join(str(r[c]) for c in COLS) + "\n")
            f.flush()
            print("\t".join(str(r[c]) for c in COLS), flush=True)
