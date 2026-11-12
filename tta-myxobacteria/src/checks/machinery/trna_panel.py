"""tRNA-Leu anticodons in the 16 blind-panel genomes (post hoc, machinery; reported only).

ARAGORN (pyaragorn) on every contig, both strands, table 11. Counts tRNA-Leu by anticodon; the TTA reader is
tRNA-Leu(UAA) (anticodon TAA). Also tmRNA count. Roles come from results/blind_panel_key.tsv.

usage: uv run --no-project --with pyaragorn,biopython python -I src/checks/machinery/trna_panel.py
"""
import collections
import csv
from pathlib import Path

import pyaragorn
from Bio import SeqIO

ROOT = Path(__file__).resolve().parents[3]


def main():
    key = list(csv.DictReader(open(ROOT / "results" / "blind_panel_key.tsv"), delimiter="\t"))
    finder = pyaragorn.RNAFinder(translation_table=11, tmrna=True, trna=True)
    out = []
    for k in sorted(key, key=lambda r: (r["role"], r["genus"], r["acc"])):
        leu, other, tm = collections.Counter(), 0, 0
        for rec in SeqIO.parse(ROOT / "blind_panel" / f"{k['id']}.fna", "fasta"):
            for g in finder.find_rna(bytes(rec.seq)):
                if isinstance(g, pyaragorn.TMRNAGene):
                    tm += 1
                elif g.amino_acid == "Leu":
                    leu[g.anticodon.upper().replace("U", "T")] += 1
                else:
                    other += 1
        row = [k["id"], k["role"], k["genus"][3:], k["acc"], k["rf2"], k["contigs"], str(sum(leu.values()) + other),
               str(tm), str(leu.get("TAA", 0)), ",".join(f"{a}:{n}" for a, n in sorted(leu.items()))]
        out.append("\t".join(row))
        print(out[-1], flush=True)
    hdr = "id\trole\tgenus\tacc\trf2\tcontigs\ttRNAs\ttmRNA\tLeu_TAA\tLeu_anticodons"
    (ROOT / "scratch" / "trna_panel" / "trna_panel.tsv").write_text(hdr + "\n" + "\n".join(out) + "\n")


if __name__ == "__main__":
    main()
