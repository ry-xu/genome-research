"""RF2 residues at chosen E. coli positions in the panel's main RF2 copy (post hoc, machinery; reported only).

Positions: 149 and 170 (Uno, Ito & Nakamura 2002: E. coli RF2 E149K-E170K released at UUA in vitro), 205-207 (SPF),
and the window 195-215. Uses the top-scoring RF2 tblastn hit per genome, as in rf_panel.py.

usage: uv run --no-project python -I src/checks/machinery/rf2_sites.py
"""
import csv
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import rf_panel  # noqa: E402

ROOT = rf_panel.ROOT
POS = [149, 170] + list(range(195, 216))


def main():
    faa = rf_panel.read_faa(rf_panel.MARKERS)
    ref = faa["RF2|RF2_ECOLI"]
    print("E. coli", ref[148], ref[169], ref[194:215], sep="\t")
    key = list(csv.DictReader(open(ROOT / "results" / "blind_panel_key.tsv"), delimiter="\t"))
    with tempfile.TemporaryDirectory() as td:
        qf = Path(td) / "q.faa"
        qf.write_text(f">RF2\n{ref}\n")
        for k in sorted(key, key=lambda r: (r["role"], r["genus"], r["acc"])):
            r = subprocess.run(["tblastn", "-query", str(qf), "-subject", str(ROOT / "blind_panel" / f"{k['id']}.fna"),
                                "-evalue", "1e-10", "-outfmt", rf_panel.FMT, "-seg", "no"],
                               capture_output=True, text=True, check=True)
            h = max((dict(zip(rf_panel.FMT.split()[1:], l.split("\t"))) for l in r.stdout.splitlines()),
                    key=lambda h: float(h["bitscore"]))
            m, qp = {}, int(h["qstart"])
            for a, b in zip(h["qseq"], h["sseq"]):
                if a != "-":
                    m[qp] = b
                    qp += 1
            s = "".join(m.get(p, "-") for p in POS)
            print(k["id"], k["role"], k["genus"][3:], k["acc"], s[0], s[1], s[2:], sep="\t", flush=True)


if __name__ == "__main__":
    main()
