"""Clade-vs-sister fixed differences along RF1 and RF2 (post hoc; is L201F next to SPF special or one of many?).

Per genome, the top tblastn hit of E. coli RF1/RF2 mapped to E. coli coordinates. A position is a fixed difference if
all 8 clade genomes share one residue that no sister carries (all 8 sisters aligned there).

usage: uv run --no-project python -I src/checks/machinery/rf_fixed_diffs.py
"""
import csv
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import rf_panel  # noqa: E402

ROOT = rf_panel.ROOT


def mapped(qf, fna):
    r = subprocess.run(["tblastn", "-query", str(qf), "-subject", str(fna), "-evalue", "1e-10", "-outfmt", rf_panel.FMT,
                        "-seg", "no"], capture_output=True, text=True, check=True)
    h = max((dict(zip(rf_panel.FMT.split()[1:], l.split("\t"))) for l in r.stdout.splitlines()),
            key=lambda h: float(h["bitscore"]))
    m, qp = {}, int(h["qstart"])
    for a, b in zip(h["qseq"], h["sseq"]):
        if a != "-":
            m[qp] = b
            qp += 1
    return m


def main():
    faa = rf_panel.read_faa(rf_panel.MARKERS)
    key = list(csv.DictReader(open(ROOT / "results" / "blind_panel_key.tsv"), delimiter="\t"))
    with tempfile.TemporaryDirectory() as td:
        for name in ("RF1|RF1_ECOLI", "RF2|RF2_ECOLI"):
            ref = faa[name]
            qf = Path(td) / "q.faa"
            qf.write_text(f">q\n{ref}\n")
            maps = {k["id"]: (k["role"], mapped(qf, ROOT / "blind_panel" / f"{k['id']}.fna")) for k in key}
            fixed, shared = [], 0
            for p in range(1, len(ref) + 1):
                c = [m.get(p) for r, m in maps.values() if r == "clade"]
                s = [m.get(p) for r, m in maps.values() if r == "control"]
                if None in c or None in s or "-" in c + s:
                    continue
                shared += 1
                if len(set(c)) == 1 and c[0] not in s:
                    fixed.append(f"{ref[p - 1]}{p}{c[0]} (sisters {''.join(sorted(set(s)))})")
            print(f"{name.split('|')[0]}: {len(ref)} residues, {shared} aligned in all 16, {len(fixed)} fixed clade "
                  f"differences: " + ", ".join(fixed), flush=True)


if __name__ == "__main__":
    main()
