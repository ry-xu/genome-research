"""Giant-virus reference for Arm B (PREREG.md, Changes 2026-11-10 21:04): GVMAGs V2 release proteins of every genome
whose GVClass `ttable` is code1, code11 or codemeta (standard stops), minus the control genomes (null30.txt), the lead
genomes and the genomes whose eRF1 motif is missing (given on the command line) -> data/giantvirus/gvref.faa and
gvref.dmnd. Streams the release tar; protein ids are kept as released (<genome>|<contig>_<n>).

usage: python -I src/gv_ref.py <exclude genome> ...
"""
import subprocess
import sys
import tarfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
GV = ROOT / "data" / "giantvirus"
KEEP = ("code1", "code11", "codemeta")


def main(extra):
    rows = [line.rstrip("\n").split("\t") for line in open(GV / "GVMAGS_V2_STATS.tsv")]
    col = rows[0].index("ttable")
    exclude = set(extra) | set((GV / "null30.txt").read_text().split())
    want = {r[0] for r in rows[1:] if r[col] in KEEP} - exclude
    print(f"{len(want)} reference genomes; {len(exclude)} excluded", flush=True)
    seen = n = 0
    with open(GV / "gvref.faa", "w") as out, tarfile.open(GV / "GVMAGS_V2_database.gz", "r|gz") as tar:
        for m in tar:
            name = Path(m.name).name
            if not (m.isfile() and name.endswith(".faa") and name[:-4] in want):
                continue
            seen += 1
            for line in tar.extractfile(m).read().decode().splitlines():
                if line.startswith(">"):
                    n += 1
                    out.write(line.split()[0] + "\n")
                elif line:
                    out.write(line.rstrip("*") + "\n")
    print(f"{seen} genomes, {n} proteins", flush=True)
    if seen != len(want):
        sys.exit(f"{len(want) - seen} reference genomes not in the release tar")
    subprocess.run(["diamond", "makedb", "--in", str(GV / "gvref.faa"), "-d", str(GV / "gvref"), "--threads", "2",
                    "--quiet"], check=True)


if __name__ == "__main__":
    main(sys.argv[1:])
