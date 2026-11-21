"""Build reference protein MSAs (MAFFT L-INS-i) per gene from table-4 translations.

usage: build_ref_msa.py refs_dir msa_dir [threads]
"""
import os
import sys
import subprocess
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import read_fasta, write_fasta  # noqa: E402

GENES = ["cox1", "cox2", "cox3", "cob", "atp6", "atp8", "atp9",
         "nad1", "nad2", "nad3", "nad4", "nad4L", "nad5", "nad6"]


def build(refs_dir, msa_dir, g):
    recs = [(n, s.replace("*", "X")) for n, s in read_fasta(os.path.join(refs_dir, f"ref_{g}.t4.faa"))]
    if not recs:
        return
    tmp = os.path.join(msa_dir, f"ref_{g}.in.faa")
    write_fasta(tmp, recs)
    out = os.path.join(msa_dir, f"ref_{g}.msa")
    with open(out, "w") as fh, open(out + ".log", "w") as lg:
        subprocess.run(["mafft", "--thread", "1", "--localpair", "--maxiterate", "1000", "--amino", tmp],
                       check=True, stdout=fh, stderr=lg)


def main():
    refs_dir, msa_dir = sys.argv[1:3]
    threads = int(sys.argv[3]) if len(sys.argv) > 3 else 4
    os.makedirs(msa_dir, exist_ok=True)
    with ThreadPoolExecutor(threads) as ex:
        list(ex.map(lambda g: build(refs_dir, msa_dir, g), GENES))


if __name__ == "__main__":
    main()
