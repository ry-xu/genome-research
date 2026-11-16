"""Follow-up S stage 1 (rule: results/rna_cgg/README.md): strandedness probe. Streams the first N spots of each run
with mate numbers (fastq-dump --split-spot --readids), maps like rna_stream.py, and reports, among filtered primary
CDS reads (rna_cgg.screen's definition) from read 1, the share on the CDS's sense strand; plus the ribosomal-protein
share of all filtered primary CDS reads.

usage: python -I src/rna_strand.py stream <outroot> <n_spots> <run>
       python -I src/rna_strand.py share <outroot> <run> [run ...]
"""
import bisect
import csv
import os
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import rna_cgg as R  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
# RNA_REF / RNA_MDIR / RNA_GTSV override the pig inputs (chicken replication, results/rna_chicken/README.md)
REF = Path(os.environ.get("RNA_REF", ROOT / "data" / "rna_cgg" / "ref" / "tanb77_24.mmi"))
MDIR = Path(os.environ.get("RNA_MDIR", ROOT / "data" / "proteomics" / "mgnify_pig"))
GTSV = Path(os.environ.get("RNA_GTSV", ROOT / "results" / "proteomics_pig" / "genomes.tsv"))
RP = re.compile(r"product=(30S|50S) ribosomal protein ")


def stream(outroot, n, run):
    out = Path(outroot) / run
    out.mkdir(parents=True, exist_ok=True)
    cmd = (f"set -o pipefail; fastq-dump --stdout --split-spot --skip-technical --readids -N 1 -X {n} {run}"
           f" 2>{out}/part.err | minimap2 -x sr -a -t 1 --sam-hit-only {REF} - 2>/dev/null"
           f" | awk -F'\\t' '$1 ~ /^@/ || ($5 >= {R.MAPQ} && int($2/4)%2 == 0 && int($2/256)%2 == 0"
           f" && int($2/2048)%2 == 0)' > {out}/part.sam")
    rc = subprocess.run(["bash", "-c", cmd]).returncode
    print(run, "exit", rc, flush=True)


def cds_strands():
    """per contig: sorted (start0, end, strand, is ribosomal protein) of CDS features"""
    out = defaultdict(list)
    for g in MDIR.glob("*.gff"):
        for line in open(g):
            if line.startswith("##FASTA"):
                break
            f = line.split("\t")
            if line[0] != "#" and len(f) >= 9 and f[2] == "CDS":
                out[f[0]].append((int(f[3]) - 1, int(f[4]), f[6], bool(RP.search(f[8]))))
    for v in out.values():
        v.sort()
    return out


def share(outroot, runs):
    prim = {g["genome"] for g in csv.DictReader(open(GTSV), delimiter="\t") if g["primary"] == "1"}
    nc = R.features(MDIR)["ncRNA"]
    cds = cds_strands()
    starts = {c: [x[0] for x in v] for c, v in cds.items()}
    print("run\tspots_streamed\tcds_reads\trp_reads\trp_share\tr1_cds_reads\tr1_sense\tr1_sense_share\tcall")
    for run in runs:
        d = Path(outroot) / run
        spots = sum(int(m) for m in re.findall(r"^Read (\d+) spots", (d / "part.err").read_text(), re.M))
        n = nrp = n1 = s1 = 0
        for line in open(d / "part.sam"):
            if line[0] == "@":
                continue
            f = line.split("\t")
            a = R.aln(f)
            if a is None or a[0].rsplit("_", 1)[0] not in prim:
                continue
            c, a0, a1, _ = a
            if c not in cds or R.overlaps(nc.get(c, []), a0, a1):
                continue
            v, st = cds[c], starts[c]
            i = bisect.bisect_left(st, a1)
            hit = [v[j] for j in range(i - 1, max(i - 50, 0) - 1, -1) if v[j][0] < a1 and v[j][1] > a0]
            if not hit:
                continue
            n += 1
            nrp += any(h[3] for h in hit)
            strands = {h[2] for h in hit}
            if len(strands) == 1 and f[0].endswith(".1"):
                n1 += 1
                s1 += ("-" if int(f[1]) & 16 else "+") == strands.pop()
        sh = s1 / n1 if n1 else float("nan")
        call = "stranded" if n1 >= 200 and (sh <= 0.10 or sh >= 0.90) else ("too few" if n1 < 200 else "unstranded")
        print(f"{run}\t{spots}\t{n}\t{nrp}\t{nrp / n if n else float('nan'):.5f}\t{n1}\t{s1}\t{sh:.4f}\t{call}",
              flush=True)


if __name__ == "__main__":
    if sys.argv[1] == "stream":
        stream(sys.argv[2], int(sys.argv[3]), sys.argv[4])
    else:
        share(sys.argv[2], sys.argv[3:])
