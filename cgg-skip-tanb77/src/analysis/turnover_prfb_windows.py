"""prfB positive control for site_turnover: 121-nt windows centred on the TGA of each CTT-TGA +1 break on an RF2 hit."""
import csv, gzip, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from candidate import read_fasta
COMP = str.maketrans("ACGTN", "TGCAN")
setf, outdir, gdir = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3])
for r in csv.DictReader(open(setf), delimiter="\t"):
    a = r["gca"]
    f = outdir / f"{a}.fs.tsv.gz"
    if not f.exists():
        continue
    seqs = None
    for x in csv.DictReader(gzip.open(f, "rt"), delimiter="\t"):
        if x["type"] != "+1" or "RF2_" not in x["sseqid"] or x["old_codon"] != "TGA":
            continue
        seqs = seqs or read_fasta(str(gdir / a[:8] / f"{a}.fna.gz"))
        s = seqs[x["contig"]] if x["strand"] == "+" else seqs[x["contig"]].translate(COMP)[::-1]
        c = x["ctx"].upper(); i = s.find(c)
        if i < 0 or s.find(c, i + 1) >= 0:
            continue
        p = i + sum(ch.islower() for ch in x["ctx"])
        w = s[p - 60:p + 61]
        if len(w) == 121 and s[p - 3:p + 3] == "CTTTGA":
            print(f">{a}|prfb|{x['contig']}|{x['strand']}|{p}\n{w}")
