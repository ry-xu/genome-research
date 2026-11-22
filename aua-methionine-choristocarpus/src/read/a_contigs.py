"""Find mito contigs in a WGS assembly: dc-megablast NC_085326.1 against it (rules: checks/chor_wgs/RULES.md).
usage: uv run python -I checks/chor_wgs/a_contigs.py ASSEMBLY.fna QUERY.fa WORKDIR  -> WORKDIR/mito.fa + printed table"""
import os, subprocess, sys
from collections import defaultdict

fna, query, wd = sys.argv[1:4]
os.makedirs(wd, exist_ok=True)
db = os.path.join(wd, "asm")
if not os.path.exists(db + ".nsq") and not os.path.exists(db + ".00.nsq"):
    subprocess.run(["makeblastdb", "-in", fna, "-dbtype", "nucl", "-out", db, "-parse_seqids"], check=True, capture_output=True)
out = os.path.join(wd, "hits.tsv")
subprocess.run(["blastn", "-task", "dc-megablast", "-query", query, "-db", db, "-evalue", "1e-20", "-num_threads", "1",
                "-max_target_seqs", "500", "-outfmt", "6 sseqid pident length qstart qend sstart send slen evalue bitscore",
                "-out", out], check=True)
aln, ident, slen, qcov = defaultdict(int), defaultdict(float), {}, defaultdict(set)
for l in open(out):
    s, pid, L, qs, qe, ss, se, sl, ev, bs = l.split("\t")
    aln[s] += int(L); ident[s] += float(pid) * int(L); slen[s] = int(sl)
    qcov[s].update(range(min(int(qs), int(qe)), max(int(qs), int(qe)) + 1))
keep = [s for s in aln if aln[s] >= 2000 and ident[s] / aln[s] >= 70]
qlen = sum(len(l.strip()) for l in open(query) if not l.startswith(">"))
for s in sorted(aln, key=lambda s: -aln[s])[:15]:
    print(s, slen[s], aln[s], round(ident[s] / aln[s], 1), "MITO" if s in keep else "")
allq = set().union(*(qcov[s] for s in keep)) if keep else set()
print("mito contigs", len(keep), "covered share of query", round(len(allq) / qlen, 3))
subprocess.run(["blastdbcmd", "-db", db, "-entry", ",".join(keep), "-out", os.path.join(wd, "mito.fa")], check=True)
