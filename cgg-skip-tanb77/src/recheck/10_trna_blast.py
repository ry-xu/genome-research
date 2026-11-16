"""tRNA cross-check: are there Arg tRNAs that the trna.parquet tables missed (G2 lists only ACG), and is there any
tRNA with a CGG-reading anticodon (CCG, or TCG/U34) in any genome? Query = every tRNA (all amino acids) from the 9
parquet tables; blastn against each genome; each genome locus is assigned to its best-scoring query; the locus
anticodon is read through the alignment at the query's anticodon position. Output: trna/LABEL.hits.tsv and a summary."""
import sys, subprocess, duckdb
from collections import defaultdict, Counter
sys.path.insert(0, "src/recheck")
from mp import read_fasta
W = "recheck_work"
G = dict(G1="GCA_946496695.1", G2="GCF_012317385.1", G3="GCA_900553755.1", G4="GCA_963931445.1", G5="GCF_005473905.2",
         G6="GCA_964631135.1", G7="GCA_009785865.1", G8="GCF_900187165.1", G9="GCA_949283575.1")
q = {}; qac = {}
for lab, acc in G.items():
    rows = duckdb.sql(f"select contig, begin, aa, anticodon, seq from 'data/atlas/products/{acc[:8]}/{acc}/trna.parquet'").fetchall()
    for (c, b, aa, ac, s) in rows:
        if aa in ("???",) or not s or len(ac) != 3: continue
        s = s.upper()
        pos = [i for i in range(26, 44) if s[i:i + 3] == ac]
        if not pos: continue
        nm = f"{lab}_{aa}_{ac}_{b}"
        q[nm] = s; qac[nm] = pos[0]
with open(f"{W}/trna/queries.fna", "w") as f:
    for k, v in q.items(): f.write(f">{k}\n{v}\n")
summary = {}
for lab in G:
    db = f"{W}/genomes/{lab}.fna"
    out = subprocess.run(["blastn", "-task", "blastn", "-word_size", "7", "-evalue", "1e-8", "-num_threads", "2",
                          "-query", f"{W}/trna/queries.fna", "-subject", db,
                          "-outfmt", "6 qseqid sseqid sstart send bitscore qstart qend qseq sseq pident length"],
                         capture_output=True, text=True).stdout
    best = {}
    for line in out.splitlines():
        qn, sn, ss, se, bs, qs, qe, qa, sa, pid, ln = line.split("\t")
        ss, se, bs, qs, qe = int(ss), int(se), float(bs), int(qs), int(qe)
        a = qac[qn]  # 0-based anticodon start in query
        if not (qs - 1 <= a and a + 3 <= qe): continue
        # walk alignment to anticodon
        qi = qs - 1; ac = ""
        for x, y in zip(qa, sa):
            if x != "-":
                if a <= qi < a + 3: ac += y
                qi += 1
            elif a < qi < a + 3:
                ac += y
        lo, hi = min(ss, se), max(ss, se)
        key = (sn, lo // 60)
        if key not in best or bs > best[key][0]:
            best[key] = (bs, qn, ac.upper(), lo, hi, float(pid), int(ln))
    with open(f"{W}/trna/{lab}.hits.tsv", "w") as f:
        f.write("contig\tlo\thi\tbest_query\tbitscore\tpident\tlen\tsubject_anticodon\n")
        for k, (bs, qn, ac, lo, hi, pid, ln) in sorted(best.items(), key=lambda x: (x[0][0], x[1][3])):
            f.write(f"{k[0]}\t{lo}\t{hi}\t{qn}\t{bs}\t{pid}\t{ln}\t{ac}\n")
    arg = Counter(ac for (bs, qn, ac, lo, hi, pid, ln) in best.values() if qn.split("_")[1] == "Arg" and ln >= 60)
    cggreaders = [(qn, ac, lo, pid, ln) for (bs, qn, ac, lo, hi, pid, ln) in best.values() if ac in ("CCG", "TCG") and ln >= 60]
    summary[lab] = (dict(arg), cggreaders)
    print(lab, "loci best-matching an Arg tRNA, by locus anticodon:", dict(arg), "| loci with CCG/TCG anticodon:", cggreaders, flush=True)
