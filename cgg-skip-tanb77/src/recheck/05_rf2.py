"""Positive control 1: the programmed +1 frameshift in prfB (RF2).
Full-length RF2 proteins (B. subtilis P28367, C. perfringens Q0TNE1, E. coli P07012) were aligned with miniprot -S.
For each genome and each RF2 query, report the best alignment, its frameshift events (query position, direction)
and the genome sequence around each event; the expected event is a +1 near query position 20-30 at CTT TGA C."""
import sys
sys.path.insert(0, "src/recheck")
from mp import read_fasta, parse_paf

W = "recheck_work"
prots = read_fasta(f"{W}/refs/rf2_full.faa")
print("label\tquery\tqspan\tqlen\tAS\tlocus\tevents(qpos:dir:context)")
for lab in sys.argv[1:]:
    genome = read_fasta(f"{W}/genomes/{lab}.fna")
    alns = list(parse_paf(f"{W}/out/{lab}.rf2.paf", genome, prots))
    if not alns:
        print(f"{lab}\t-\tno alignment"); continue
    for q in sorted(set(a["qn"] for a in alns)):
        a = max((x for x in alns if x["qn"] == q), key=lambda x: x["AS"])
        evs = []
        for (e, d, n, op), qp in zip(a["events"], a["evq"]):
            ctx = a["seg"][max(0, e - 9):e].lower() + "|" + a["seg"][e:e + n].upper() + "|" + a["seg"][e + n:e + n + 9].lower()
            evs.append(f"{qp}:{'+1' if d == 1 else '-1'}:{op}{n}:{ctx}")
        print(f"{lab}\t{q.split('|')[2]}\t{a['qs']}-{a['qe']}\t{a['qlen']}\t{a['AS']}\t{a['contig']}:{a['ts']}-{a['te']}{a['strand']}\t{' '.join(evs) or 'none'}")
