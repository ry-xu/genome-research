"""Collect proteins within W bp of every rRNA operon unit (65 order genomes + controls).
Units: loci of out/loci.tsv grouped by operon id (unpaired loci alone), frac >= 0.3; control operons from out/ctl_*.fna ids.
Writes out/units.tsv (unit, acc, genus, label, contig, start, end, strand, contig_len, room_left, room_right)
       out/flank_genes.tsv (unit, gene_id, side, dist, rel_strand, aa_len)
       out/flank.faa
Side 'up' = 5' of the operon in its own orientation (before 16S), 'down' = 3' (after 23S/5S).
"""
import csv, glob, os
from collections import defaultdict
from Bio import SeqIO

W = 10000
FAA = "lab/work/frozen-core/results/lsu/bac/cca/kill_div/faa"
clen = {r.id: len(r.seq) for r in SeqIO.parse("data/aq65.fna", "fasta")}
for g in ["ec", "vv", "tb", "tt"]:
    for r in SeqIO.parse(f"data/ctl/{g}.fna", "fasta"):
        clen[f"{g}|{r.id}"] = len(r.seq)

loci = [r for r in csv.DictReader(open("out/loci.tsv"), delimiter="\t") if float(r["frac"]) >= 0.3]
units = defaultdict(list)
for r in loci:
    units[r["operon"] or r["id"]].append(r)

U = []
for u, rs in units.items():
    lab = "div" if any(r["label"] == "div" for r in rs) else "norm"
    U.append(dict(unit=u, acc=rs[0]["acc"], family=rs[0]["family"], genus=rs[0]["genus"], label=lab,
                  contig=rs[0]["contig"], start=min(int(r["start"]) for r in rs), end=max(int(r["end"]) for r in rs),
                  strand=rs[0]["strand"], members=",".join(r["id"] for r in rs)))

# controls: pair 16S and 23S on same contig/strand within 2 kb
ctl = defaultdict(list)
for fn, kind in [("out/ctl_ssu.fna", "16S"), ("out/ctl_lsu.fna", "23S")]:
    for r in SeqIO.parse(fn, "fasta"):
        org, _, rest = r.id.split("_", 2)
        ctg, pos = rest.rsplit("_", 1)
        st = pos[-1]; s = int(pos[:-1]); e = s + len(r.seq) - 1
        ctl[(org, ctg, st)].append((s, e, kind, r.id))
for (org, ctg, st), lst in ctl.items():
    lst.sort()
    grp = []
    for x in lst:
        if grp and x[0] - grp[-1][-1][1] < 2500:
            grp[-1].append(x)
        else:
            grp.append([x])
    for g in grp:
        U.append(dict(unit=f"{org}_{g[0][0]}", acc=org, family="control", genus=org, label="ctl", contig=f"{org}|{ctg}",
                      start=g[0][0], end=g[-1][1], strand=st, members=",".join(x[3] for x in g)))

prot = {}
def load(fn, prefix_acc):
    for r in SeqIO.parse(fn, "fasta"):
        p = r.id.split("|")
        acc, ctg, s, e, st = p[0], p[1], int(p[2]), int(p[3]), p[4]
        key = f"{acc}|{ctg}"
        prot.setdefault(key, []).append((s, e, st, r.id, str(r.seq).rstrip("*")))

accs = {u["acc"] for u in U}
for a in accs:
    fn = f"{FAA}/{a}.faa" if not a in ("ec", "vv", "tb", "tt") else f"out/ctl_faa/{a}.faa"
    if os.path.exists(fn):
        load(fn, a)

with open("out/units.tsv", "w") as fu, open("out/flank_genes.tsv", "w") as fg, open("out/flank.faa", "w") as ff:
    fu.write("unit\tacc\tfamily\tgenus\tlabel\tcontig\tstart\tend\tstrand\tcontig_len\troom_up\troom_down\tmembers\n")
    fg.write("unit\tgene_id\tside\tdist\trel_strand\taa_len\n")
    seen = set()
    for u in U:
        L = clen.get(u["contig"], 0)
        left, right = u["start"] - 1, L - u["end"]
        up, down = (left, right) if u["strand"] == "+" else (right, left)
        fu.write("\t".join(str(u[k]) for k in ["unit", "acc", "family", "genus", "label", "contig", "start", "end", "strand"])
                 + f"\t{L}\t{min(up, W)}\t{min(down, W)}\t{u['members']}\n")
        for s, e, st, gid, seq in prot.get(u["contig"], []):
            if e < u["start"] - W or s > u["end"] + W:
                continue
            if e < u["start"]:
                side_lr, dist = "L", u["start"] - e
            elif s > u["end"]:
                side_lr, dist = "R", s - u["end"]
            else:
                side_lr, dist = "IN", 0
            if side_lr == "IN":
                side = "in"
            else:
                side = ("up" if side_lr == "L" else "down") if u["strand"] == "+" else ("down" if side_lr == "L" else "up")
            rel = "same" if st == u["strand"] else "opp"
            fg.write(f"{u['unit']}\t{gid}\t{side}\t{dist}\t{rel}\t{len(seq)}\n")
            if gid not in seen:
                ff.write(f">{gid}\n{seq}\n"); seen.add(gid)
print(len(U), "units")
