"""UBA6777 read test on Tara run ERR599038 (the run GCA_002453135 / DKGF01 was binned from).

usage: python3 -I tara_ux.py HITS.sam OUT_PREFIX
HITS.sam: minimap2 -ax sr --sam-hit-only --secondary=no against ref3.fa = MAG|<DKGF01 contigs> + CEN|<census 23S, one per
genome, UBA6777 entries dropped> + CEN|FAM_<3 normal 23S loci of the family (JACKEL01, JACKEK01 x2)>. Identity = 1 - NM / (aligned query bases + deleted bases).
Reports:
  1. MAG coverage: reads per kb at >=0.98 for contigs >=5 kb (median = the host's coverage).
  2. Divergent 23S (MAG|DKGF01000210.1:6356-10856, minus strand) reads per kb at >=0.98, and reads spanning each end
     of the 23S (+-15 nt around 6356 and 10856) at >=0.98: if the junctions are read-supported, both are near the median.
  3. 23S reads elsewhere: reads on CEN at >=0.90, best genome taxon. Writes OUT_PREFIX.cen.names (every read name with any
     CEN hit; mates share names under --split-spot -Z), for extraction and assembly.
"""
import collections
import gzip
import re
import statistics
import sys

DIV = ("MAG|DKGF01000210.1", 6356, 10856)
# control positions on the same contig, outside the 23S-5S unit (contig length checked at run time)
CTRL = (2000, 3000, 4000, 5000, 12500, 14000, 15500, 17000, 18500, 20000)
TAX = "data/external/codon-atlas/ref/bac120_metadata.tsv.gz"


def ident(cigar, nm):
    m = sum(int(n) for n, op in re.findall(r"(\d+)([MID])", cigar) if op in "MID")
    return 1 - nm / m if m else 0.0, re.findall(r"(\d+)([MIDNSHP=X])", cigar)


def ref_end(pos, ops):
    return pos + sum(int(n) for n, op in ops if op in "MDN=X") - 1


def main(sam, out):
    clen, per_ctg, div, ends = {}, collections.Counter(), 0, collections.Counter()
    cen_reads, cen_best = set(), {}
    for line in open(sam):
        if line[0] == "@":
            if line.startswith("@SQ"):
                f = dict(x.split(":", 1) for x in line.split("\t")[1:])
                clen[f["SN"]] = int(f["LN"])
            continue
        f = line.split("\t", 11)
        flag, rname = int(f[1]), f[2]
        if flag & 4 or flag & 256 or flag & 2048:
            continue
        nm = int(re.search(r"NM:i:(\d+)", f[11]).group(1))
        idn, ops = ident(f[5], nm)
        pos = int(f[3])
        end = ref_end(pos, ops)
        if rname.startswith("MAG|"):
            if idn >= 0.98:
                per_ctg[rname] += 1
                if rname == DIV[0]:
                    if pos <= DIV[2] and end >= DIV[1]:
                        div += 1
                    for x in (DIV[1], DIV[2]) + CTRL:
                        if pos <= x - 15 and end >= x + 15:
                            ends[x] += 1
        else:
            cen_reads.add(f[0])  # every census hit goes to the assembly, whatever its identity
            acc = rname.split("|")[1]
            if idn >= 0.90 and (f[0] not in cen_best or idn > cen_best[f[0]][1]):
                cen_best[f[0]] = (acc, idn)
    dens = {c: per_ctg[c] / clen[c] * 1000 for c in clen if c.startswith("MAG|") and clen[c] >= 5000}
    med = statistics.median(dens.values())
    q = statistics.quantiles(dens.values(), n=10)
    print(f"MAG contigs >=5 kb: {len(dens)}; reads/kb at >=0.98: median {med:.1f}, deciles 10-90% {q[0]:.1f}-{q[-1]:.1f}")
    span = DIV[2] - DIV[1] + 1
    print(f"divergent 23S: {div} reads, {div / span * 1000:.1f}/kb ({div / span * 1000 / med:.2f}x median)")
    print(f"div contig {DIV[0]} overall: {dens.get(DIV[0], 0):.1f}/kb")
    print(f"reads spanning 23S ends (+-15 nt): {DIV[1]}: {ends[DIV[1]]}, {DIV[2]}: {ends[DIV[2]]}")
    ok = [x for x in CTRL if x + 15 < clen[DIV[0]]]
    print(f"same contig, control positions ({len(ok)} within contig): " + ", ".join(f"{x}: {ends[x]}" for x in ok) +
          f"; median {statistics.median(ends[x] for x in ok)}")
    need = {a for a, _ in cen_best.values()}
    tax = {}
    for l in gzip.open(TAX, "rt"):
        a = l.split("\t", 1)[0]
        if a in need:
            tax[a] = next(x for x in l.split("\t") if x.startswith("d__"))
    by = collections.Counter()
    hi = collections.Counter()
    for r, (a, idn) in cen_best.items():
        t = ";".join(tax.get(a, "?").split(";")[1:3])
        by[t] += 1
        if idn >= 0.97:
            hi[";".join(tax.get(a, "?").split(";")[1:6])] += 1
    print(f"reads on census 23S at >=0.90: {len(cen_best)}; by phylum/class (top 12):")
    for t, n in by.most_common(12):
        print(f"  {n}\t{t}")
    print("Myxococcota* at >=0.97 identity (family level):")
    for t, n in hi.most_common():
        if "Myxococcota" in t:
            print(f"  {n}\t{t}")
    with open(out + ".cen.names", "w") as o:
        o.write("\n".join(sorted(cen_reads)) + "\n")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
