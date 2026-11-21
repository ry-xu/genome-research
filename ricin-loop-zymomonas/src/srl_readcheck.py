"""Read-level check of an SRL call in a metagenome-assembled genome. usage:
  python -I src/srl_readcheck.py ACC PAF OUT.tsv
ACC is a census genome; PAF is minimap2 -c --cs=long output of its source-metagenome reads against that genome's census
23S (data/lsu/bac_lsu_all.fna). The SRL dodecamer (E. coli 2654-2665) is located in the census 23S as the unique match
of AGTA[CT]G[ACG][AG]A[AG]GA followed by CC; it must occur once. For every alignment that spans the dodecamer +/- 10 nt
without an indel there, the read's 12 bases are tallied, split by identity outside the site (>= 0.98 = same population;
0.90-0.98 = relatives). Prints and writes the tallies.
"""
import os
import collections
import re
import sys

W = os.environ.get("FROZEN_CORE", "frozen-core")  # lab working tree of the census
PAD = 10


def census_23s(acc):
    seq, take = [], False
    with open(f"{W}/data/lsu/bac_lsu_all.fna") as f:
        for line in f:
            if line.startswith(">"):
                if take:
                    break
                take = line[1:].split()[0] == acc
            elif take:
                seq.append(line.strip())
    return "".join(seq).upper()


def read_bases(cs, tstart, lo, hi):
    """Read bases over target [lo, hi) from a long cs string; None if an indel touches the window."""
    pos, out = tstart, {}
    for op, s in re.findall(r"([=*+\-~:])([A-Za-z0-9]+)", cs):
        if op == "=":
            for i, b in enumerate(s.upper()):
                out[pos + i] = b
            pos += len(s)
        elif op == "*":
            out[pos] = s[1].upper()
            pos += 1
        elif op == "-":
            if pos + len(s) > lo - PAD and pos < hi + PAD:
                return None
            pos += len(s)
        elif op == "+":
            if lo - PAD <= pos <= hi + PAD:
                return None
        elif op == "~":
            return None
        if pos > hi + PAD:
            break
    if all(p in out for p in range(lo - PAD, hi + PAD)):
        return "".join(out[p] for p in range(lo, hi))
    return None


def idbin(x):
    return "≥0.98" if x >= 0.98 else "0.90-0.98" if x >= 0.90 else "<0.90"


def main():
    acc, paf, outp = sys.argv[1:4]
    ref = census_23s(acc)
    hits = [m.start() for m in re.finditer(r"(?=AGTA[CT]G[ACG][AG]A[AG]GACC)", ref)]
    assert len(hits) == 1, f"dodecamer matches in census 23S: {len(hits)}"
    lo = hits[0]; hi = lo + 12
    assembled = ref[lo:hi]
    # v2 (after audit results/lsu/bac/reads/audit/AUDIT.md): identity is computed OUTSIDE the 12 site bases, so a read
    # carrying the other bases is not pushed below 0.98 by them; the alignment must cover >= 90% of the read; reads and
    # fragments (mates .1/.2 of one spot, counted once if they agree) are both reported; alignments touching the site but
    # not spanning it +/- PAD are counted as dropped.
    reads = collections.Counter(); frags = collections.defaultdict(set); places = collections.defaultdict(set)
    n_aln = n_short = n_partial = 0
    # molecule = the set of placements of all alignments of a spot (both mates), after audit2: one-end placements
    # double count a pair's two mates and merge unrelated pairs that share one end
    spot_places = collections.defaultdict(set)
    for line in open(paf):
        x = line.split("\t", 9)
        spot_places[x[0].rsplit(".", 1)[0]].add((int(x[7]), int(x[8]), x[4]))
    for line in open(paf):
        x = line.rstrip("\n").split("\t")
        n_aln += 1
        qlen, qs, qe = int(x[1]), int(x[2]), int(x[3])
        tstart, tend = int(x[7]), int(x[8])
        if tend <= lo or tstart >= hi:
            continue
        if tstart > lo - PAD or tend < hi + PAD:
            n_partial += 1
            continue
        if (qe - qs) < 0.9 * qlen:
            n_short += 1
            continue
        cs = next((t[5:] for t in x[12:] if t.startswith("cs:Z:")), None)
        if cs is None:
            continue
        b = read_bases(cs, tstart, lo, hi)
        if b is None:
            key = ("indel_or_gap", idbin(int(x[9]) / int(x[10])))
        else:
            site_match = sum(c == a for c, a in zip(b, assembled))
            key = (b, idbin((int(x[9]) - site_match) / (int(x[10]) - 12)))
        reads[key] += 1
        frags[key].add(x[0].rsplit(".", 1)[0])
        spot = x[0].rsplit(".", 1)[0]  # fastq-dump --readids: RUN.SPOT.MATE
        places[key].add((spot.split(".")[0], frozenset(spot_places[spot])))  # same run and both-mate placements: PCR/optical duplicates
    with open(outp, "w") as out:
        out.write(f"# {acc} assembled {assembled} at census-23S index {lo}-{hi} (0-based, not E. coli numbering); "
                  f"alignments {n_aln}; touching site but not spanning +/-{PAD}: {n_partial}; aligned <90% of read: {n_short}\n")
        out.write("dodecamer\tidentity_outside_site\treads\tfragments\tmolecules\tmolecules_by_run\n")
        for k, n in sorted(reads.items(), key=lambda kv: (kv[0][1], -kv[1])):
            byrun = collections.Counter(r for r, _ in places[k])
            out.write(f"{k[0]}\t{k[1]}\t{n}\t{len(frags[k])}\t{len(places[k])}\t"
                      + ";".join(f"{r}:{c}" for r, c in sorted(byrun.items())) + "\n")
    print(open(outp).read())


if __name__ == "__main__":
    main()
