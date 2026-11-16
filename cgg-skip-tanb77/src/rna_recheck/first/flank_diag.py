"""Diagnostic (not part of the rule): is the low observation count low coverage, or exact matching losing reads?
  .venv/bin/python -I flank_diag.py <run> <max_spots>
Streams the first <max_spots> spots (fastq-dump -X) and counts, per Set A / Set C site (status ok), the spots in which
a read holds the site's left 20-nt flank or right 20-nt flank (either orientation), and of those the spots with a
keep or lack probe match. A spot whose reads hold a flank but no probe is a read that ends inside the window, has an
error or variant in it, or has some other indel. Flank 20-mers are not uniqueness-checked (diagnostic only).
Output: diag/<run>.tsv (site_id, set, flank_spots, probe_spots, keeps, lacks) and diag/<run>.summary
"""
import csv
import os
import subprocess
import sys
from collections import defaultdict

import ahocorasick

R = "results/rna_cgg/recheck/"
COMP = str.maketrans("ACGTN", "TGCAN")


def rc(s):
    return s.translate(COMP)[::-1]


def main():
    run, xmax = sys.argv[1], int(sys.argv[2])
    sites = [r for r in csv.DictReader(open(R + "sites_recheck.tsv"), delimiter="\t")
             if r["status"] == "ok" and r["set"] in ("A", "C")]
    A = ahocorasick.Automaton()
    pat = defaultdict(set)
    for i, r in enumerate(sites):
        k, l = r["keep_probe"], r["lack_probe"]
        for kind, p in (("F", k[:20]), ("F", k[-20:]), ("K", k), ("L", l)):
            for q in (p, rc(p)):
                pat[q].add((i, kind))
    for q, v in pat.items():
        A.add_word(q, tuple(v))
    A.make_automaton()
    os.makedirs(R + "diag", exist_ok=True)
    err = open(R + f"diag/{run}.err", "w")
    proc = subprocess.Popen(["fastq-dump", "--stdout", "--split-spot", "--skip-technical", "-X", str(xmax), run],
                            stdout=subprocess.PIPE, stderr=err, text=True, bufsize=1 << 20)
    flank, probe, keeps, lacks = (defaultdict(int) for _ in range(4))
    spots = reads = flank_reads = 0
    cur, hit = None, defaultdict(set)

    def flush():
        for i, kinds in hit.items():
            if "F" in kinds:
                flank[i] += 1
                if kinds & {"K", "L"}:
                    probe[i] += 1
            keeps[i] += "K" in kinds and "L" not in kinds
            lacks[i] += "L" in kinds and "K" not in kinds

    it = iter(proc.stdout)
    for h in it:
        seq = next(it)
        next(it)
        next(it)
        name = h.split(None, 1)[0]
        reads += 1
        if name != cur:
            if hit:
                flush()
                hit = defaultdict(set)
            cur = name
            spots += 1
        any_f = False
        for _, v in A.iter(seq.rstrip()):
            for i, kind in v:
                hit[i].add(kind)
                any_f |= kind == "F"
        flank_reads += any_f
    if hit:
        flush()
    rcode = proc.wait()
    with open(R + f"diag/{run}.tsv", "w") as o:
        o.write("site_id\tset\tflank_spots\tprobe_spots\tkeeps\tlacks\n")
        for i in sorted(set(flank) | set(keeps) | set(lacks)):
            o.write(f'{sites[i]["site_id"]}\t{sites[i]["set"]}\t{flank[i]}\t{probe[i]}\t{keeps[i]}\t{lacks[i]}\n')
    tot = {s: [sum(v[i] for i in v if sites[i]["set"] == s) for v in (flank, probe, keeps, lacks)] for s in "AC"}
    with open(R + f"diag/{run}.summary", "w") as o:
        o.write(f"run\t{run}\nspots\t{spots}\nreads\t{reads}\nreads_with_a_flank\t{flank_reads}\nexit\t{rcode}\n")
        for s in "AC":
            o.write(f"{s}\tflank_spots={tot[s][0]}\tprobe_spots={tot[s][1]}\tkeeps={tot[s][2]}\tlacks={tot[s][3]}\n")
    print(run, spots, reads, flank_reads, tot, "exit", rcode, flush=True)


if __name__ == "__main__":
    main()
