"""Re-check, step 3: alignment-free keep/lack counts for one SRA run, or one spot range of it.
Follow-up S copy (recheck_S/): default output dir is recheck_S/counts/; sites still from recheck/sites_recheck.tsv.
  .venv/bin/python -I count_run.py <run> [--from A] [--to B] [--out <dir>]
        streams fastq-dump --stdout --split-spot --skip-technical [-N A] [-X B] <run>; spots A..B (1-based, inclusive).
        Amendment 4: the main set streams spots 1..SPOTS (SPOTS = field 3 of main_runs.txt), possibly as several
        ranges (units.py). --max-spots N is the same as --from 1 --to N.
  .venv/bin/python -I count_run.py <name> --fastq <path> [--from A --to B] [--out <dir>]
        reads a FASTQ file instead (known-answer test); --from/--to then only label the range and are checked against
        the spot numbers in the read names.
A read (either orientation) containing a site's keep probe exactly is a keeps observation, its lack probe a lacks
observation (Aho-Corasick over keep, lack and their reverse complements). One observation per spot and site: reads
of a spot share a name in --split-spot output and come consecutively; if the spot's reads give both classes for a site
it is discordant and dropped. Nothing but counts is stored.
Spot orientation: "m1plus" if read 1 of the spot reads along the genome's plus strand (read 1 matched a plus-strand
probe, or read 2 matched a reverse-complement probe); used only to look at library strandedness (RNA vs DNA reads).
Output: counts/<stem>.tsv (site_id, set, keeps, lacks, discordant, keeps_m1plus, lacks_m1plus; sites with any count)
and counts/<stem>.summary (run, from, to, spots, reads, first_spot, last_spot, exit, seconds), where <stem> is <run>
for a whole run and <run>.<A>-<B> for a range. first_spot/last_spot are the spot numbers in the first and last read
names (the number after the last "." of the name; -1 if not a number), so aggregate.py can check the range streamed.
The .tsv is written only if the stream exits 0.
"""
import csv
import os
import subprocess
import sys
import time
from collections import defaultdict

import ahocorasick

R = "results/rna_cgg/recheck/"  # sites_recheck.tsv (unchanged)
S = "results/rna_cgg/recheck_S/"  # Follow-up S: counts go to S/counts
COMP = str.maketrans("ACGTN", "TGCAN")


def rc(s):
    return s.translate(COMP)[::-1]


def automaton():
    sites = [r for r in csv.DictReader(open(R + "sites_recheck.tsv"), delimiter="\t")
             if r["status"] == "ok" and r["set"] in ("A", "C")]
    pat = {}
    for i, r in enumerate(sites):
        for cls, p in ((0, r["keep_probe"]), (1, r["lack_probe"])):
            for ori, q in ((0, p), (1, rc(p))):
                v = pat.setdefault(q, (i, cls, ori))
                assert v[:2] == (i, cls), ("probe shared by two sites/classes", r["site_id"], v)
    A = ahocorasick.Automaton()
    for q, v in pat.items():
        A.add_word(q, v)
    A.make_automaton()
    return sites, A


def main():
    run, args = sys.argv[1], dict(zip(sys.argv[2::2], sys.argv[3::2]))
    fq = args.get("--fastq")
    odir = args.get("--out", S + "counts") + "/"
    os.makedirs(odir, exist_ok=True)
    if "--max-spots" in args:
        assert "--from" not in args and "--to" not in args
        args["--from"], args["--to"] = "1", args["--max-spots"]
    lo, hi = int(args.get("--from", 1)), (int(args["--to"]) if "--to" in args else None)
    assert lo >= 1 and (hi is None or hi >= lo), (lo, hi)
    ranged = "--from" in args or "--to" in args
    assert not ranged or hi is not None, "a range needs --to"
    stem = f"{run}.{lo}-{hi}" if ranged else run
    out = odir + f"{stem}.tsv"
    sites, A = automaton()
    if fq:
        proc, fh = None, open(fq)
    else:
        err = open(odir + f"{stem}.err", "w")
        lim = ["-N", str(lo), "-X", str(hi)] if ranged else []
        proc = subprocess.Popen(["fastq-dump", "--stdout", "--split-spot", "--skip-technical"] + lim + [run],
                                stdout=subprocess.PIPE, stderr=err, text=True, bufsize=1 << 20)
        fh = proc.stdout
    keeps, lacks, disc = defaultdict(int), defaultdict(int), defaultdict(int)
    kplus, lplus = defaultdict(int), defaultdict(int)
    spots = reads = 0
    first_id = last_id = None
    cur, hits, k = None, {}, 0
    t0 = time.time()

    def flush():
        for i, (c, o) in hits.items():
            if c == 0:
                keeps[i] += 1
                kplus[i] += o == 0
            elif c == 1:
                lacks[i] += 1
                lplus[i] += o == 0
            else:
                disc[i] += 1

    it = iter(fh)
    for h in it:
        seq = next(it)
        next(it)
        next(it)
        assert h[0] == "@", h
        name = h.split(None, 1)[0]
        reads += 1
        if name != cur:
            if hits:
                flush()
                hits = {}
            cur, k = name, 0
            spots += 1
            sid = name.rsplit(".", 1)[-1]
            last_id = int(sid) if sid.isdigit() else -1
            if first_id is None:
                first_id = last_id
        else:
            k += 1
        for _, (i, c, o) in A.iter(seq.rstrip()):
            prev = hits.get(i)
            if prev is None:
                hits[i] = (c, o ^ (k & 1))
            elif prev[0] != c:
                hits[i] = (2, prev[1])
    if hits:
        flush()
    rcode = proc.wait() if proc else 0
    tmp = out + ".tmp"
    with open(tmp, "w") as o:
        o.write("site_id\tset\tkeeps\tlacks\tdiscordant\tkeeps_m1plus\tlacks_m1plus\n")
        for i in sorted(set(keeps) | set(lacks) | set(disc)):
            o.write(f'{sites[i]["site_id"]}\t{sites[i]["set"]}\t{keeps[i]}\t{lacks[i]}\t{disc[i]}\t{kplus[i]}\t'
                    f'{lplus[i]}\n')
    with open(odir + f"{stem}.summary", "w") as o:
        o.write(f"run\t{run}\nfrom\t{lo if ranged else ''}\nto\t{hi if ranged else ''}\nspots\t{spots}\nreads\t{reads}\n"
                f"first_spot\t{first_id}\nlast_spot\t{last_id}\nexit\t{rcode}\nseconds\t{time.time() - t0:.0f}\n")
    if rcode == 0:
        os.replace(tmp, out)
    print(stem, "spots", spots, "reads", reads, "ids", first_id, last_id, "exit", rcode, f"{time.time() - t0:.0f}s",
          flush=True)
    sys.exit(rcode)


if __name__ == "__main__":
    main()
