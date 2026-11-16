"""Known-answer test of count_run.py (alignment-free counter), before any real read is counted.
  .venv/bin/python -I known_answer.py
Simulated spots at 1,000 random Set A and 1,000 random Set C sites (status ok), seed 1:
  - paired spots (fragment 150-400 nt, 150-nt mates, mate 2 reverse-complemented; half the fragments from the minus
    strand), keep or lack version (lack: one base of the site run removed), substitutions at 0.005 per base;
  - per site: 2 keep spots, 2 lack spots, 1 single-end keep or lack read, 1 discordant spot (mate 1 keep, mate 2 lack,
    short fragment so both mates cover the site).
Expected observations are computed from coordinates, not strings: a mate observes a site if the whole probe window
lies in the mate and no substitution falls in it; a site whose window holds the removed base of another site's lack
version is not observed. The FASTQ is piped into count_run.py (nothing written but test/ counts), then compared.
Pass: every site's (keeps, lacks, discordant, keeps_m1plus, lacks_m1plus) equals the expected, for all ok sites
touched.
Split test (Amendment 4 ranges): the same spots are then counted as 3 spot ranges (count_run.py --from/--to, into
test/split/), and aggregate.py must join them to exactly the per-site counts of the unsplit run. Two broken copies
(a missing middle range; a range whose summary claims one spot fewer) must be refused by aggregate.py.
"""
import csv
import os
import random
import shutil
import subprocess
import sys
from bisect import bisect_left
from collections import defaultdict

R = "results/rna_cgg/recheck/"
MDIR = "data/proteomics/mgnify_pig/"
COMP = str.maketrans("ACGTN", "TGCAN")
F = 20


def rc(s):
    return s.translate(COMP)[::-1]


def fasta(path):
    seqs, name, buf = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf).upper()
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf).upper()
    return seqs


def main():
    rnd = random.Random(1)
    ok = [r for r in csv.DictReader(open(R + "sites_recheck.tsv"), delimiter="\t")
          if r["status"] == "ok" and r["set"] in ("A", "C")]
    by_contig = defaultdict(list)
    for r in ok:
        by_contig[(r["genome"], r["contig"])].append((int(r["run_start0"]) - F, int(r["run_end0"]) + F, r["site_id"]))
    for v in by_contig.values():
        v.sort()
    starts = {k: [w[0] for w in v] for k, v in by_contig.items()}
    targets = rnd.sample([r for r in ok if r["set"] == "A"], 1000) + rnd.sample([r for r in ok if r["set"] == "C"], 1000)
    genomes = {g: fasta(MDIR + f"{g}.fna") for g in {r["genome"] for r in targets}}
    expected = defaultdict(lambda: [0, 0, 0, 0, 0])
    fq = []
    n = 0

    def mate(x, positions, lo, hi):
        """bases and original positions of fragment slice [lo, hi), with substitutions"""
        pos = positions[lo:hi]
        bases = list("".join(x[p] for p in pos))
        err = set()
        for k in range(len(bases)):
            if rnd.random() < 0.005:
                bases[k] = rnd.choice([b for b in "ACGT" if b != bases[k]])
                err.add(pos[k])
        return "".join(bases), pos, err

    def observed(key, pos, err, deleted, own, rs):
        """sites observed by a mate covering original positions pos (the last base `deleted` of site `own`'s run
        [rs, deleted] removed). Removing any base of a run gives the same string, so a window holding the run's end
        but not its start is intact under the labelling where the run's first base is the one removed."""
        out = {}
        lo, hi = min(pos), max(pos)
        v, st = by_contig[key], starts[key]
        j = bisect_left(st, lo - 200)
        while j < len(v) and v[j][0] <= hi:
            ws, we, sid = v[j]
            j += 1
            w = set(range(ws, we))
            p_use, e_use, cls = pos, err, 0
            if deleted is not None and deleted in w:
                if sid == own:
                    w.discard(deleted)
                    cls = 1
                elif rs in w:
                    continue  # the whole run is in the window: neither probe can match
                else:
                    shift = lambda p: p + 1 if rs <= p < deleted else p  # noqa: E731
                    p_use, e_use = [shift(p) for p in pos], {shift(p) for p in err}
            if w <= set(p_use) and not (w & e_use):
                out[sid] = cls
        return out

    for t in targets:
        key = (t["genome"], t["contig"])
        x = genomes[t["genome"]][t["contig"]]
        rs, re = int(t["run_start0"]), int(t["run_end0"])
        ws, we = rs - F, re + F
        plan = [("keep", "keep", "pe"), ("keep", "keep", "pe"), ("lack", "lack", "pe"), ("lack", "lack", "pe"),
                (rnd.choice(["keep", "lack"]), None, "se"), ("keep", "lack", "disc")]
        for v1, v2, kind in plan:
            flen = rnd.randint(150, 400) if kind == "pe" else (150 if kind == "se" else rnd.randint(150, 175))
            flen = min(flen, len(x) - 2)
            a = max(0, we - flen + (5 if kind == "disc" else -30))
            b = min(ws - (5 if kind == "disc" else -30), len(x) - flen - 1)
            f0 = rnd.randint(a, max(a, b))
            f0 = max(0, min(f0, len(x) - flen - 1))
            minus = rnd.random() < 0.5
            n += 1
            mates = []
            for vi, (ver, end) in enumerate(((v1, 0), (v2, 1))):
                if ver is None:
                    continue
                d = re - 1 if ver == "lack" else None
                positions = [p for p in range(f0, f0 + flen + 1) if p != d][:flen]
                if minus:
                    positions = positions[::-1]
                lo, hi = (0, 150) if end == 0 else (len(positions) - 150, len(positions))
                lo, hi = max(lo, 0), min(hi, len(positions))
                bases, pos, err = mate(x, positions, lo, hi)
                if minus:  # positions are already in reverse order: complement gives the minus strand
                    bases = bases.translate(COMP)
                seq = bases if end == 0 else rc(bases)
                mates.append((seq, observed(key, pos, err, d, t["site_id"], rs)))
            spot = {}
            for _, obs in mates:
                for sid, c in obs.items():
                    spot[sid] = c if spot.get(sid, c) == c else 2
            for sid, c in spot.items():
                expected[sid][c] += 1
                if c < 2 and not minus:  # read 1 runs along the plus strand
                    expected[sid][3 + c] += 1
            for seq, _ in mates:
                fq.append(f"@sim.{n} {n} length={len(seq)}\n{seq}\n+\n{'I' * len(seq)}\n")
    p = subprocess.run([R + ".venv/bin/python", "-I", R + "count_run.py", "sim", "--fastq", "/dev/stdin",
                        "--out", R + "test"], input="".join(fq), text=True, capture_output=True)
    print(p.stdout.strip(), p.stderr.strip()[-500:])
    got = {r["site_id"]: [int(r[k]) for k in ("keeps", "lacks", "discordant", "keeps_m1plus", "lacks_m1plus")]
           for r in csv.DictReader(open(R + "test/sim.tsv"), delimiter="\t")}
    sset = {r["site_id"]: r["set"] for r in ok}
    tot = defaultdict(lambda: [0] * 10)
    bad = []
    for sid in set(got) | set(expected):
        e, g = expected.get(sid, [0] * 5), got.get(sid, [0] * 5)
        for k in range(5):
            tot[sset[sid]][k] += e[k]
            tot[sset[sid]][5 + k] += g[k]
        if e != g:
            bad.append((sid, e, g))
    for s in sorted(tot):
        print(s, "expected keeps/lacks/disc/keeps_m1plus/lacks_m1plus", tot[s][:5], "got", tot[s][5:])
    print("spots", n, "sites touched", len(set(got) | set(expected)), "sites differing", len(bad))
    for b in bad[:20]:
        print("  ", b)
    print("PASS" if not bad else "FAIL")
    if bad:
        sys.exit(1)
    split_ok = split_test(fq, n)
    print("split test", "PASS" if split_ok else "FAIL")
    sys.exit(0 if split_ok else 1)


def agg(cdir):
    p = subprocess.run([R + ".venv/bin/python", "-I", R + "aggregate.py", "--counts", cdir, "--partial"],
                       text=True, capture_output=True)
    return p.returncode, (p.stderr.strip().splitlines() or [""])[-1]


def split_test(fq, n):
    sdir = R + "test/split/"
    os.makedirs(sdir, exist_ok=True)
    cuts = [(1, n // 3), (n // 3 + 1, 2 * n // 3), (2 * n // 3 + 1, n)]
    for a, b in cuts:
        part = "".join(r for r in fq if a <= int(r.split(None, 1)[0].rsplit(".", 1)[1]) <= b)
        p = subprocess.run([R + ".venv/bin/python", "-I", R + "count_run.py", "sim", "--fastq", "/dev/stdin",
                            "--from", str(a), "--to", str(b), "--out", sdir], input=part, text=True,
                           capture_output=True)
        print(" ", p.stdout.strip())
    ok = True
    rc0, msg0 = agg(R + "test")
    rc1, msg1 = agg(sdir)
    whole = open(R + "test/per_site_counts.tsv").read()
    joined = open(sdir + "per_site_counts.tsv").read()
    print("  aggregate unsplit exit", rc0, "| split exit", rc1, "| per-site counts identical:", whole == joined)
    ok &= rc0 == 0 and rc1 == 0 and whole == joined
    # broken copies, built in the scratch directory: aggregate must refuse them
    scratch = "scratch/"
    for label, keep, edit in (("missing middle range", [0, 2], None), ("one spot short", [0, 1, 2], 1)):
        bdir = scratch + "split_broken_" + label.replace(" ", "_") + "/"
        os.makedirs(bdir, exist_ok=True)
        for j in keep:
            a, b = cuts[j]
            for ext in (".tsv", ".summary"):
                shutil.copy(sdir + f"sim.{a}-{b}{ext}", bdir + f"sim.{a}-{b}{ext}")
        if edit is not None:
            a, b = cuts[edit]
            f = bdir + f"sim.{a}-{b}.summary"
            t = open(f).read().replace(f"spots\t{b - a + 1}\n", f"spots\t{b - a}\n", 1)
            open(f, "w").write(t)
        rcb, msgb = agg(bdir)
        print(f"  broken copy ({label}): aggregate exit {rcb}: {msgb[:150]}")
        ok &= rcb != 0
    return ok


if __name__ == "__main__":
    main()
