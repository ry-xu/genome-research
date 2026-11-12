"""Check 3 of PREREG.md "From an Arm B pass to a result" (operational rules: Changes, 2026-11-11): do the
reads behind a candidate assembly carry X at its Arm B event sites?

Sites: each event of the genome whose window holds X, at X's window offset nearest e*+1 (ties: the earlier).
Control sites: the codon of query residue floor(ql/2) in the same events. Positions come from the gene id
(contig_begin_end_strand, 1-based, as term_codon.call_genes writes it) and are checked against the assembly.
Reads: the WGS runs of data/reads/plan.tsv for the genome (amplicon runs and runs without fastq skipped; Illumina
when a genome has both), streamed from the NCBI SRA mirror (fastq-dump, 4 spot ranges in parallel; ENA refuses the
lab's IP) into minimap2 (-x sr / -x map-ont, --frag=no: each read on its own, --secondary=no, --sam-hit-only),
never stored; only alignments overlapping a site are kept (data/reads/<acc>/<run>.part<i>.sam). Streaming stops
once >= 95% of all sites have >= 30 counted alignments.
Counted per site: primary alignments, MAPQ >= 20, NM <= 3% (Illumina) / 10% (Nanopore) of aligned bases, all three
codon bases aligned (no deletion, no insertion inside the codon), base quality >= 20 (Illumina) / >= 10 (Nanopore)
at all three; codon read on the gene's strand.
Spot ranges use the spot count SRA itself reports (vdb-dump --info SEQ; ENA's read_count can count both mates).
A run read to its end must have fastq-dump report every spot of each range written; any pipeline or reader error
aborts (fixes after the independent audit, scratch/audit_reads/).
Support: depth >= 5 and >= 90% of counted codons are X. Check 3 passes with >= 10 supporting sites.
Depth counts read pairs: one counted alignment per read name per site (fix after audit 2, 2026-11-11 12:25; before,
overlapping mates both counted; results/reads_pairs/ re-tallies the kept SAMs that way, src/reads_retally.py).

usage: uv run --no-project python -I src/reads_check.py <gtdb_acc> [max_runs]
       -> data/reads/<acc>/{sites.tsv,ref.mmi,<run>.sam,<run>.done}, results/reads/<acc>.tsv and .txt
"""
import bisect
import collections
import csv
import fcntl
import gzip
import os
import re
import signal
import threading
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from armB_plan import path  # noqa: E402

G = ROOT / "results" / "gtdb"
MEM = ROOT / "results" / "clade_members"  # GTDB members of the TTA clade outside the Arm B labels (2026-11-11)
OUTD = ROOT / "results" / "reads"
MAPQ, DEPTH, FRAC, NSITES = 20, 5, 0.9, 10
STOP_DEPTH, STOP_SHARE = 30, 0.95
BQ = {"sr": 20, "map-ont": 10}
COMP = str.maketrans("ACGTN", "TGCAN")
CIG = re.compile(r"(\d+)([MIDNSHP=X])")
NMRE = re.compile(r"\bNM:i:(\d+)")
MAXDIV = {"sr": 0.03, "map-ont": 0.10}


def tsv(fn):
    return list(csv.DictReader(open(fn), delimiter="\t"))


def read_fasta(fn):
    op = gzip.open if str(fn).endswith(".gz") else open
    name, buf, out = None, [], {}
    with op(fn, "rt") as f:
        for line in f:
            if line.startswith(">"):
                if name:
                    out[name] = "".join(buf).upper()
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name:
        out[name] = "".join(buf).upper()
    return out


def codon_pos(gene, k):
    """Genome positions (1-based, gene order) of query residue k's codon, contig, strand."""
    contig, b, e, st = gene.rsplit("_", 3)
    b, e, st = int(b), int(e), int(st)
    if st == 1:
        s = b + 3 * (k - 1)
        return contig, st, [s, s + 1, s + 2]
    s = e - 3 * (k - 1)
    return contig, st, [s, s - 1, s - 2]


def genome_codon(seqs, contig, st, pos):
    nt = "".join(seqs[contig][p - 1] for p in pos)
    return nt if st == 1 else nt.translate(COMP)


def make_sites(acc, x, d, mem=None):
    if mem:
        label, evs = mem["label"], [MEM / mem.get("events", "armB_events.tsv")]
    else:
        label = {r["gtdb_acc"]: r["label"] for r in tsv(ROOT / "data" / "gtdb_armB" / "labels.tsv")}[acc]
        evs = [G / f"armB_{i}_events.tsv" for i in (0, 1)]
    if mem and "local_path" in mem:  # catalogue genomes: own path
        fn = Path(mem["local_path"]) if mem["local_path"] else sys.exit(f"{acc}: catalogue row without local_path")
    else:
        fn = path(acc[3:])
    seqs = read_fasta(fn)
    rows = []
    for ev in evs:
        for r in tsv(ev):
            if r["label"] != label:
                continue
            win = r["window"].split(",")
            offs = [j for j, c in enumerate(win) if c == x]
            if not offs:
                continue
            j = min(offs, key=lambda j: (abs(j - 2), j))
            es, ql = int(r["e_star"]), int(r["ql"])
            for kind, k in (("X", es - 1 + j), ("control", ql // 2)):
                contig, st, pos = codon_pos(r["gene"], k)
                c = genome_codon(seqs, contig, st, pos)
                if kind == "X" and c != x:
                    sys.exit(f"site check failed: {r['gene']} residue {k} is {c}, not {x}")
                rows.append(dict(kind=kind, gene=r["gene"], family=r["family"], residue=k, offset=j - 1 if kind == "X"
                                 else "", contig=contig, strand=st, p1=pos[0], p2=pos[1], p3=pos[2], codon=c))
    with open(d / "sites.tsv", "w") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    return fn, rows


def ref_span(pos, cigar):
    return pos + sum(int(n) for n, op in CIG.findall(cigar) if op in "MDN=X") - 1


class Index:
    def __init__(self, sites):
        self.by = collections.defaultdict(list)
        for i, s in enumerate(sites):
            lo, hi = min(int(s["p1"]), int(s["p3"])), max(int(s["p1"]), int(s["p3"]))
            self.by[s["contig"]].append((lo, hi, i))
        for v in self.by.values():
            v.sort()
        self.starts = {c: [t[0] for t in v] for c, v in self.by.items()}

    def hits(self, contig, a, b):
        v = self.by.get(contig)
        if not v:
            return []
        st = self.starts[contig]
        return [i for lo, hi, i in v[bisect.bisect_left(st, a):bisect.bisect_right(st, b)] if hi <= b]


def counted(f, preset):
    """Primary, MAPQ >= 20, NM <= MAXDIV of aligned bases (f: SAM fields split on tabs)."""
    if int(f[1]) & 0x904 or int(f[4]) < MAPQ:
        return False
    m = NMRE.search(f[-1]) or NMRE.search("\t".join(f[11:]))
    al = sum(int(n) for n, op in CIG.findall(f[5]) if op in "M=X")
    return bool(m) and al > 0 and int(m.group(1)) <= MAXDIV[preset] * al


def sra_spots(run):
    """Spot count of an SRA run as the mirror read by fastq-dump reports it (vdb-dump --info, line 'SEQ : n')."""
    out = subprocess.run(["vdb-dump", "--info", run], capture_output=True, text=True).stdout
    m = re.search(r"^SEQ\s*:\s*([\d,]+)\s*$", out, re.M)
    if not m:
        sys.exit(f"{run}: no spot count from vdb-dump --info")
    return int(m.group(1).replace(",", ""))


def stream(run, nspots, preset, mmi, d, idx, cov, streams=4):
    """fastq-dump (NCBI SRA mirror) | minimap2 | keep alignments over sites, in `streams` spot ranges read by
    threads; -> (alignments read, stopped early). Kept lines go to <d>/<run>.part<i>.sam."""
    per = -(-nspots // streams)
    procs, outs = [], []
    for i in range(streams):
        a, b = i * per + 1, min((i + 1) * per, nspots)
        part = d / f"{run}.part{i}"
        cmd = (f"fastq-dump --stdout --split-spot --skip-technical -N {a} -X {b} {run} 2>{part}.fqd.err | "
               f"minimap2 -a -x {preset} --frag=no --secondary=no --sam-hit-only -t 1 {mmi} - 2>{part}.mm2.log")
        procs.append(subprocess.Popen(["/bin/bash", "-o", "pipefail", "-c", cmd], stdout=subprocess.PIPE, text=True,
                                      bufsize=1 << 20, start_new_session=True))
        outs.append(d / f"{run}.part{i}.sam")
    lock, stop, n, err = threading.Lock(), threading.Event(), [0] * streams, []

    def reader(i):
        try:
            with open(outs[i], "w") as out:
                for line in procs[i].stdout:
                    if stop.is_set():
                        break
                    if line[0] == "@":
                        out.write(line)
                        continue
                    n[i] += 1
                    f = line.rstrip("\n").split("\t")
                    if int(f[1]) & 4 or f[5] == "*":
                        continue
                    pos = int(f[3])
                    h = idx.hits(f[2], pos, ref_span(pos, f[5]))
                    if not h:
                        continue
                    out.write(line)
                    if counted(f, preset):
                        with lock:
                            for j in h:
                                cov[j] += 1
        except Exception as e:  # a lost reader must not pass for a finished range
            err.append(f"part{i}: {e!r}")
            stop.set()

    def kill():
        for p in procs:
            if p.poll() is None:
                try:
                    os.killpg(p.pid, signal.SIGTERM)  # our own pipelines' process groups only
                except (ProcessLookupError, PermissionError):
                    pass

    ths = [threading.Thread(target=reader, args=(i,), daemon=True) for i in range(streams)]
    for t in ths:
        t.start()
    t0, early, need = time.time(), False, len(cov)
    while any(t.is_alive() for t in ths):
        time.sleep(30)
        if err:
            kill()
            break
        if not any(t.is_alive() for t in ths):
            break
        done = sum(c >= STOP_DEPTH for c in cov)
        print(f"  {run}: {sum(n)} alignments, {done}/{need} sites at >= {STOP_DEPTH}, {time.time() - t0:.0f} s",
              flush=True)
        if done >= STOP_SHARE * need:
            early = True
            stop.set()
            kill()
            break
    for t in ths:
        t.join()
    rc = [p.wait() for p in procs]
    if err:
        sys.exit(f"{run}: reader error {err}")
    if not early:
        if any(rc):
            sys.exit(f"{run}: pipelines exited {rc} (see {d}/{run}.part*.fqd.err)")
        for i in range(streams):
            a, b = i * per + 1, min((i + 1) * per, nspots)
            w = re.findall(r"Written (\d+) spots", (d / f"{run}.part{i}.fqd.err").read_text())
            if sum(map(int, w)) != b - a + 1:
                sys.exit(f"{run}: part{i} wrote {w} spots, expected {b - a + 1}")
    return sum(n), early


def tally(sams, sites, preset):
    bq = BQ[preset]
    idx = Index(sites)
    counts = [collections.Counter() for _ in sites]
    # one count per read name (read pair) per site: --frag=no maps mates separately, and overlapping mates of short
    # inserts both covered a site (audit 2026-11-11 12:24, scratch/audit_reads_cat/); the first counted mate is kept
    seen = [set() for _ in sites]
    for sam in sams:
        for line in open(sam):
            if line[0] == "@":
                continue
            f = line.rstrip("\n").split("\t")
            contig, pos, cig, seq, qual = f[2], int(f[3]), f[5], f[9], f[10]
            if not counted(f, preset):
                continue
            h = idx.hits(contig, pos, ref_span(pos, cig))
            if not h:
                continue
            # ref position -> (read base, qual) for aligned bases; ref positions followed by an insertion marked
            rp, qp, al, ins_after = pos, 0, {}, set()
            for nn, op in CIG.findall(cig):
                nn = int(nn)
                if op in "M=X":
                    for t in range(nn):
                        al[rp + t] = (seq[qp + t], ord(qual[qp + t]) - 33 if qual != "*" else 99)
                    rp += nn
                    qp += nn
                elif op in "IS":
                    if op == "I":
                        ins_after.add(rp - 1)
                    qp += nn
                elif op in "DN":
                    rp += nn
            for i in h:
                if f[0] in seen[i]:
                    continue
                s = sites[i]
                ps = sorted((int(s["p1"]), int(s["p2"]), int(s["p3"])))
                if any(p not in al for p in ps) or ps[0] in ins_after or ps[1] in ins_after:
                    counts[i]["<indel>"] += 1
                    seen[i].add(f[0])
                    continue
                if any(al[p][1] < bq for p in ps):
                    continue
                nt = "".join(al[p][0] for p in ps)
                counts[i][nt if int(s["strand"]) == 1 else nt.translate(COMP)[::-1]] += 1
                seen[i].add(f[0])
    return counts


def main(acc, max_runs=99):
    cand, mem = [r for r in tsv(G / "candidates.tsv") if r["gtdb_acc"] == acc], None
    if not cand:
        cand = [r for n in ("candidates.tsv", "cat_candidates.tsv") if (MEM / n).exists()
                for r in tsv(MEM / n) if r["gtdb_acc"] == acc]
        mem = cand[0] if cand else None
    assert len(cand) == 1, f"{acc}: {len(cand)} candidate codons; outputs are named by genome only"
    x = cand[0]["codon"]
    d = ROOT / "data" / "reads" / acc
    d.mkdir(parents=True, exist_ok=True)
    lock = open(d / ".lock", "w")  # two runs on one genome truncate each other's SAMs (2026-11-11 09:39)
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        sys.exit(f"{acc}: another reads_check is running on this genome")
    fn, sites = make_sites(acc, x, d, mem)
    nx = sum(s["kind"] == "X" for s in sites)
    print(f"{acc} {x}: {nx} X sites, {len(sites) - nx} control sites (all X sites match the assembly)", flush=True)
    plans = [p for p in (ROOT / "data" / "reads" / n for n in ("plan.tsv", "plan_members.tsv", "plan_cat.tsv")) if p.exists()]
    plan = [r for p in plans for r in tsv(p) if r["gtdb_acc"] == acc and r["fastq_ftp"]
            and r["library_strategy"] != "AMPLICON"]
    if not plan:
        sys.exit(f"{acc}: no WGS run with fastq in plan.tsv")
    short = ("ILLUMINA", "ELEMENT")  # Element AVITI short reads: same preset and thresholds as Illumina (2026-11-11 12:09)
    if any(r["instrument_platform"] in short for r in plan):  # short reads when both exist
        plan = [r for r in plan if r["instrument_platform"] in short]
    plan = list({r["run_accession"]: r for r in plan}.values())[:max_runs]
    plat = {r["instrument_platform"] for r in plan}
    preset = "map-ont" if plat == {"OXFORD_NANOPORE"} else "sr" if plat <= set(short) else sys.exit(f"mixed {plat}")
    mmi = d / "ref.mmi"
    if not mmi.exists():
        subprocess.run(["minimap2", "-x", preset, "-d", str(mmi), str(fn)], check=True, capture_output=True)
    idx = Index(sites)
    cov = [0] * len(sites)
    for r in plan:
        run, done = r["run_accession"], d / f"{r['run_accession']}.done"
        if done.exists():
            print(f"  {run}: already done ({done.read_text().strip()})", flush=True)
            cov = [a + b for a, b in zip(cov, tally_cov(sorted(d.glob(f"{run}.part*.sam")), sites, preset))]
        else:
            t0, nspots = time.time(), sra_spots(run)
            if nspots != int(r["read_count"]):
                print(f"  {run}: SRA {nspots} spots, ENA read_count {r['read_count']}; using SRA", flush=True)
            n, early = stream(run, nspots, preset, mmi, d, idx, cov)
            done.write_text(f"{n} alignments read, early stop {early}, {time.time() - t0:.0f} s, "
                            f"{time.strftime('%Y-%m-%d %H:%M %Z')}\n")
            print(f"  {run}: {done.read_text().strip()}", flush=True)
        if sum(c >= STOP_DEPTH for c in cov) >= STOP_SHARE * len(sites):
            break
    runs = [r["run_accession"] for r in plan if (d / f"{r['run_accession']}.done").exists()]
    sams = [p for run in runs for p in sorted(d.glob(f"{run}.part*.sam"))]
    counts = tally(sams, sites, preset)
    report(acc, x, sites, counts, runs, preset)


def tally_cov(sams, sites, preset):
    """Counted alignments overlapping each site in the kept SAMs (the stop rule's depth)."""
    idx = Index(sites)
    cov = [0] * len(sites)
    for sam in sams:
        for line in open(sam):
            if line[0] == "@":
                continue
            f = line.rstrip("\n").split("\t")
            if not counted(f, preset):
                continue
            for i in idx.hits(f[2], int(f[3]), ref_span(int(f[3]), f[5])):
                cov[i] += 1
    return cov


def report(acc, x, sites, counts, runs, preset):
    OUTD.mkdir(parents=True, exist_ok=True)
    rows = []
    for s, c in zip(sites, counts):
        dep = sum(v for k, v in c.items() if k != "<indel>")
        top = c.most_common(3)
        fx = c[s["codon"]] / dep if dep else 0.0
        rows.append(dict(s, depth=dep, indel=c["<indel>"], frac_assembly_codon=f"{fx:.3f}",
                         supports=int(dep >= DEPTH and fx >= FRAC), top=" ".join(f"{k}:{v}" for k, v in top)))
    with open(OUTD / f"{acc}.tsv", "w") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    lines = [f"{acc} {x}; runs {','.join(runs)}; minimap2 -x {preset}; {time.strftime('%Y-%m-%d %H:%M %Z')}"]
    for kind in ("X", "control"):
        rs = [r for r in rows if r["kind"] == kind]
        cov = [r for r in rs if r["depth"] >= DEPTH]
        sup = sum(r["supports"] for r in rs)
        deps = sorted(r["depth"] for r in rs)
        lines.append(f"{kind:8s} sites {len(rs)}; depth >= {DEPTH}: {len(cov)}; >= {FRAC:.0%} assembly codon: {sup}; "
                     f"covered but not supporting: {len(cov) - sup}; median depth {deps[len(deps) // 2] if deps else 0}")
    sx = sum(r["supports"] for r in rows if r["kind"] == "X")
    lines.append(f"check 3: {sx} supporting X sites -> {'PASS' if sx >= NSITES else 'not passed'} (rule >= {NSITES})")
    bad = [r for r in rows if r["kind"] == "X" and r["depth"] >= DEPTH and not r["supports"]]
    lines += [f"  not supporting: {r['gene']} offset {r['offset']} depth {r['depth']} {r['top']}" for r in bad[:20]]
    txt = "\n".join(lines) + "\n"
    (OUTD / f"{acc}.txt").write_text(txt)
    print(txt, end="")


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 99)
