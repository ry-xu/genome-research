"""RNA test of the CGG skip's mechanism (rule: results/rna_cgg/README.md, fixed 2026-11-12 07:00): is the C of
in-frame CGG present in TANB77 messenger RNA? Run with the atlas Python (proteomics_db imports codon_atlas).
  sites   <mgnify_dir> <genomes.tsv> <pig_check_dir> <out.tsv>
          Set A: C of in-frame CGG met by skip walks of catalogue CDSs in primary genomes; C: the C of in-frame
          CGA/CGC/CGT on the same walks; A0: CGG in non-primary genomes; B: A sites of MGYG000302328 within 3 nt of a
          clean +1 break whose old codon is CGG. Positions are 0-based on the contig's plus strand.
  inject  <sites.tsv> <mgnify_dir> <n_per_set> <seed> <out.fastq>
          150-nt reads from Set A and C sites, half missing one base of the site run, 0.5% substitutions, either strand;
          the truth is in the read name (inj|<site_id>|<keeps|lacks>|k).
  count   <sites.tsv> <out_prefix> <sam...> [--totals-only] [--per-read]
          observations per site (keeps / lacks / other), one per spot; writes <out_prefix>.obs.tsv and prints totals
  screen  <genomes.tsv> <sam...>   reads passing the read filters on primary genomes (prints the count)
  score   <genomes.tsv> <mgnify_dir> <root> <runs.tsv> <out.tsv>
          screen scores per run per M spots: all filtered primary reads, and CDS reads (amendment 2; ranks the table)
  decide  <obs.tsv> <seed>         D per set, bootstrap over sites, and the decision
"""
import bisect
import csv
import os
import random
import re
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from proteomics_db import gff, rc, read_fasta, walk  # noqa: E402
from summarize import flag  # noqa: E402

FLANK, WIN, EDGE, MAPQ, MIN_ID = 20, 5, 30, 20, 0.97
CTRL = {"CGA", "CGC", "CGT"}
CIG = re.compile(r"(\d+)([MIDNSHP=X])")


def run_of(x, p):
    """maximal run of identical bases holding position p: [start, end)"""
    a, b = p, p + 1
    while a > 0 and x[a - 1] == x[p]:
        a -= 1
    while b < len(x) and x[b] == x[p]:
        b += 1
    return a, b


def sites(mdir, gtsv, pdir, out):
    gen = list(csv.DictReader(open(gtsv), delimiter="\t"))
    breaks = []
    # RNA_SETB: the genome whose clean +1 CGG breaks define Set B (chicken replication: GCA_963931445.1)
    fs = Path(pdir) / f"{os.environ.get('RNA_SETB', 'MGYG000302328')}.fs.tsv.gz"
    if fs.exists():
        import gzip
        rows = list(csv.DictReader(gzip.open(fs, "rt"), delimiter="\t"))
        flag(rows)
        breaks = [(r["contig"], int(r["gpos"])) for r in rows if r["clean"] and r["type"] == "+1"
                  and r["old_codon"] == "CGG"]
    # amendment 3: a site is dropped if its observation window (site run + FLANK nt each side) overlaps any ncRNA
    # feature of the GFF (rRNA, tRNA, tmRNA, RNase P...); reads there are not mRNA. walk_has_cgg = 1 if any walk that
    # meets the position meets a CGG (order-independent).
    nc = features(mdir)["ncRNA"]
    rows, n, dropped = [], 0, defaultdict(int)
    for g in gen:
        name, primary = g["genome"], g["primary"] == "1"
        seqs = read_fasta(Path(mdir) / f"{name}.fna")
        site = {}
        for c, a, b, st, cid in gff(Path(mdir) / f"{name}.gff"):
            x = seqs[c]
            s = x[a - 1:] if st == "+" else rc(x[:b])
            _, _, ss, _, pos = walk(s, "skip")
            hits = [(i, "CGG") for _, i in ss] + [(i, s[i:i + 3]) for i in pos if s[i:i + 3] in CTRL]
            for i, codon in hits:
                gp = a - 1 + i if st == "+" else b - 1 - i
                if gp < EDGE or gp >= len(x) - EDGE:
                    continue
                if (c, gp) in site:
                    site[(c, gp)][-1] |= bool(ss)
                else:
                    site[(c, gp)] = [c, gp, st, codon, cid, bool(ss)]
        for c, gp, st, codon, cid, has in site.values():
            x = seqs[c]
            r0, r1 = run_of(x, gp)
            if overlaps(nc.get(c, []), r0 - FLANK, r1 + FLANK):
                dropped[(codon == "CGG", primary)] += 1
                continue
            if codon == "CGG":
                sets = ["A"] if primary else ["A0"]
                if primary and any(bc == c and abs(bp - gp) <= 3 for bc, bp in breaks):
                    sets.append("B")
            else:
                sets = ["C"] if primary else []
            for st_ in sets:
                n += 1
                rows.append([f"s{n}", st_, name, c, gp, st, codon, r0, r1, cid, int(has)])
    with open(out, "w") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["site_id", "set", "genome", "contig", "cpos", "strand", "codon", "run_start", "run_end", "cds",
                    "walk_has_cgg"])
        w.writerows(rows)
    cnt = defaultdict(int)
    for r in rows:
        cnt[r[1]] += 1
    print(dict(cnt), "clean +1 CGG breaks used for B:", len(breaks), "dropped for ncRNA (is CGG, primary):",
          dict(dropped))


def load_sites(path, sets=None):
    by = defaultdict(list)
    for r in csv.DictReader(open(path), delimiter="\t"):
        if sets and r["set"] not in sets:
            continue
        r["run_start"], r["run_end"], r["cpos"] = int(r["run_start"]), int(r["run_end"]), int(r["cpos"])
        by[r["contig"]].append(r)
    for v in by.values():
        v.sort(key=lambda r: r["run_start"])
    return by, {c: [r["run_start"] for r in v] for c, v in by.items()}


def inject(stsv, mdir, nper, seed, out):
    rnd = random.Random(seed)
    by, _ = load_sites(stsv, {"A", "C"})
    allsites = [r for v in by.values() for r in v]
    seqs = {}
    for g in {r["genome"] for r in allsites}:
        seqs.update(read_fasta(Path(mdir) / f"{g}.fna"))
    used = []
    with open(out, "w") as fh:
        for st in ("A", "C"):
            pool = [r for r in allsites if r["set"] == st]
            for r in rnd.sample(pool, min(nper, len(pool))):
                x = seqs[r["contig"]]
                used.append([r["site_id"], st])
                for k, truth in enumerate(("keeps", "lacks")):
                    o = rnd.randint(25, 115)
                    a = r["run_start"] - o
                    seg = x[max(a, 0):a + 151]
                    if a < 0 or len(seg) < 151:
                        continue
                    seg = list(seg)
                    if truth == "lacks":
                        del seg[o]  # one base of the run
                    seg = seg[:150]
                    seg = [rnd.choice([b for b in "ACGT" if b != c]) if rnd.random() < 0.005 else c for c in seg]
                    read = "".join(seg)
                    if rnd.random() < 0.5:
                        read = rc(read)
                    fh.write(f"@inj|{r['site_id']}|{truth}|{k}\n{read}\n+\n{'I' * len(read)}\n")
    with open(out.rsplit(".", 1)[0] + ".sites_used.tsv", "w") as fh:
        fh.write("site_id\tset\n")
        fh.writelines(f"{a}\t{b}\n" for a, b in used)


def aln(f):
    """(contig, ref start, ref end, list of (kind, ref pos, length)) for a read passing the read filters, else None"""
    flagv = int(f[1])
    if flagv & 0x904 or int(f[4]) < MAPQ or f[5] == "*":
        return None
    nm = next((int(t[5:]) for t in f[11:] if t.startswith("NM:i:")), None)
    pos, ev, cols = int(f[3]) - 1, [], 0
    for n, op in CIG.findall(f[5]):
        n = int(n)
        if op in "M=X":
            pos += n
            cols += n
        elif op == "D":
            ev.append(("D", pos, n))
            pos += n
            cols += n
        elif op == "I":
            ev.append(("I", pos, n))
            cols += n
        elif op == "N":
            return None
    if nm is None or cols == 0 or 1 - nm / cols < MIN_ID:
        return None
    return f[2], int(f[3]) - 1, pos, ev


def classify(site, a0, a1, ev):
    r0, r1 = site["run_start"], site["run_end"]
    if a0 > r0 - FLANK or a1 < r1 + FLANK:
        return None
    near = [e for e in ev if r0 - WIN <= e[1] <= r1 + WIN]
    if not near:
        return "keeps"
    if len(near) == 1 and near[0][0] == "D" and near[0][2] == 1 and r0 <= near[0][1] < r1:
        return "lacks"
    return "other"


def count(stsv, prefix, sams, totals_only, per_read):
    by, starts = load_sites(stsv)
    spot = {}  # (sam stem, read name, site_id) -> class
    pr = open(f"{prefix}.per_read.tsv", "w") if per_read else None
    for sam in sams:
        stem = Path(sam).parent.name
        for line in open(sam):
            if line[0] == "@":
                continue
            f = line.split("\t")
            a = aln(f)
            if a is None or a[0] not in by:
                continue
            c, a0, a1, ev = a
            v, st = by[c], starts[c]
            for i in range(max(bisect.bisect_left(st, a0 + FLANK - 200), 0), bisect.bisect_right(st, a1)):
                k = classify(v[i], a0, a1, ev)
                if k is None:
                    continue
                key = (stem, f[0], v[i]["site_id"])
                spot[key] = k if spot.get(key, k) == k else "other"
                if pr:
                    pr.write(f"{f[0]}\t{v[i]['site_id']}\t{v[i]['set']}\t{k}\n")
    obs = defaultdict(lambda: [0, 0, 0])
    for (_, _, sid), k in spot.items():
        obs[sid][("keeps", "lacks", "other").index(k)] += 1
    info = {r["site_id"]: r for v in by.values() for r in v}
    tot = defaultdict(lambda: [0, 0, 0, 0])
    for sid, o in obs.items():
        t = tot[info[sid]["set"]]
        t[0] += o[0] + o[1]
        t[1] += o[2]
        t[2] += (o[0] + o[1]) > 0
    for st, t in sorted(tot.items()):
        print(f"set {st}: observations {t[0]} over {t[2]} sites, other {t[1]}")
    if totals_only:
        return
    with open(f"{prefix}.obs.tsv", "w") as fh:
        fh.write("site_id\tset\tgenome\tcodon\tkeeps\tlacks\tother\n")
        for sid, o in sorted(obs.items(), key=lambda kv: int(kv[0][1:])):
            r = info[sid]
            fh.write(f"{sid}\t{r['set']}\t{r['genome']}\t{r['codon']}\t{o[0]}\t{o[1]}\t{o[2]}\n")


def features(mdir):
    """per contig: sorted (start0, end) of CDS and of ncRNA features (rRNA, tRNA, tmRNA...), from the catalogue GFFs"""
    out = {"CDS": defaultdict(list), "ncRNA": defaultdict(list)}
    for g in Path(mdir).glob("*.gff"):
        for line in open(g):
            if line.startswith("##FASTA"):
                break
            f = line.split("\t")
            if line[0] == "#" or len(f) < 9 or f[2] not in out:
                continue
            out[f[2]][f[0]].append((int(f[3]) - 1, int(f[4])))
    for d in out.values():
        for v in d.values():
            v.sort()
    return out


def overlaps(iv, a0, a1):
    """does [a0, a1) overlap any interval of the sorted list iv (features are < 100 kb)"""
    i = bisect.bisect_left(iv, (a1,))
    return any(iv[j][1] > a0 for j in range(i - 1, max(i - 50, 0) - 1, -1) if iv[j][0] < a1)


def screen(gtsv, sams, feats=None):
    """reads passing the read filters on primary genomes: (all, CDS-only = overlapping a CDS and no ncRNA feature)"""
    prim = {g["genome"] for g in csv.DictReader(open(gtsv), delimiter="\t") if g["primary"] == "1"}
    n = n_cds = 0
    for sam in sams:
        for line in open(sam):
            if line[0] == "@":
                continue
            f = line.split("\t")
            a = aln(f)
            if a is None or a[0].rsplit("_", 1)[0] not in prim:
                continue
            n += 1
            if feats:
                c, a0, a1, _ = a
                n_cds += overlaps(feats["CDS"].get(c, []), a0, a1) and not overlaps(feats["ncRNA"].get(c, []), a0, a1)
    return n, n_cds


def score(gtsv, mdir, root, runs_tsv, out):
    """screen score per run, per million spots streamed: all filtered reads on primary genomes, and the amendment-2
    score (CDS reads: overlapping a CDS and no ncRNA feature), which ranks the table"""
    feats = features(mdir)
    rows = []
    for r in csv.reader(open(runs_tsv), delimiter="\t"):
        run, spots = r[0], int(r[4])
        d = Path(root) / run
        ns = list(d.glob("part_*.n"))
        if not ns:
            rows.append([run, r[1], r[2], r[3], spots, "NA", "NA", "NA", "NA", "NA", "NA"])
            continue
        n_reads = sum(float(x.read_text() or 0) for x in ns)
        streamed = sum(int(m) for x in ns
                       for m in re.findall(r"^Read (\d+) spots", x.with_suffix(".err").read_text(), re.M))
        if streamed != min(1_000_000, spots):
            print(f"{run}: streamed {streamed} spots, ENA read_count {spots}", file=sys.stderr)
        hits, cds = screen(gtsv, [x.with_suffix(".sam") for x in ns], feats)
        rows.append([run, r[1], r[2], r[3], spots, int(n_reads), streamed, hits, f"{hits / streamed * 1e6:.1f}", cds,
                     f"{cds / streamed * 1e6:.1f}"])
    rows.sort(key=lambda x: -1 if x[10] == "NA" else -float(x[10]))
    with open(out, "w") as fh:
        fh.write("run\tstudy\tinstrument\tlayout\tspots\treads_streamed\tspots_streamed\tfiltered_primary\tscore_all"
                 "\tcds_reads\tscore_cds\n")
        fh.writelines("\t".join(map(str, x)) + "\n" for x in rows)


def injeval(fastq, per_read):
    """known-answer test: calls at each injected read's own site, against the rule's pass criteria"""
    made = defaultdict(int)
    sset = {}
    for r in csv.DictReader(open(fastq.rsplit(".", 1)[0] + ".sites_used.tsv"), delimiter="\t"):
        sset[r["site_id"]] = r["set"]
    for i, line in enumerate(open(fastq)):
        if i % 4 == 0:
            _, sid, truth, _ = line[1:].strip().split("|")
            made[(sset[sid], truth)] += 1
    call = defaultdict(lambda: defaultdict(int))
    seen = set()
    for line in open(per_read):
        name, sid, st, k = line.rstrip("\n").split("\t")
        _, tsid, truth, _ = name.split("|")
        if sid != tsid or name in seen:
            continue
        seen.add(name)
        call[(st, truth)][k] += 1
    ok = True
    for st in ("A", "C"):
        rate = {}
        for truth, wrong in (("keeps", "lacks"), ("lacks", "keeps")):
            c = call[(st, truth)]
            n_obs = c["keeps"] + c["lacks"]
            mis = c[wrong] / n_obs if n_obs else 1
            rate[truth] = n_obs / made[(st, truth)]
            ok &= mis <= 0.01
            print(f"set {st} {truth}-reads: made {made[(st, truth)]}, observed {n_obs} ({rate[truth]:.3f}), "
                  f"called {wrong} {c[wrong]} ({mis:.4f}), other {c['other']}")
        ok &= rate["lacks"] >= 0.8 * rate["keeps"]
        print(f"set {st} call-rate ratio lacks/keeps {rate['lacks'] / rate['keeps']:.3f}")
    print("PASS" if ok else "FAIL")


def decide(obs, seed):
    import numpy as np
    rows = list(csv.DictReader(open(obs), delimiter="\t"))
    rng = np.random.default_rng(seed)
    res = {}
    for st in ("A", "C", "A0", "B"):
        k = np.array([int(r["keeps"]) for r in rows if r["set"] == st and int(r["keeps"]) + int(r["lacks"]) > 0])
        l = np.array([int(r["lacks"]) for r in rows if r["set"] == st and int(r["keeps"]) + int(r["lacks"]) > 0])
        if len(k) == 0:
            print(f"set {st}: no observations")
            continue
        d = l.sum() / (k.sum() + l.sum())
        idx = rng.integers(0, len(k), size=(10000, len(k)))
        bd = l[idx].sum(1) / (k[idx].sum(1) + l[idx].sum(1))
        lo, hi = np.quantile(bd, [0.025, 0.975])
        res[st] = (d, lo, hi, int(k.sum() + l.sum()), len(k))
        print(f"set {st}: D {d:.4f} (95% {lo:.4f}-{hi:.4f}), lacks {l.sum()} of {k.sum() + l.sum()} observations "
              f"over {len(k)} sites")
    if "A" not in res or "C" not in res:
        print("decision: not powered")
        return
    d, lo, hi, na, sa = res["A"]
    dc, nc = res["C"][0], res["C"][3]
    if na < 200 or sa < 30 or nc < 1000:
        print("decision: not powered")
    elif hi < 0.10:
        print("decision: transcripts keep the C")
    elif d >= 0.25 and lo >= max(5 * dc, 0.05):
        print("decision: transcripts lack the C (genomic check needed before any claim)")
    else:
        print("decision: inconclusive")


if __name__ == "__main__":
    cmd, a = sys.argv[1], sys.argv[2:]
    if cmd == "sites":
        sites(*a)
    elif cmd == "inject":
        inject(a[0], a[1], int(a[2]), int(a[3]), a[4])
    elif cmd == "count":
        flags = {x for x in a if x.startswith("--")}
        a = [x for x in a if not x.startswith("--")]
        count(a[0], a[1], a[2:], "--totals-only" in flags, "--per-read" in flags)
    elif cmd == "screen":
        print(screen(a[0], a[1:]))  # all filtered reads on primary genomes
    elif cmd == "score":
        score(*a)
    elif cmd == "injeval":
        injeval(*a)
    elif cmd == "decide":
        decide(a[0], int(a[1]))
