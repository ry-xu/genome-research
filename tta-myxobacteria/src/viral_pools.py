"""V-B: pools of viral genomes for pooled Arm B (PREREG.md, "Phages and viruses").

  clusters  viral RFs (results/viral_rf_hits.tsv) joined when >= 50% identical over >= 50% of both (DIAMOND
            all-vs-all, single linkage) -> results/viral_rf_clusters.tsv; genomes with no MGV readthrough call
            (INPHARED, RefSeq) get one from Pyrodigal coding density (PREREG Changes 21:13)
  pools     one pool per RF cluster and MGV readthrough call (TAG / TGA / none), one genome per unit (the longest;
            unit = MGV vOTU, INPHARED genus or accession, RefSeq organism),
            and the calibration pools, size-matched to the RF pools (same number of genomes, each drawn from the
            genomes nearest in length):
              ord   MGV Complete/High-quality genomes with no readthrough call (positives under forced 4 and 16,
                    nulls under 11, injection hosts)
              tag   MGV genomes called as reading TAG through (negatives under forced 16)
              tga   MGV genomes called as reading TGA through (negatives under forced 4)
            -> data/phage/pools/<name>.fna, data/phage/pools/pools.tsv, data/phage/pools/args_<set>.txt

usage: python -I src/viral_pools.py clusters      (needs pyrodigal)
       python -I src/viral_pools.py pools <n_per_set>
"""
import gzip
import random
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PH = ROOT / "data" / "phage"
OUT = PH / "pools"
HITS = ROOT / "results" / "viral_rf_hits.tsv"
CLU = ROOT / "results" / "viral_rf_clusters.tsv"
MGV_INFO = PH / "mgv" / "mgv_contig_info.tsv"
MGV_FNA = PH / "mgv" / "mgv_selected.fna"
INP_TAB = PH / "inphared" / "14Apr2025_data.tsv.gz"
INP_FNA = PH / "inphared" / "14Apr2025_genomes.fa.gz"
REFSEQ_FNA = PH / "refseq" / "rf_genomes.fna"
SEED = 20261110
XS = [a + b + c for a in "TCAG" for b in "TCAG" for c in "TCAG" if a + b + c not in ("TAA", "TAG", "TGA", "ATG")]


def fasta(fn):
    op = gzip.open if str(fn).endswith(".gz") else open
    head, buf = None, []
    with op(fn, "rt") as f:
        for line in f:
            if line.startswith(">"):
                if head is not None:
                    yield head, "".join(buf)
                head, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if head is not None:
        yield head, "".join(buf)


def mgv_info():
    info = {}
    with open(MGV_INFO) as f:
        next(f)
        for line in f:
            r = line.rstrip("\n").split("\t")
            info[r[0]] = dict(votu=r[1], length=int(r[2]), quality=r[3], readthrough=r[9])
    return info


def inp_info():
    info = {}
    with gzip.open(INP_TAB, "rt") as f:
        next(f)
        for line in f:
            r = line.rstrip("\n").split("\t")
            info[r[0]] = dict(length=int(r[3]), host=r[14], genus=r[16])
    return info


DENSITY_GAIN = 0.10


def density_call(seq):
    """TGA / TAG if Pyrodigal single-mode coding density under table 4 / 16 is >= 10 points above table 11's (the
    larger gain if both), else NULL; also returns the three densities."""
    import pyrodigal
    b = seq.upper().encode()
    dens = {}
    for t in (11, 4, 16):
        f = pyrodigal.GeneFinder(meta=False)
        f.train(b, translation_table=t)
        dens[t] = sum(g.end - g.begin + 1 for g in f.find_genes(b)) / len(b)
    gain = {"TGA": dens[4] - dens[11], "TAG": dens[16] - dens[11]}
    best = max(gain, key=gain.get)
    return (best if gain[best] >= DENSITY_GAIN else "NULL"), dens


def clusters():
    rows = [line.rstrip("\n").split("\t") for line in open(HITS)]
    parent = {r[1]: r[1] for r in rows}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    with tempfile.TemporaryDirectory() as td:
        q, d, h = Path(td) / "rf.faa", Path(td) / "rf", Path(td) / "rf.tsv"
        q.write_text("".join(f">{r[1]}\n{r[5]}\n" for r in rows))
        subprocess.run(["diamond", "makedb", "--in", str(q), "-d", str(d), "--quiet"], check=True)
        subprocess.run(["diamond", "blastp", "-q", str(q), "-d", str(d), "-o", str(h), "--quiet", "--id", "50",
                        "--query-cover", "50", "--subject-cover", "50", "--max-target-seqs", "0", "--outfmt", "6",
                        "qseqid", "sseqid"], check=True)
        for line in open(h):
            a, b = line.split()
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[max(ra, rb)] = min(ra, rb)
    mi, ii = mgv_info(), inp_info()
    roots = sorted({find(r[1]) for r in rows})
    need = {r[0] for r in rows if r[0] not in mi}
    calls = {}
    for src in (INP_FNA, REFSEQ_FNA):
        if src.exists():
            for h, s in fasta(src):
                if h in need and h not in calls:
                    calls[h] = density_call(s)
                    print(h, calls[h][0], " ".join(f"{t}:{d:.3f}" for t, d in calls[h][1].items()), flush=True)
    with open(CLU, "w") as out:
        out.write("cluster\tacc\tgene\tclass\tpident\tsource\tunit\tlength\treadthrough\n")
        for r in sorted(rows, key=lambda r: (roots.index(find(r[1])), r[0])):
            src = "mgv" if r[0] in mi else "inphared" if r[0] in ii else "refseq"
            if src == "mgv":
                unit, ln, rt = mi[r[0]]["votu"], mi[r[0]]["length"], mi[r[0]]["readthrough"]
            elif src == "inphared":
                g = ii[r[0]]["genus"]
                unit, ln, rt = (g if g != "Unclassified" else r[0]), ii[r[0]]["length"], calls.get(r[0], ("NA",))[0]
            else:
                unit, ln, rt = r[0], 0, calls.get(r[0], ("NA",))[0]
            out.write(f"rf{roots.index(find(r[1])) + 1}\t{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{src}\t{unit}\t{ln}\t{rt}\n")
    print(f"{len(rows)} viral RFs in {len(roots)} clusters -> {CLU}")


TABLE_RT = {"TAG": 16, "TGA": 4}  # forced table of an RF pool by its MGV readthrough call, as its negatives


def write_pool(name, members):
    """members = [(genome id, seq)] -> pool FASTA with contigs named <genome>, total bp."""
    fn = OUT / f"{name}.fna"
    with open(fn, "w") as f:
        for gid, s in members:
            f.write(f">{gid}\n{s}\n")
    return fn, sum(len(s) for _, s in members)


def pick(cands, lengths, rng, used):
    """For each target length, a random genome among the 5 unused candidates (one per vOTU) nearest in length."""
    out = []
    for L in lengths:
        free = [c for c in cands if c[0] not in used and c[1] not in used]
        free.sort(key=lambda c: abs(c[2] - L))
        g = rng.choice(free[:5])
        used.update((g[0], g[1]))
        out.append(g[0])
    return out


def pools(n):
    OUT.mkdir(parents=True, exist_ok=True)
    rng = random.Random(SEED)
    mi = mgv_info()
    clu = [line.rstrip("\n").split("\t") for line in open(CLU)][1:]
    want = {r[1] for r in clu} | {g for g, v in mi.items() if v["quality"] in ("Complete", "High-quality")}
    seqs = {}
    for src in (MGV_FNA, INP_FNA, REFSEQ_FNA):
        if src.exists():
            for h, s in fasta(src):
                if h in want:
                    seqs[h] = s.upper()
    rf_votus = {mi[r[1]]["votu"] for r in clu if r[1] in mi}
    rows, shapes, args = [], [], {k: [] for k in ("rf", "pos4", "pos16", "neg16", "neg4", "null", "inj")}
    by = {}
    for r in clu:
        if r[1] not in seqs:
            print(r[0], r[1], "no genome sequence; skipped", file=sys.stderr)
            continue
        c = r[0] + (f"_{r[8]}" if r[8] in TABLE_RT else "")
        units = by.setdefault(c, {})
        if r[6] not in units or len(seqs[r[1]]) > len(seqs[units[r[6]][0]]):
            units[r[6]] = (r[1], r[6], r[8])
    for c, units in by.items():
        mem = sorted(units.values())
        t = TABLE_RT.get(mem[0][2], 11)
        fn, tot = write_pool(c, [(g, seqs[g]) for g, _, _ in mem])
        rows.append((c, "rf", t, len(mem), tot, ",".join(g for g, _, _ in mem)))
        shapes.append([len(seqs[g]) for g, _, _ in mem])
        args["rf"].append(f"{fn}:{t}")
    # calibration genomes: one per vOTU, never an RF vOTU
    def cand(rt):
        seen, out = set(), []
        for g in sorted(seqs):
            v = mi.get(g)
            if v is None or v["votu"] in rf_votus or v["votu"] in seen or v["readthrough"] != rt:
                continue
            if rt == "NULL" and v["quality"] not in ("Complete", "High-quality"):
                continue
            seen.add(v["votu"])
            out.append((g, v["votu"], len(seqs[g])))
        return out
    sets = {"ord": cand("NULL"), "tag": cand("TAG"), "tga": cand("TGA")}
    print({k: len(v) for k, v in sets.items()}, "calibration genomes (one per vOTU)")
    plan = [("ord", "null", 11), ("ord", "pos", None), ("tag", "neg16", 16), ("tga", "neg4", 4), ("ord", "inj", 11)]
    for kind, role, t in plan:
        used = set()
        for k in range(n):
            shape = shapes[k % len(shapes)]
            try:
                gs = pick(sets[kind], shape, rng, used)
            except IndexError:
                print(role, "ran out of genomes at pool", k, file=sys.stderr)
                break
            name = f"{role}_{k + 1:02d}"
            fn, tot = write_pool(name, [(g, seqs[g]) for g in gs])
            if role == "pos":
                args["pos4"].append(f"{fn}:4")
                args["pos16"].append(f"{fn}:16")
                rows.append((name, role, "4,16", len(gs), tot, ",".join(gs)))
            elif role == "inj":
                X, s, q = rng.choice(XS), (0.03, 0.1)[k % 2], (1.0, 0.0)[(k // 2) % 2]
                args["inj"].append(f"inject:{fn}:{X}:{s}:{q}:{rng.randrange(10 ** 6)}:11")
                rows.append((name, f"inj {X} s={s} q={q}", 11, len(gs), tot, ",".join(gs)))
            else:
                args[role].append(f"{fn}:{t}")
                rows.append((name, role, t, len(gs), tot, ",".join(gs)))
    with open(OUT / "pools.tsv", "w") as f:
        f.write("pool\trole\ttable\tgenomes\tbp\tmembers\n")
        for r in rows:
            f.write("\t".join(map(str, r)) + "\n")
    for k, v in args.items():
        (OUT / f"args_{k}.txt").write_text("".join(a + "\n" for a in v))
    print({k: len(v) for k, v in args.items()}, "runs ->", OUT)


if __name__ == "__main__":
    if sys.argv[1] == "clusters":
        clusters()
    else:
        pools(int(sys.argv[2]))
