"""Marker-protein AAI screen for the Polyangiaceae TTA clade (g__JADLFU01/JANWXX01/DASTPP01/JAXLYQ01/CAYNEK01).

Method
  1. pyrodigal (single mode, table 11) gene calls. Reference/panel genomes: trained on the whole genome. Candidates:
     trained on the first 3 Mb of contigs (speed), genes called on all contigs.
  2. Core markers: proteins of the reference GCA_020849655.1 (JADLFU01, CheckM2 97%) >= 100 aa that hit (diamond
     blastp, e<=1e-10, query & subject cover >= 50%) exactly one protein in >= 85% of the 94 DB genomes
     (10 clade references + 84 outgroup genomes, one GTDB R232 species rep per non-clade f__Polyangiaceae genus).
     Marker DB = each marker's best hit in each of the 94 genomes.
  3. Candidate proteins vs marker DB (same filters). For each marker the candidate's best-scoring protein is kept;
     AAI(candidate, g) = mean percent identity of those proteins' hits to genome g, over markers with a hit in g.
     Reported: nearest clade reference and its AAI, nearest outgroup and its AAI, number of markers.

usage:
  uv run --no-project --with pyrodigal python -I src/catalogues/screen.py predict OUT_FAA_DIR FULL|FAST genome.fna.gz...
  uv run --no-project python -I src/catalogues/screen.py build FAA_DIR DB_DIR
  uv run --no-project python -I src/catalogues/screen.py run DB_DIR FAA_DIR OUT.tsv [exclude_genomes_csv]
"""
import gzip
import subprocess
import sys
from collections import defaultdict
from multiprocessing import Pool
from pathlib import Path

DIAMOND = "/opt/homebrew/bin/diamond"
CLADE = {"GCA_020849655.1": "JADLFU01", "GCA_963169165.1": "JADLFU01", "GCA_963169955.1": "JADLFU01",
         "GCA_964605685.1": "JANWXX01", "GCA_025453595.1": "JANWXX01", "GCA_025057345.1": "JANWXX01",
         "GCA_033578855.1": "JANWXX01", "GCA_035762175.1": "DASTPP01", "GCA_046351315.1": "JAXLYQ01",
         "GCA_964396175.1": "CAYNEK01"}
SEED = "GCA_020849655.1"
FILT = ["--evalue", "1e-10", "--query-cover", "50", "--subject-cover", "50"]


def read_fa(p):
    out, cur, name = [], [], None
    op = gzip.open if str(p).endswith(".gz") else open
    with op(p, "rt") as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    out.append((name, "".join(cur)))
                name, cur = line[1:].split()[0], []
            else:
                cur.append(line.strip())
    if name is not None:
        out.append((name, "".join(cur)))
    return out


def gname(p):
    n = Path(p).name
    for s in (".fna.gz", ".fa.gz", ".fasta.gz", ".fna", ".faa", ".fa", ".fasta"):
        if n.endswith(s):
            return n[: -len(s)]
    return n


def predict_one(args):
    import pyrodigal
    p, outdir, mode = args
    g = gname(p)
    out = Path(outdir) / f"{g}.faa"
    if out.exists():
        return g, -1
    try:
        S = [s.upper().encode() for _, s in read_fa(p) if len(s) >= 500]
    except Exception as e:
        print("bad", p, e, flush=True)
        return g, 0
    if sum(map(len, S)) < 100000:
        return g, 0
    gf = pyrodigal.GeneFinder()
    if mode == "FAST":
        tr, tot = [], 0
        for s in sorted(S, key=len, reverse=True):
            tr.append(s)
            tot += len(s)
            if tot > 3_000_000:
                break
        gf.train(*tr)
    else:
        gf.train(*S)
    n = 0
    tmp = out.with_suffix(".part")
    with open(tmp, "w") as o:
        for i, s in enumerate(S):
            for gene in gf.find_genes(s):
                prot = gene.translate(include_stop=False).replace("*", "X")
                if len(prot) >= 60:
                    n += 1
                    o.write(f">{g}|{n}\n{prot}\n")
    tmp.rename(out)
    return g, n


def predict(outdir, mode, files, procs=3):
    Path(outdir).mkdir(parents=True, exist_ok=True)
    with Pool(procs) as pool:
        for g, n in pool.imap_unordered(predict_one, [(f, outdir, mode) for f in files]):
            print("pred", g, n, flush=True)


def diamond(query, db, k, threads=3):
    cmd = [DIAMOND, "blastp", "-q", str(query), "-d", str(db), "-k", str(k), "--threads", str(threads), *FILT,
           "--quiet", "-f", "6", "qseqid", "sseqid", "pident", "bitscore"]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    for line in p.stdout:
        q, s, pid, bs = line.rstrip("\n").split("\t")
        yield q, s, float(pid), float(bs)
    p.wait()


def build(faa_dir, db_dir):
    faa_dir, db_dir = Path(faa_dir), Path(db_dir)
    db_dir.mkdir(parents=True, exist_ok=True)
    allf = db_dir / "all94.faa"
    gs = sorted(gname(f) for f in faa_dir.glob("*.faa"))
    with open(allf, "w") as o:
        for f in sorted(faa_dir.glob("*.faa")):
            o.write(open(f).read())
    subprocess.run([DIAMOND, "makedb", "--in", str(allf), "-d", str(db_dir / "all94"), "--quiet"], check=True)
    seedq = db_dir / "seed.faa"
    with open(seedq, "w") as o:
        for n, s in read_fa(faa_dir / f"{SEED}.faa"):
            if len(s) >= 100:
                o.write(f">{n}\n{s}\n")
    hits = defaultdict(lambda: defaultdict(list))  # seedprot -> genome -> [(bs, sid)]
    for q, s, pid, bs in diamond(seedq, db_dir / "all94", 2000):
        hits[q][s.split("|")[0]].append((bs, s))
    ng = len(gs)
    markers = []
    for q, d in hits.items():
        single = sum(1 for g in gs if len(d.get(g, [])) == 1)
        if single >= 0.85 * ng:
            markers.append(q)
    print("genomes", ng, "seed proteins with hits", len(hits), "markers", len(markers), flush=True)
    # marker DB: best hit per genome; a protein assigned to at most one marker (keep the highest-scoring)
    best = {}
    for m in markers:
        for g, lst in hits[m].items():
            bs, sid = max(lst)
            if sid not in best or bs > best[sid][0]:
                best[sid] = (bs, m)
    seqs = dict(read_fa(allf))
    with open(db_dir / "markers.faa", "w") as o:
        for sid, (bs, m) in best.items():
            o.write(f">{sid}|{m.split('|')[1]}\n{seqs[sid]}\n")
    subprocess.run([DIAMOND, "makedb", "--in", str(db_dir / "markers.faa"), "-d", str(db_dir / "markers"), "--quiet"],
                   check=True)
    open(db_dir / "genomes.txt", "w").write("\n".join(gs) + "\n")
    print("marker DB proteins", len(best))


def run(db_dir, faa_dir, out_tsv, exclude=""):
    db_dir, faa_dir = Path(db_dir), Path(faa_dir)
    excl = set(x for x in exclude.split(",") if x)
    files = sorted(faa_dir.glob("*.faa"))
    q = db_dir / f"query_{Path(out_tsv).stem}.faa"
    with open(q, "w") as o:
        for f in files:
            o.write(open(f).read())
    # per query genome: marker -> (bs, qprot); per qprot: dbgenome -> best pid
    qbest = defaultdict(dict)          # qgenome -> marker -> (bs, qprot)
    phit = defaultdict(dict)           # qprot -> dbgenome -> (bs, pid)
    for qq, s, pid, bs in diamond(q, db_dir / "markers", 500):
        sg, _, m = s.split("|")
        if sg in excl:
            continue
        qg = qq.split("|")[0]
        if sg == qg:
            continue  # never score a genome against itself
        cur = qbest[qg].get(m)
        if cur is None or bs > cur[0]:
            qbest[qg][m] = (bs, qq)
        d = phit[qq]
        if sg not in d or bs > d[sg][0]:
            d[sg] = (bs, pid)
    q.unlink()
    rows = []
    for f in files:
        qg = gname(f)
        ms = qbest.get(qg, {})
        aai = defaultdict(list)
        top_clade = 0
        for m, (bs, qp) in ms.items():
            d = phit[qp]
            for sg, (b2, pid) in d.items():
                aai[sg].append(pid)
            if d:
                top = max(d.items(), key=lambda x: x[1][0])[0]
                top_clade += top in CLADE
        st = {g: (sum(v) / len(v), len(v)) for g, v in aai.items() if len(v) >= max(20, 0.5 * len(ms))}
        cl = sorted(((a, n, g) for g, (a, n) in st.items() if g in CLADE), reverse=True)
        og = sorted(((a, n, g) for g, (a, n) in st.items() if g not in CLADE), reverse=True)
        rows.append(dict(genome=qg, n_markers=len(ms),
                         clade_ref=cl[0][2] if cl else "", clade_ref_genus=CLADE.get(cl[0][2], "") if cl else "",
                         clade_aai=round(cl[0][0], 2) if cl else 0, clade_n=cl[0][1] if cl else 0,
                         out_ref=og[0][2] if og else "", out_aai=round(og[0][0], 2) if og else 0,
                         out_n=og[0][1] if og else 0,
                         top_hit_clade_frac=round(top_clade / len(ms), 3) if ms else 0))
    new = not Path(out_tsv).exists()
    with open(out_tsv, "a") as o:
        if new:
            o.write("\t".join(rows[0]) + "\n" if rows else "")
        for r in rows:
            o.write("\t".join(str(v) for v in r.values()) + "\n")
    return rows


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "predict":
        predict(sys.argv[2], sys.argv[3], sys.argv[4:])
    elif cmd == "build":
        build(sys.argv[2], sys.argv[3])
    elif cmd == "run":
        run(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5] if len(sys.argv) > 5 else "")
