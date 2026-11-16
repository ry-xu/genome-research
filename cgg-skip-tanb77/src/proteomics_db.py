"""TANB77 protein sets for a metaproteomics test of the CGG skip (for another agent's PXD068444 search).
v3, after an independent review (scratch/proteomics_review.md). Every catalogue CDS is walked from its 5' end with
code 11 until a stop, under three readings of an in-frame CGG:
  skip      the C is skipped and GGN read as Gly, then the walk continues in the +1 frame
  cgg_gly   CGG read as Gly in the old frame (what a frame-blind code tool predicts)
  standard  CGG read as Arg (the catalogue's own reading)
Skip proteins sharing a stop codon are merged (longest kept), since the catalogue splits a gene at its skip sites.
Competitors at each skip site: the stop-to-stop ORF in the site's own frame read through the CGG as Gly or Arg
(cgg_gly, standard), and the ordinary ORF in the +1 frame through the GGN (plus1). Control junctions: for each site,
the nearest other in-frame Arg codon (AGA AGG CGA CGC CGT) of the same skip protein, given the same skip of its first
base (no skip is expected there, so their detection rate is the false-positive rate for junction peptides).
Site peptides: every tryptic peptide (after K/R, not before P; 0-2 missed cleavages; 7-35 aa) covering
  skip, control              the junction residue with >= 2 residues on each side
  cgg_gly, standard, plus1   the site residue
A skip or control peptide counts only if (I = L) no frame-keeping protein yields it as a tryptic peptide: no catalogue
CDS or site ORF under cgg_gly or standard, and no stop-to-stop stretch of a six-frame translation (code 11, and also
CGG = Gly for the 82 representatives) of any of the 341 TANB77 genomes in the catalogue; if no local competitor peptide
is within 10 ppm of its mass (with 0-2 13C isotope errors either way; for a control, its skip protein and +1 ORF); and
if it does not end at a contig edge. A control peptide must also not be yielded by any skip protein. A cgg_gly or
standard peptide counts only if no skip protein yields it and it does not end at a contig edge.
usage: python -I src/proteomics_db.py <mgnify_dir> <genomes.txt> <conspecific_dir> <outdir>
inputs: <mgnify_dir>/<g>.fna, <g>.gff (catalogue CDS), <g>.faa (catalogue proteins, to check the standard reading);
        <conspecific_dir>/*.gff.gz (the other TANB77 genomes of the catalogue, with ##FASTA).
CGG reader: a tRNA-Arg with anticodon CCG or TCG by ARAGORN (codon_atlas.trna, codon-hunt's screen settings).
Primary genomes: an ACG tRNA-Arg is found (a fuller tRNA set) and no CCG or TCG.
outputs: tanb77_search.faa (all readings, competitors and controls, deduplicated by sequence: the database to search),
         tanb77_{skip,cgg_gly,standard}.faa (one reading each), sites.tsv, controls.tsv, peptides.tsv (one row per site
         peptide, for joining search results by sequence with I = L), genomes.tsv
"""
import csv
import gzip
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from codon_atlas import trna

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
BASES = "TCAG"
AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODE = {a + b + c: AA[16 * i + 4 * j + k] for i, a in enumerate(BASES) for j, b in enumerate(BASES)
        for k, c in enumerate(BASES)}
GLY_CODE = dict(CODE, CGG="G")
STARTS = {"ATG", "GTG", "TTG"}
CTRL_CODONS = {"AGA", "AGG", "CGA", "CGC", "CGT"}
MIN_PEP, MAX_PEP, MISSED, PPM, WORKERS = 7, 35, 2, 10, 2
MODES = ("skip", "cgg_gly", "standard")
# monoisotopic residue masses, C carbamidomethylated
MASS = {"G": 57.02146, "A": 71.03711, "S": 87.03203, "P": 97.05276, "V": 99.06841, "T": 101.04768, "C": 160.03065,
        "L": 113.08406, "I": 113.08406, "N": 114.04293, "D": 115.02694, "Q": 128.05858, "K": 128.09496,
        "E": 129.04259, "M": 131.04049, "H": 137.05891, "F": 147.06841, "R": 156.10111, "Y": 163.06333,
        "W": 186.07931}
WATER, C13 = 18.010565, 1.003355


def rc(s):
    return s.translate(COMP)[::-1]


def fasta_lines(lines):
    seqs, name, buf = {}, None, []
    for line in lines:
        if line[0] == ">":
            if name:
                seqs[name] = "".join(buf).upper()
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf).upper()
    return seqs


def read_fasta(path):
    with open(path) as fh:
        return fasta_lines(fh)


def gff_fasta(path):
    """sequences after ##FASTA in a gzipped GFF (MGnify all_genomes)"""
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
        return fasta_lines(fh)


def gff(path):
    cds = []
    for line in open(path):
        if line.startswith("##FASTA"):
            break
        if line[0] == "#":
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        attrs = dict(x.split("=", 1) for x in f[8].split(";") if "=" in x)
        if f[2] == "CDS":
            cds.append((f[0], int(f[3]), int(f[4]), f[6], attrs.get("ID", f"{f[0]}_{f[3]}")))
    return cds


def walk(s, mode):
    """translate s from 0 until a stop; returns (protein, nt consumed, sites as (residue, nt of the C), complete,
    nt position of each residue's codon)"""
    p, pos, sites, i = [], [], [], 0
    while i + 3 <= len(s):
        c = s[i:i + 3]
        if c == "CGG":
            sites.append((len(p), i))
            if mode != "standard":
                pos.append(i + 1 if mode == "skip" else i)
                p.append("G")
                i += 4 if mode == "skip" else 3
                continue
        aa = CODE.get(c, "X")
        if aa == "*":
            return "".join(p), i + 3, sites, True, pos
        p.append(aa)
        pos.append(i)
        i += 3
    return "".join(p), i, sites, False, pos


def orf_start(s, j):
    """first codon of the stop-to-stop ORF holding codon j (in j's frame, not before s[0])"""
    q = j - 3
    while q >= 0 and CODE.get(s[q:q + 3]) != "*":
        q -= 3
    return q + 3


def met(p, start_codon):
    """a GTG or TTG start reads as M"""
    return "M" + p[1:] if start_codon and p[:1] in "VL" else p


def cuts(p):
    return [0] + [j + 1 for j in range(len(p) - 1) if p[j] in "KR" and p[j + 1] != "P"] + [len(p)]


def windows(p, lo, hi):
    """tryptic peptides of p (0-MISSED missed cleavages, MIN_PEP-MAX_PEP aa) covering residues lo..hi:
    (peptide, start, end, missed)"""
    if lo < 0 or hi >= len(p):
        return []
    c, out = cuts(p), []
    for a in range(len(c) - 1):
        if c[a] > lo:
            break
        for b in range(a + 1, min(a + 2 + MISSED, len(c))):
            if c[b] - 1 >= hi and MIN_PEP <= c[b] - c[a] <= MAX_PEP:
                out.append((p[c[a]:c[b]], c[a], c[b], b - a - 1))
    return out


def digest(p):
    """tryptic peptides of p (0-MISSED missed cleavages, MIN_PEP-MAX_PEP aa), I written as L"""
    p = p.replace("I", "L")
    c = cuts(p)
    return {p[c[a]:c[b]] for a in range(len(c) - 1) for b in range(a + 1, min(a + 2 + MISSED, len(c)))
            if MIN_PEP <= c[b] - c[a] <= MAX_PEP}


def mass(p):
    return sum(MASS[a] for a in p) + WATER if all(a in MASS for a in p) else None


def near(m, others):
    """within PPM of another peptide's mass, allowing 0-2 13C isotope errors either way"""
    tol = PPM * 1e-6 * m
    return any(abs(m - x - n * C13) <= tol for x in others for n in range(-2, 3))


def reader(seqs):
    """anticodons of the genome's tRNA-Arg genes (ARAGORN)"""
    return {r[5].replace("U", "T") for r in trna.call(seqs.items()) if r[4].startswith("Arg")}


TARGETS = set()


def _init(targets):
    global TARGETS
    TARGETS = targets


def six_frame(job):
    """TARGETS met as tryptic peptides of the stop-to-stop stretches of a genome's six frames"""
    path, gly = job
    seqs = gff_fasta(path) if str(path).endswith(".gz") else read_fasta(path)
    hits = set()
    for x in seqs.values():
        for y in (x, rc(x)):
            for code in (CODE, GLY_CODE) if gly else (CODE,):
                for f in range(3):
                    t = "".join([code.get(y[i:i + 3], "X") for i in range(f, len(y) - 2, 3)])
                    for piece in t.split("*"):
                        if len(piece) >= MIN_PEP:
                            hits |= digest(piece) & TARGETS
    return hits


def peps(p, lo, hi, complete):
    """site peptides with an edge flag (ends where an incomplete walk ends)"""
    return [(x, mc, int(not complete and e == len(p))) for x, _, e, mc in windows(p, lo, hi)]


def genome(g, mdir, fa, put, frame_world, skip_world):
    """one genome's readings, sites and controls; returns (site records, control records, genomes.tsv row)"""
    seqs = read_fasta(mdir / f"{g}.fna")
    cds = gff(mdir / f"{g}.gff")
    cat = {k: v.rstrip("*") for k, v in read_fasta(mdir / f"{g}.faa").items()} if (mdir / f"{g}.faa").exists() else {}
    best, same, n_cmp, later = {}, 0, 0, []
    for c, a, b, st, cid in cds:
        x = seqs[c]
        s = x[a - 1:] if st == "+" else rc(x[:b])
        sc = s[:3] in STARTS
        prot = {m: walk(s, m) for m in MODES}
        if cid in cat:
            n_cmp += 1
            same += prot["standard"][0][1:] == cat[cid][1:]
        for m in ("standard", "cgg_gly"):  # every CDS under both frame-keeping readings
            pm = met(prot[m][0], sc)
            fa[m].write(f">{m.upper()}|{g}|{cid}\n{pm}\n")
            frame_world.append(pm)
            if m == "standard" or prot["standard"][2]:
                later.append((f"{m.upper()}|{g}|{cid}", pm))
        if not prot["standard"][2]:
            skip_world.append(met(prot["standard"][0], sc))
        used, complete = prot["skip"][1], prot["skip"][3]
        stop = (c, st, a - 1 + used if st == "+" else b - used)  # genomic end of the skip walk
        if not complete:
            stop += (cid,)
        if stop not in best or len(prot["skip"][0]) > len(best[stop][3][0]):
            best[stop] = (cid, s, sc, prot["skip"], c, a, b, st)
    sites, controls, n_prot = [], [], 0
    for cid, s, sc, (sk, _, ss, complete, pos), c, a, b, st in best.values():
        sk = met(sk, sc)
        tag = f"{g}|{cid}"
        fa["skip"].write(f">SKIP|{tag} sites={len(ss)} complete={int(complete)}\n{sk}\n")
        put(f"SKIP|{tag}", sk)
        skip_world.append(sk)
        n_prot += bool(ss)
        site_res = {i for i, _ in ss}
        ctrl_used = set()
        for k, (i, j) in enumerate(ss):
            q = orf_start(s, j)
            alt = {m: walk(s[q:], m) for m in ("cgg_gly", "standard")}
            alt = {m: (met(w[0], sc and q == 0), w[3]) for m, w in alt.items()}
            ia = (j - q) // 3
            assert alt["cgg_gly"][0][ia] == "G" and alt["standard"][0][ia] == "R", (g, cid, k)
            q1 = orf_start(s, j + 1)
            w1 = walk(s[q1:], "standard")
            p1, i1 = (met(w1[0], sc and q1 == 0), w1[3]), (j + 1 - q1) // 3
            for name, (p, _) in (("STANDARD_ORF", alt["standard"]), ("CGG_GLY_ORF", alt["cgg_gly"]),
                                 ("PLUS1_ORF", p1)):
                frame_world.append(p)
                later.append((f"{name}|{tag}|{k + 1}", p))
            comp = [mass(x) for p, _ in (alt["cgg_gly"], alt["standard"], p1) for x in digest(p)]
            comp = [m for m in comp if m]
            rec = {"genome": g, "protein": cid, "site": k + 1, "residue": i + 1,
                   "locus": f"{c}:{a + j if st == '+' else b - j}:{st}",
                   "skip": peps(sk, i - 2, i + 2, complete),
                   "cgg_gly": peps(alt["cgg_gly"][0], ia, ia, alt["cgg_gly"][1]),
                   "standard": peps(alt["standard"][0], ia, ia, alt["standard"][1]),
                   "plus1": peps(p1[0], i1, i1, p1[1]) if i1 < len(p1[0]) else []}
            rec["near"] = {x for x, _, _ in rec["skip"] if mass(x) is None or near(mass(x), comp)}
            sites.append(rec)
            # control: the nearest other in-frame Arg codon of this skip protein, skipped the same way
            opts = [r for r in range(2, len(sk) - 2) if r not in site_res and r not in ctrl_used
                    and s[pos[r]:pos[r] + 3] in CTRL_CODONS]
            if not opts:
                continue
            r = min(opts, key=lambda r: (abs(r - i), r))
            ctrl_used.add(r)
            wc = walk(s[pos[r] + 1:], "skip")
            ctrl = sk[:r] + wc[0]
            qc = orf_start(s, pos[r] + 1)
            pc = walk(s[qc:], "standard")[0]
            cm = [mass(x) for p in (sk, pc) for x in digest(p)]
            cm = [m for m in cm if m]
            put(f"CONTROL|{tag}|{k + 1}", ctrl)
            cp = peps(ctrl, r - 2, r + 2, wc[3])
            controls.append({"genome": g, "protein": cid, "site": k + 1, "residue": r + 1,
                             "codon": s[pos[r]:pos[r] + 3], "control": cp,
                             "near": {x for x, _, _ in cp if mass(x) is None or near(mass(x), cm)}})
    for name, p in later:
        put(name, p)
    rd = reader(seqs)
    primary = int("ACG" in rd and not rd & {"CCG", "TCG"})
    row = [g, len(cds), len(best), n_prot, len(sites), ",".join(sorted(rd)), int(bool(rd & {"CCG", "TCG"})), primary,
           f"{same}/{n_cmp}"]
    for x in sites + controls:
        x["primary"] = primary
    print(*row, len(controls), file=sys.stderr, flush=True)
    return sites, controls, row


def main():
    mdir, glist, cdir, out = Path(sys.argv[1]), sys.argv[2], Path(sys.argv[3]), Path(sys.argv[4])
    out.mkdir(parents=True, exist_ok=True)
    fa = {m: open(out / f"tanb77_{m}.faa", "w") for m in MODES}
    search, seen = open(out / "tanb77_search.faa", "w"), set()

    def put(name, p):
        if p not in seen:
            seen.add(p)
            search.write(f">{name}\n{p}\n")

    names = [x.strip() for x in open(glist) if x.strip()]
    frame_world, skip_world, sites, controls, grows = [], [], [], [], []
    for g in names:
        s, c, row = genome(g, mdir, fa, put, frame_world, skip_world)
        sites += s
        controls += c
        grows.append(row)
    for f in list(fa.values()) + [search]:
        f.close()
    # which junction peptides does a frame-keeping reading yield? (catalogue CDS and site ORFs, then six frames)
    junction = {x.replace("I", "L") for r in sites for x, _, _ in r["skip"]} | \
               {x.replace("I", "L") for r in controls for x, _, _ in r["control"]}
    framed = set()
    for p in frame_world:
        framed |= digest(p) & junction
    jobs = [(mdir / f"{g}.fna", True) for g in names] + [(p, False) for p in sorted(cdir.glob("*.gff.gz"))]
    with ProcessPoolExecutor(WORKERS, initializer=_init, initargs=(junction,)) as ex:
        for n, h in enumerate(ex.map(six_frame, jobs, chunksize=4)):
            framed |= h
            if n % 50 == 0:
                print("six-frame", n, len(jobs), len(framed), file=sys.stderr, flush=True)
    frame_kept = {x.replace("I", "L") for r in sites for k in ("cgg_gly", "standard", "plus1") for x, _, _ in r[k]}
    want, skipped = frame_kept | junction, set()
    for p in skip_world:
        skipped |= digest(p) & want
    print("framed", len(framed), "skipped", len(skipped), file=sys.stderr, flush=True)

    def why(x, edge, kind, near_set):
        y = x.replace("I", "L")
        if edge:
            return "contig_edge"
        if kind in ("skip", "control") and y in framed:
            return "frame_keeping_yields"
        if kind in ("cgg_gly", "standard", "control") and y in skipped:
            return "skip_yields"
        if kind in ("skip", "control") and x in near_set:
            return "near_mass"
        return ""

    pw = csv.writer(open(out / "peptides.tsv", "w"), delimiter="\t", lineterminator="\n")
    pw.writerow(["peptide", "reading", "genome", "protein", "site", "missed", "primary_genome", "dup_site", "counts",
                 "why_not"])
    sw = csv.writer(open(out / "sites.tsv", "w"), delimiter="\t", lineterminator="\n")
    sw.writerow(["genome", "protein", "site", "residue", "locus", "primary_genome", "dup_site", "skip_peptides",
                 "cgg_gly_peptides", "standard_peptides", "plus1_peptides", "skip_ok", "cgg_gly_ok", "standard_ok",
                 "matched"])
    seen_locus = set()
    for r in sites:
        dup = int((r["genome"], r["locus"]) in seen_locus)
        seen_locus.add((r["genome"], r["locus"]))
        ok = {}
        for kind in ("skip", "cgg_gly", "standard", "plus1"):
            good = []
            for x, mc, edge in r[kind]:
                w = why(x, edge, kind, r["near"]) if kind != "plus1" else ("contig_edge" if edge else "")
                pw.writerow([x, kind, r["genome"], r["protein"], r["site"], mc, r["primary"], dup, int(not w), w])
                if not w:
                    good.append(x)
            ok[kind] = good
        sw.writerow([r["genome"], r["protein"], r["site"], r["residue"], r["locus"], r["primary"], dup] +
                    [";".join(ok[k]) for k in ("skip", "cgg_gly", "standard", "plus1")] +
                    [int(bool(ok[k])) for k in ("skip", "cgg_gly", "standard")] +
                    [int(all(ok[k] for k in ("skip", "cgg_gly", "standard")))])
    cw = csv.writer(open(out / "controls.tsv", "w"), delimiter="\t", lineterminator="\n")
    cw.writerow(["genome", "protein", "site", "residue", "codon", "primary_genome", "control_peptides", "control_ok"])
    for r in controls:
        good = []
        for x, mc, edge in r["control"]:
            w = why(x, edge, "control", r["near"])
            pw.writerow([x, "control", r["genome"], r["protein"], r["site"], mc, r["primary"], 0, int(not w), w])
            if not w:
                good.append(x)
        cw.writerow([r["genome"], r["protein"], r["site"], r["residue"], r["codon"], r["primary"], ";".join(good),
                     int(bool(good))])
    gw = csv.writer(open(out / "genomes.tsv", "w"), delimiter="\t", lineterminator="\n")
    gw.writerow(["genome", "cds", "skip_proteins", "proteins_with_sites", "sites", "arg_trna_anticodons", "cgg_reader",
                 "primary", "standard_matches_catalogue"])
    gw.writerows(grows)


if __name__ == "__main__":
    main()
