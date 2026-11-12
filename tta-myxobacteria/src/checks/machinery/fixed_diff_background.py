"""Is RF2 unusually changed in the TTA clade? Fixed clade-vs-sister differences per gene, genome-wide (post hoc).

Uses the red team's own gene calls and DIAMOND all-vs-all of the 16 panel genomes (src/checks/redteam/genes, allvall.tsv).
Families: genes of a reference sister with a reciprocal best hit in all 15 other genomes. Each ortholog is aligned
to the reference protein (global, BLOSUM62) and mapped to reference coordinates. A site counts if all 16 have a residue
there; it is a fixed difference if all 8 clade genomes share one residue that no sister carries. RF1 and RF2 are found
by blastp of E. coli RF1/RF2 (census markers.faa) against the reference proteins.

Polarity: three Polyangiales outgroups outside Polyangiaceae (Sandaracinus amylolyticus GCF_000737325.1, JACKDX01
GCA_964640175.1, JALVDI01 GCA_965227405.1; Pyrodigal calls, reciprocal best blastp hit to the reference gene). A fixed
difference is clade-derived if the outgroups that align there agree on one residue (at least 2 of them) that a sister
carries; sister-side if that residue is the clade's; otherwise unresolved.

usage: uv run --no-project --with biopython,pyrodigal python -I src/checks/machinery/fixed_diff_background.py [reference]
"""
import collections
import gzip
import subprocess
import sys
import tempfile
from pathlib import Path

import pyrodigal
from Bio import Align, SeqIO
from Bio.Align import substitution_matrices

HERE = Path(__file__).resolve().parent
RT = HERE.parent / "redteam"
MARKERS = HERE.parents[3] / "tag-bacteria" / "data" / "reference" / "markers.faa.gz"
GENOMES = HERE.parents[2] / "data" / "genomes"
OUTGROUPS = ("GCF_000737325.1", "GCA_964640175.1", "GCA_965227405.1")
REF = sys.argv[1] if len(sys.argv) > 1 else "S_035318505"


def read_faa(fn):
    seqs, name = {}, None
    for line in (gzip.open(fn, "rt") if str(fn).endswith(".gz") else open(fn)):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        elif name:
            seqs[name].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def call_genes(acc):
    recs = list(SeqIO.parse(gzip.open(GENOMES / f"{acc}.fna.gz", "rt"), "fasta"))
    finder = pyrodigal.GeneFinder(meta=False)
    finder.train(*(bytes(r.seq) for r in recs))
    out = {}
    for r in recs:
        for i, g in enumerate(finder.find_genes(bytes(r.seq))):
            out[f"O_{acc}|{r.id}_{i}"] = g.translate(translation_table=11).rstrip("*")
    return out


def best_hits(q, s, td):
    """Best blastp hit per query of dict q against dict s -> {query: subject}."""
    qf, sf = Path(td) / "bq.faa", Path(td) / "bs.faa"
    qf.write_text("".join(f">{k}\n{v}\n" for k, v in q.items()))
    sf.write_text("".join(f">{k}\n{v}\n" for k, v in s.items()))
    subprocess.run(["diamond", "makedb", "--in", str(sf), "-d", str(sf), "--quiet"], check=True)
    r = subprocess.run(["diamond", "blastp", "-q", str(qf), "-d", str(sf), "--quiet", "-k", "1", "-e", "1e-5",
                        "--threads", "2", "--outfmt", "6", "qseqid", "sseqid", "bitscore"],
                       capture_output=True, text=True, check=True)
    best = {}
    for line in r.stdout.splitlines():
        a, b, bits = line.split("\t")
        if a not in best or float(bits) > best[a][0]:
            best[a] = (float(bits), b)
    return {a: b for a, (_, b) in best.items()}


def families():
    labels = [l.split("\t")[0] for l in open(RT / "panel.tsv")]
    best = {}  # (query, subject genome) -> (bits, subject)
    with open(RT / "allvall.tsv") as f:
        for line in f:
            a = line.split("\t", 10)
            q, s, bits = a[0], a[1], float(a[9])
            g = s.split("|")[0]
            if g == q.split("|")[0]:
                continue
            k = (q, g)
            if k not in best or bits > best[k][0]:
                best[k] = (bits, s)
    faa = {}
    for l in labels:
        faa.update(read_faa(RT / "genes" / f"{l}.faa"))
    fams = {}
    for q in (k for k in faa if k.startswith(REF + "|")):
        orth = {REF: q}
        for g in labels:
            if g == REF:
                continue
            s = best.get((q, g), (0, None))[1]
            if s and best.get((s, REF), (0, None))[1] == q:
                orth[g] = s
        if len(orth) == len(labels):
            fams[q] = orth
    return labels, faa, fams


def main():
    labels, faa, fams = families()
    clade = [l for l in labels if l.startswith("C_")]
    sisters = [l for l in labels if l.startswith("S_")]
    ref_all = {k: v for k, v in faa.items() if k.startswith(REF + "|")}
    og = {}  # reference gene -> {outgroup: protein}
    with tempfile.TemporaryDirectory() as td:
        for acc in OUTGROUPS:
            prot = call_genes(acc)
            fwd = best_hits({k: ref_all[k] for k in fams}, prot, td)
            rev = best_hits({v: prot[v] for v in set(fwd.values())}, ref_all, td)
            for q, s in fwd.items():
                if rev.get(s) == q:
                    og.setdefault(q, {})[acc] = prot[s]
        m = read_faa(MARKERS)
        qf, sf = Path(td) / "q.faa", Path(td) / "s.faa"
        qf.write_text("".join(f">{k}\n{m[k]}\n" for k in ("RF1|RF1_ECOLI", "RF2|RF2_ECOLI")))
        sf.write_text("".join(f">{k}\n{faa[k]}\n" for k in fams))
        r = subprocess.run(["blastp", "-query", str(qf), "-subject", str(sf), "-evalue", "1e-20", "-outfmt",
                            "6 qseqid sseqid bitscore"], capture_output=True, text=True, check=True)
    tag = {}
    for line in sorted(r.stdout.splitlines(), key=lambda l: -float(l.split("\t")[2])):
        qq, s, _ = line.split("\t")
        if qq not in tag.values() and s not in tag:
            tag[s] = qq
    al = Align.PairwiseAligner()
    al.mode = "global"
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score, al.extend_gap_score = -10, -0.5
    al.end_gap_score = 0

    def columns(ref, seqs):
        cols = collections.defaultdict(dict)
        for g, seq in seqs.items():
            a = al.align(ref, seq)[0]
            for (rs, re_), (qs, qe) in zip(*a.aligned):
                for k in range(re_ - rs):
                    cols[rs + k][g] = seq[qs + k]
        return cols

    rows, rf2_sites = [], []
    for q, orth in fams.items():
        ref = faa[q]
        cols = columns(ref, {g: faa[gene] for g, gene in orth.items()})
        ocols = columns(ref, og.get(q, {}))
        sites = fixed = cd = sd = 0
        for p, d in cols.items():
            if len(d) < len(labels):
                continue
            sites += 1
            c = {d[g] for g in clade}
            sv = {d[g] for g in sisters}
            if len(c) == 1 and not c & sv:
                fixed += 1
                o = collections.Counter(ocols.get(p, {}).values()).most_common()
                side = "unresolved"
                if o and o[0][1] >= 2 and (len(o) == 1 or o[1][1] < o[0][1]):
                    side = "clade" if o[0][0] in sv else "sister" if o[0][0] in c else "unresolved"
                cd += side == "clade"
                sd += side == "sister"
                if tag.get(q) == "RF2|RF2_ECOLI":
                    rf2_sites.append(f"ref {p + 1}: clade {next(iter(c))} sisters {''.join(sorted(sv))} "
                                     f"outgroups {''.join(sorted(ocols.get(p, {}).values()))} -> {side}")
        rows.append((q, len(ref), sites, fixed, cd, sd, len(og.get(q, {}))))
    out = HERE / f"fixed_diff_background_{REF}.tsv"
    with open(out, "w") as f:
        f.write("ref_gene\tlen\tsites\tfixed\tclade_derived\tsister_side\toutgroups\trate\tclade_rate\ttag\n")
        for q, n, s, x, cd, sd, no in sorted(rows, key=lambda r: -(r[4] / r[2] if r[2] else 0)):
            f.write(f"{q}\t{n}\t{s}\t{x}\t{cd}\t{sd}\t{no}\t{x / s if s else 0:.4f}\t{cd / s if s else 0:.4f}\t"
                    f"{tag.get(q, '')}\n")
    ok = [r for r in rows if r[2] >= 100 and r[6] >= 2]
    tot_s = sum(r[2] for r in ok)
    print(f"reference {REF}: {len(fams)} families in all 16; {len(ok)} with >= 100 sites and >= 2 outgroups; "
          f"pooled fixed {sum(r[3] for r in ok)}/{tot_s}, clade-derived {sum(r[4] for r in ok)}, sister-side "
          f"{sum(r[5] for r in ok)}")
    for i, name in ((3, "fixed"), (4, "clade-derived")):
        rates = sorted(r[i] / r[2] for r in ok)
        print(f"  {name}: median rate {rates[len(rates) // 2]:.4f}")
        for r in ok:
            if r[0] in tag:
                rank = 1 + sum(v > r[i] / r[2] for v in rates)
                print(f"  {tag[r[0]]} {name}: {r[i]}/{r[2]} = {r[i] / r[2]:.4f}, rank {rank} of {len(ok)}")
    print("RF2 fixed differences (reference coordinates):")
    print("\n".join("  " + s for s in rf2_sites))


if __name__ == "__main__":
    main()
