"""Candidate stages 1-2 for a flagged genome (README Gate), and its motif sites for the raw-read check (stage 3).
usage: python -I src/candidate.py <outdir> <name> <fasta> [gencode=11]  -> report on stdout,
       <outdir>/cand/<name>.sites.tsv (another agent's reads_check.py breaks.tsv: qseqid strand sseqid shift ctx, plus route,
       gpos, protein)
Stage 1, shifted windows: the flagged motif is counted again with every break moved SHIFT nt up- and downstream along
the gene strand (same frame). Route S: among the breaks that pass every prog test except the abandoned-frame stop,
distinct subject proteins whose first frame-A stop from the tie window lies within A_STOP codons and reads the motif.
Route K: among clean breaks (one per subject protein, as in summarize.py), the fraction whose 30-nt window holds the
4-mer. Passes when each shifted fraction is at most half the fraction at the breaks.
Stage 2, mobile elements: share of the motif proteins whose Swiss-Prot name is a transposase, IS / insertion element,
integrase, phage or viral protein, RF2 or DnaX (passes when < 0.5).
"""
import csv
import gzip
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from summarize import A_STOP, flag, summarize  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
SHIFT = 60
MOBILE = re.compile(r"transpos|insertion element|\bIS\d|integrase|phage|prophage|capsid|terminase|portal protein|"
                    r"tail fiber|release factor 2|DNA polymerase III subunit (tau|gamma)|\bdnaX", re.I)
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def rc(s):
    return s.translate(COMP)[::-1]


def read_fasta(path):
    op = gzip.open if str(path).endswith(".gz") else open
    seqs, name, buf = {}, None, []
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf).upper()
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf).upper()
    return seqs


def names():
    out = {}
    for line in gzip.open(ROOT / "data/sprot_names.tsv.gz", "rt"):
        sid, desc, org = line.rstrip("\n").split("\t")
        out[sid] = (desc, org)
    return out


def gene_strand(seqs, r):
    """gene-strand sequence and the break's index p in it (fs_scan's gpos is p+1 on +, len-p on -)"""
    s = seqs[r["contig"]]
    if r["strand"] == "+":
        return s, int(r["gpos"]) - 1
    return rc(s), len(s) - int(r["gpos"])


def s_motif(s, p, k_lo, k_hi, stops):
    """P codon + '-' + first frame-A stop from the tie window's left edge, if within A_STOP codons of its right edge"""
    for m in range(k_lo, k_hi + A_STOP + 1):
        if p + 3 * m < 3:
            continue
        c = s[p + 3 * m:p + 3 * m + 3]
        if len(c) < 3:
            return None
        if c in stops:
            return s[p + 3 * m - 3:p + 3 * m] + "-" + c
    return None


def main():
    outdir, name, fasta = Path(sys.argv[1]), sys.argv[2], sys.argv[3]
    code = sys.argv[4] if len(sys.argv) > 4 else "11"
    stops = {"TAA", "TAG", "TGA"} - ({"TGA"} if code in ("4", "25", "10") else set()) - (
        {"TAA", "TAG"} if code == "6" else set()) - ({"TAG"} if code in ("15", "16") else set())
    summ = summarize(outdir, name)
    fs = list(csv.DictReader(gzip.open(outdir / f"{name}.fs.tsv.gz", "rt"), delimiter="\t"))
    flag(fs)
    seqs = read_fasta(fasta)
    nm = names()
    sites = []
    print(f"{name}\tloci {summ['loci']}\troute_s {summ['route_s']} {summ['top_motif']} {summ['top_genes']} "
          f"share {summ['top_share']}\troute_k {summ['route_k']} {summ['kmer']} n {summ['kmer_n']} lp {summ['kmer_lp']}")
    if summ["route_s"]:
        typ, mot = summ["top_motif"].split(":")
        good = [r for r in fs if r["good"] and r["type"] == typ]
        cnt = {}
        for r in good:
            s, p = gene_strand(seqs, r)
            for d in (0, -SHIFT, SHIFT):
                if 0 <= p + d and p + d + 3 * (int(r["k_hi"]) + A_STOP) + 3 <= len(s):
                    hit = s_motif(s, p + d, int(r["k_lo"]), int(r["k_hi"]), stops) == mot
                    cnt.setdefault(d, {}).setdefault(r["sseqid"], False)
                    cnt[d][r["sseqid"]] |= hit
                    if d == 0 and hit:
                        sites.append(("S", r))
        fr = {d: sum(v.values()) / max(len(v), 1) for d, v in cnt.items()}
        ok = all(fr.get(d, 0) <= fr.get(0, 0) / 2 for d in (-SHIFT, SHIFT))
        print(f"stage1 S {typ}:{mot}\tat breaks {sum(cnt[0].values())}/{len(cnt[0])} = {fr[0]:.3f}\t"
              + "\t".join(f"shift {d:+d}: {sum(cnt[d].values())}/{len(cnt[d])} = {fr[d]:.3f}" for d in (-SHIFT, SHIFT)
                          if d in cnt) + f"\t{'PASS' if ok else 'FAIL'}")
    if summ["route_k"]:
        km = summ["kmer"]
        one = {}
        for r in fs:
            if r["clean"]:
                one.setdefault(r["sseqid"], r)
        cnt = {0: [0, 0], -SHIFT: [0, 0], SHIFT: [0, 0]}
        for r in one.values():
            if len(r["win_seq"]) != 30:
                continue
            s, p = gene_strand(seqs, r)
            w0 = p + 3 * int(r["k_lo"]) - 9
            assert s[w0:w0 + 30] == r["win_seq"], (name, r["contig"], r["gpos"])
            for d in cnt:
                if 0 <= w0 + d and w0 + d + 30 <= len(s):
                    cnt[d][1] += 1
                    cnt[d][0] += km in s[w0 + d:w0 + d + 30]
            if km in r["win_seq"]:
                sites.append(("K", r))
        fr = {d: a / max(b, 1) for d, (a, b) in cnt.items()}
        ok = all(fr[d] <= fr[0] / 2 for d in (-SHIFT, SHIFT))
        print(f"stage1 K {km}\t" + "\t".join(f"shift {d:+d}: {a}/{b} = {fr[d]:.3f}" for d, (a, b) in cnt.items())
              + f"\t{'PASS' if ok else 'FAIL'}")
    for route in sorted({x for x, _ in sites}):
        prot = {r["sseqid"] for x, r in sites if x == route}
        mob = [p for p in prot if MOBILE.search(" ".join(nm.get(p, ("", ""))))]
        print(f"stage2 {route}\tmotif proteins {len(prot)}\tmobile/RF2/dnaX {len(mob)} = {len(mob) / max(len(prot), 1):.2f}"
              f"\t{'PASS' if len(mob) < len(prot) / 2 else 'FAIL'}")
        for p in sorted(prot):
            print(f"  {route}\t{p}\t{' | '.join(nm.get(p, ('?', '?')))}")
    (outdir / "cand").mkdir(exist_ok=True)
    with open(outdir / "cand" / f"{name}.sites.tsv", "w") as fh:
        fh.write("qseqid\tstrand\tsseqid\tshift\tctx\troute\tgpos\tprotein\n")
        for x, r in sites:
            if len(r["ctx"]) == 40:
                fh.write(f"{r['contig']}\t{r['strand']}\t{r['sseqid']}\t{r['type']}\t{r['ctx'][5:35]}\t{x}\t{r['gpos']}\t"
                         f"{nm.get(r['sseqid'], ('?',))[0]}\n")


if __name__ == "__main__":
    main()
