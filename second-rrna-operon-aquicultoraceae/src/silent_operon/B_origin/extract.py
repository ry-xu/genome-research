"""Extract 16S and 23S loci (nhmmer vs Rfam HMMs) from the 65 order-level genomes; label div/norm by blastn
identity to the DASXVE01 divergent and normal references; pair each 23S with a 16S upstream on the same strand.
Usage: python -I src/extract.py   (run from B_origin)
Writes out/loci.tsv, out/lsu_all.fna, out/ssu_all.fna
"""
import subprocess, csv, os, sys
from Bio import SeqIO
from Bio.Seq import Seq

KD = "lab/work/frozen-core/results/lsu/bac/cca/kill_div"
GEN = "data/aq65.fna"
TAX = "lab/work/frozen-core/results/lsu/bac/cca/aq/genomes.tsv"

tax = {}
for r in csv.DictReader(open(TAX), delimiter="\t"):
    tax[r["acc"][3:]] = (r["family"][3:], r["genus"][3:])

seqs = SeqIO.index(GEN, "fasta")


def read_tbl(fn, mlen):
    hits = []
    for line in open(fn):
        if line.startswith("#"):
            continue
        f = line.split()
        tgt, hf, ht, af, at, strand = f[0], int(f[4]), int(f[5]), int(f[6]), int(f[7]), f[11]
        lo, hi = min(af, at), max(af, at)
        hits.append(dict(contig=tgt, hmm_from=hf, hmm_to=ht, start=lo, end=hi, strand=strand, frac=(ht - hf + 1) / mlen))
    return hits


def merge(hits, maxgap=1500):
    """merge nhmmer fragments of one gene (same contig+strand, hmm ranges in order, genomic gap < maxgap)"""
    hits = sorted(hits, key=lambda h: (h["contig"], h["strand"], h["start"]))
    out = []
    for h in hits:
        p = out[-1] if out else None
        if p and p["contig"] == h["contig"] and p["strand"] == h["strand"] and h["start"] - p["end"] < maxgap:
            # model order must increase along the transcript
            ok = (h["hmm_from"] > p["hmm_to"] - 50) if h["strand"] == "+" else (p["hmm_from"] > h["hmm_to"] - 50)
            if ok:
                p["end"] = max(p["end"], h["end"]); p["hmm_from"] = min(p["hmm_from"], h["hmm_from"])
                p["hmm_to"] = max(p["hmm_to"], h["hmm_to"]); p["frac"] += h["frac"]; p["nfrag"] += 1
                continue
        h["nfrag"] = 1
        out.append(dict(h))
    return out


def blast_label(fa, refs):
    """identity over greedy non-overlapping (query and subject) HSPs, best bitscore first"""
    out = {}
    for lab, ref in refs.items():
        r = subprocess.run(["blastn", "-query", fa, "-subject", ref, "-outfmt",
                            "6 qseqid length pident bitscore qstart qend sstart send", "-max_hsps", "30", "-evalue", "1e-5"],
                           capture_output=True, text=True, check=True).stdout
        hs = {}
        for line in r.splitlines():
            q, L, p, b, qs, qe, ss, se = line.split("\t")
            hs.setdefault(q, []).append((float(b), int(L), float(p), int(qs), int(qe), min(int(ss), int(se)), max(int(ss), int(se))))
        for q, lst in hs.items():
            lst.sort(reverse=True)
            took = []
            for h in lst:
                if all(h[4] < t[3] or h[3] > t[4] for t in took) and all(h[6] < t[5] or h[5] > t[6] for t in took):
                    took.append(h)
            L = sum(t[1] for t in took)
            out.setdefault(q, {})[lab] = (L, sum(t[1] * t[2] for t in took) / L)
    return out


rows = []
for kind, tbl, mlen, refs in [("23S", "out/lsu_nhmmer.tbl", 2925, {"div": f"{KD}/s_dx_div23.fa", "norm": f"{KD}/s_dx_norm23.fa"}),
                               ("16S", "out/ssu_nhmmer.tbl", 1533, {"div": f"{KD}/s_dx_div16.fa", "norm": f"{KD}/s_dx_norm16.fa"})]:
    hits = merge(read_tbl(tbl, mlen))
    recs = []
    for i, h in enumerate(hits):
        s = seqs[h["contig"]].seq[h["start"] - 1:h["end"]]
        if h["strand"] == "-":
            s = s.reverse_complement()
        acc = h["contig"].split("|")[0]
        h["id"] = f"{kind}_{i:03d}"
        h["acc"] = acc
        h["kind"] = kind
        h["family"], h["genus"] = tax.get(acc, ("?", "?"))
        h["len"] = len(s)
        h["seq"] = str(s)
        recs.append(h)
    fa = f"out/tmp_{kind}.fna"
    with open(fa, "w") as fh:
        for h in recs:
            fh.write(f">{h['id']}\n{h['seq']}\n")
    lab = blast_label(fa, refs)
    for h in recs:
        d = lab.get(h["id"], {})
        dl, dp = d.get("div", (0, 0)); nl, np_ = d.get("norm", (0, 0))
        h["div_pid"], h["norm_pid"] = round(dp, 1), round(np_, 1)
        h["div_aln"], h["norm_aln"] = dl, nl
        # divergent: closer to div ref than norm ref by >= 5 points over >= 200 nt
        md, mn = dl * dp / 100, nl * np_ / 100
        h["label"] = "div" if (md > mn and md >= 600) else "norm"
        rows.append(h)
    os.remove(fa)

# pair 16S -> 23S into operons (same contig, same strand, 16S upstream, gap < 2 kb)
for h in rows:
    h["operon"] = ""
by_contig = {}
for h in rows:
    by_contig.setdefault(h["contig"], []).append(h)
op = 0
for c, hs in by_contig.items():
    for l in [x for x in hs if x["kind"] == "23S"]:
        for s in [x for x in hs if x["kind"] == "16S" and x["strand"] == l["strand"]]:
            gap = (l["start"] - s["end"]) if l["strand"] == "+" else (s["start"] - l["end"])
            if 0 < gap < 2000:
                op += 1
                l["operon"] = s["operon"] = f"op{op:02d}"
                l["its"] = s["its"] = gap

cols = ["id", "kind", "acc", "family", "genus", "contig", "start", "end", "strand", "hmm_from", "hmm_to", "frac", "nfrag", "len",
        "label", "div_pid", "div_aln", "norm_pid", "norm_aln", "operon"]
with open("out/loci.tsv", "w") as fh:
    fh.write("\t".join(cols) + "\n")
    for h in rows:
        fh.write("\t".join(str(round(h[c], 3) if isinstance(h[c], float) else h[c]) for c in cols) + "\n")
for kind, fn in [("23S", "out/lsu_all.fna"), ("16S", "out/ssu_all.fna")]:
    with open(fn, "w") as fh:
        for h in rows:
            if h["kind"] == kind:
                fh.write(f">{h['id']}|{h['acc']}|{h['genus']}|{h['label']}\n{h['seq']}\n")
print(len(rows), "loci;", sum(1 for h in rows if h["label"] == "div"), "div")
