"""Find core mito genes in one panel genome by tblastn of reference proteins; build exon-chained gene models.

usage: find_genes.py refs_dir msa_dir genome.fna outdir
Writes outdir/<acc>.genes.json : {gene: {"strand":..,"codons":[[gpos, codon, ref, qpos], ...], "segments":[...]}}
and outdir/<acc>.hsps.tsv (all HSPs) and outdir/<acc>.models.tsv (summary).
Reference queries are the MASKED translations (contested codons -> X), so gene models do not
depend on how the references' contested codons are read.
"""
import os
import sys
import json
import subprocess
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import read_fasta, revcomp, T4  # noqa: E402

GENES = ["cox1", "cox2", "cox3", "cob", "atp6", "atp8", "atp9",
         "nad1", "nad2", "nad3", "nad4", "nad4L", "nad5", "nad6"]
EVAL = 1e-3
MINBITS = 30.0
CLUSTER_GAP = 10000
MIN_REMAINDER = 8


def run_tblastn(query, genome, out):
    cmd = ["tblastn", "-query", query, "-subject", genome, "-db_gencode", "4", "-seg", "no",
           "-evalue", str(EVAL), "-max_hsps", "200",
           "-outfmt", "6 qseqid sseqid qstart qend sstart send evalue bitscore qseq sseq sframe"]
    with open(out, "w") as fh:
        subprocess.run(cmd, check=True, stdout=fh)


def hsp_codons(h, gseq):
    """List of (gpos0_first_nt_in_transcript_orientation, codon, qpos or None) for subject residues."""
    s, e = h["sstart"], h["send"]
    out = []
    q = h["qstart"]
    k = 0
    for qa, sa in zip(h["qseq"], h["sseq"]):
        if sa == "-":
            q += 1
            continue
        if s < e:
            p = s - 1 + 3 * k
            cod = gseq[p:p + 3]
            gp = p
        else:
            p = s - 3 * k  # 1-based end of codon on plus strand
            cod = revcomp(gseq[p - 3:p])
            gp = p - 1
        aa = T4.get(cod, "X")
        if aa != sa and not (sa == "X") and not (aa == "*" and sa == "*"):
            raise ValueError(f"translation mismatch {h['qseqid']} {cod} {aa} vs {sa}")
        out.append((gp, cod, q if qa != "-" else None))
        if qa != "-":
            q += 1
        k += 1
    return out


def ref_columns(msa_path):
    out = {}
    for n, row in read_fasta(msa_path):
        out[n.split("|")[0]] = [j for j, ch in enumerate(row) if ch != "-"]
    return out


def nt_range(gp, strand):
    return (gp, gp + 3) if strand == "+" else (gp - 2, gp + 1)


def main():
    refs_dir, msa_dir, genome, outdir = sys.argv[1:5]
    os.makedirs(outdir, exist_ok=True)
    name, gseq = read_fasta(genome)[0]
    acc = name.split()[0]
    gseq = gseq.upper()
    q_all = os.path.join(outdir, f"{acc}.queries.faa")
    with open(q_all, "w") as fh:
        for g in GENES:
            for n, s in read_fasta(os.path.join(refs_dir, f"ref_{g}.mask.faa")):
                fh.write(f">{n}\n{s}\n")
    hsp_file = os.path.join(outdir, f"{acc}.hsps.tsv")
    run_tblastn(q_all, genome, hsp_file)
    hsps = []
    for line in open(hsp_file):
        f = line.rstrip("\n").split("\t")
        h = dict(qseqid=f[0], qstart=int(f[2]), qend=int(f[3]), sstart=int(f[4]), send=int(f[5]),
                 evalue=float(f[6]), bits=float(f[7]), qseq=f[8], sseq=f[9])
        h["ref"], h["gene"] = h["qseqid"].split("|")
        h["strand"] = "+" if h["sstart"] < h["send"] else "-"
        if h["bits"] < MINBITS:
            continue
        h["cods"] = hsp_codons(h, gseq)
        hsps.append(h)
    hsps.sort(key=lambda h: -h["bits"])
    occupied = {}  # nt position -> gene
    accepted = defaultdict(list)
    for h in hsps:
        other = 0
        keep = []
        for gp, cod, qp in h["cods"]:
            a, b = nt_range(gp, h["strand"])
            owners = {occupied.get(x) for x in range(a, b)} - {None}
            if owners - {h["gene"]}:
                other += 1
            if owners:
                keep.append(None)
            else:
                keep.append((gp, cod, qp))
        if other > 10:
            continue  # region already explained by another gene
        # keep longest contiguous run of unoccupied codons
        runs, cur = [], []
        for k in keep:
            if k is None:
                if cur:
                    runs.append(cur)
                cur = []
            else:
                cur.append(k)
        if cur:
            runs.append(cur)
        for run in runs:
            if len(run) < MIN_REMAINDER:
                continue
            seg = dict(ref=h["ref"], gene=h["gene"], strand=h["strand"], bits=h["bits"], evalue=h["evalue"],
                       cods=run, lo=min(nt_range(c[0], h["strand"])[0] for c in run),
                       hi=max(nt_range(c[0], h["strand"])[1] for c in run))
            accepted[h["gene"]].append(seg)
            for gp, cod, qp in run:
                a, b = nt_range(gp, h["strand"])
                for x in range(a, b):
                    occupied[x] = h["gene"]
    Lg = len(gseq)
    models = {}
    summ = []
    for g in GENES:
        segs = accepted.get(g, [])
        if not segs:
            summ.append((g, "absent", 0, 0, "", ""))
            continue
        colmap = ref_columns(os.path.join(msa_dir, f"ref_{g}.msa"))
        for s in segs:
            s["cols"] = [colmap[s["ref"]][qp - 1] if qp is not None else None for gp, cod, qp in s["cods"]]
            cc = [c for c in s["cols"] if c is not None]
            s["cmin"], s["cmax"] = (min(cc), max(cc)) if cc else (None, None)
        segs = [s for s in segs if s["cmin"] is not None]
        # cluster by strand and proximity (circular genome: merge clusters touching both ends)
        clusters = []
        for strand in "+-":
            ss = sorted([s for s in segs if s["strand"] == strand], key=lambda s: s["lo"])
            cur, cl = [], []
            for s in ss:
                if cur and s["lo"] - max(c["hi"] for c in cur) > CLUSTER_GAP:
                    cl.append(cur)
                    cur = []
                cur.append(s)
            if cur:
                cl.append(cur)
            for c in cl:
                for s in c:
                    s["ukey"] = s["lo"]
            if len(cl) > 1 and cl[0][0]["lo"] < CLUSTER_GAP and Lg - max(s["hi"] for s in cl[-1]) < CLUSTER_GAP:
                for s in cl[0]:
                    s["ukey"] = s["lo"] + Lg
                    s["wrapped"] = True
                cl = [cl[-1] + cl[0]] + cl[1:-1]
            clusters.extend(cl)
        clusters.sort(key=lambda c: -sum(s["bits"] for s in c))
        best = clusters[0]
        strand = best[0]["strand"]
        best.sort(key=lambda s: s["ukey"], reverse=(strand == "-"))
        # colinearity: max-bits chain with increasing column ranges (allow small overlap)
        n = len(best)
        dp = [s["bits"] for s in best]
        prev = [-1] * n
        for i in range(n):
            for j in range(i):
                if best[j]["cmin"] < best[i]["cmin"] and best[j]["cmax"] < best[i]["cmax"] and dp[j] + best[i]["bits"] > dp[i]:
                    dp[i] = dp[j] + best[i]["bits"]
                    prev[i] = j
        i = max(range(n), key=lambda k: dp[k])
        chain = []
        while i >= 0:
            chain.append(best[i])
            i = prev[i]
        chain.reverse()
        dropped_noncolinear = [s for s in best if s not in chain]
        # trim overlaps in column space between consecutive segments (drop both sides)
        for a, b in zip(chain, chain[1:]):
            if a["cmax"] >= b["cmin"]:
                lo_b, hi_a = b["cmin"], a["cmax"]
                # tail of a: from first codon with col >= lo_b
                ka = next((k for k, c in enumerate(a["cols"]) if c is not None and c >= lo_b), len(a["cols"]))
                # head of b: through last codon with col <= hi_a
                kb = max((k for k, c in enumerate(b["cols"]) if c is not None and c <= hi_a), default=-1)
                a["trim_tail"], b["trim_head"] = ka, kb + 1
        cods = []
        segrep = []
        for s in chain:
            k0 = s.get("trim_head", 0)
            k1 = s.get("trim_tail", len(s["cods"]))
            run = s["cods"][k0:k1]  # already in transcript order (walked along the HSP)
            cods.extend([[gp, cod, s["ref"], qp] for gp, cod, qp in run])
            segrep.append(dict(ref=s["ref"], lo=s["lo"], hi=s["hi"], bits=s["bits"], n=len(run),
                               cols=[s["cmin"], s["cmax"]], trimmed=len(s["cods"]) - len(run)))
        models[g] = dict(strand=strand, codons=cods, segments=segrep,
                         dropped_noncolinear=[dict(ref=s["ref"], lo=s["lo"], hi=s["hi"], bits=s["bits"]) for s in dropped_noncolinear],
                         dropped_clusters=[[dict(ref=s["ref"], lo=s["lo"], hi=s["hi"], bits=s["bits"]) for s in c] for c in clusters[1:]])
        summ.append((g, strand, len(cods), len(chain),
                     f"{min(s['lo'] for s in chain)+1}-{max(s['hi'] for s in chain)}",
                     ";".join(f"{s['lo']+1}-{s['hi']}:{s['ref']}:{s['bits']:.0f}:c{s['cmin']}-{s['cmax']}:trim{len(s['cods'])-len(s['cods'][s.get('trim_head',0):s.get('trim_tail',len(s['cods']))])}" for s in chain)
                     + (" | dropped_noncolinear=" + ",".join(f"{s['lo']+1}-{s['hi']}" for s in dropped_noncolinear) if dropped_noncolinear else "")))
    json.dump(models, open(os.path.join(outdir, f"{acc}.genes.json"), "w"))
    with open(os.path.join(outdir, f"{acc}.models.tsv"), "w") as fh:
        fh.write("gene\tstrand\tn_codons\tn_segments\tspan\tsegments(lo-hi:ref:bits)\n")
        for r in summ:
            fh.write("\t".join(map(str, r)) + "\n")
    for r in summ:
        print(acc, *r[:5], sep="\t")


if __name__ == "__main__":
    main()
