"""Close-relative ortholog test (red team).

usage: ortho_end.py <allvall.tsv> <genes_dir> <panel.tsv> <outdir>

Reciprocal best hits (RBH, DIAMOND bitscore) between the 16 panel genomes, own Pyrodigal calls (table 11, TTA = Leu).
For an ordered pair (A = ender, B = runner) of RBH orthologs, both genes complete (no partial flag):
  A's alignment reaches within 2 residues of A's C-terminus and covers >= 60% of A; e* = B residue aligned to A's
  last residue (qe + (slen - send)); B "runs on" when B continues >= 5 residues past e*; B "same end" when B ends
  within 0 residues of e* (its stop codon is at e*+1).
Per runner genome B and reference group R (sister / clade, never B itself), each B gene is used once, with its
highest-bitscore qualifying RBH partner in R.
Writes:
  <outdir>/runon_events.tsv  one row per (B gene, R): codons e*-1..e*+3, pident, partner, A stop codon
  <outdir>/end_codon_summary.tsv per (B, R): N events, run-on events, same-end events, TTA at e*+1, E from downstream
           (e*+4..e*+30 of the same run-on events) and from interior aligned codons, Poisson log10 p, TTA rank
  <outdir>/interior.tsv per (B, R): aligned interior codons (residues qs+10..qe-10 of RBH alignments), TTA / CTA / TTG
           counts there; <outdir>/interior_tta_sites.tsv: every interior aligned TTA (gene, residue, distance to end)
"""
import math
import sys
from collections import defaultdict

STOPS = {"TAA", "TAG", "TGA"}
BASES = "TCAG"
CODONS = [a + b + c for a in BASES for b in BASES for c in BASES]
SENSE = [c for c in CODONS if c not in STOPS]


def pois_sf_log10(k, mu):
    """log10 P(X >= k), X ~ Poisson(mu)."""
    if k <= 0:
        return 0.0
    # sum terms from k upward in log space
    lt = -mu + k * math.log(mu) - math.lgamma(k + 1)
    s, t, i = 0.0, 0.0, k
    terms = [lt]
    while True:
        i += 1
        lt = lt + math.log(mu) - math.log(i)
        terms.append(lt)
        if lt < terms[0] - 40 or i > k + 10000:
            break
    m = max(terms)
    return (m + math.log(sum(math.exp(x - m) for x in terms))) / math.log(10)


def main(hits_fn, gdir, panel_fn, outdir):
    panel = [l.rstrip("\n").split("\t") for l in open(panel_fn)]
    role = {r[0]: ("clade" if r[2] == "clade" else "sister") for r in panel}
    labs = [r[0] for r in panel]
    genes = {}
    for lab in labs:
        with open(f"{gdir}/{lab}.genes.tsv") as f:
            next(f)
            for line in f:
                gid, ctg, b, e, st, pb, pe, nt = line.rstrip("\n").split("\t")
                genes[f"{lab}|{gid}"] = (pb == "0" and pe == "0" and len(nt) % 3 == 0, nt)
    # best hit per (query, target genome)
    best = {}
    rows = {}
    with open(hits_fn) as f:
        for line in f:
            r = line.rstrip("\n").split("\t")
            q, s = r[0], r[1]
            gq, gs = q.split("|")[0], s.split("|")[0]
            if gq == gs:
                continue
            bits = float(r[9])
            k = (q, gs)
            if k not in best or bits > best[k][0]:
                best[k] = (bits, s)
                rows[(q, s)] = r
    rbh = []  # (B gene, A gene, row of B query vs A subject)
    for (q, gs), (bits, s) in best.items():
        gq = q.split("|")[0]
        back = best.get((s, gq))
        if back and back[1] == q and (q, s) in rows:
            rbh.append((q, s, rows[(q, s)]))
    print(f"RBH ordered pairs: {len(rbh)}", file=sys.stderr)
    with open(f"{outdir}/rbh.tsv", "w") as f:
        f.write("q\ts\tpident\tqs\tqe\tss\tse\tql\tsl\tbits\n")
        for q, s, r in rbh:
            f.write(f"{q}\t{s}\t{r[2]}\t{r[4]}\t{r[5]}\t{r[6]}\t{r[7]}\t{r[10]}\t{r[11]}\t{r[9]}\n")

    # per (B gene, group) choose best qualifying partner
    chosen = {}
    interior_best = {}
    for q, s, r in rbh:
        gq, gs = q.split("|")[0], s.split("|")[0]
        grp = role[gs]
        okq, ntq = genes[q]
        oka, nta = genes[s]
        pid, qs, qe, ss, se, bits, ql, sl = float(r[2]), int(r[4]), int(r[5]), int(r[6]), int(r[7]), float(r[9]), \
            int(r[10]), int(r[11])
        # interior: best partner per group regardless of ends
        key = (q, grp)
        if okq and (key not in interior_best or bits > interior_best[key][0]):
            interior_best[key] = (bits, qs, qe, ql, s, pid)
        if not (okq and oka):
            continue
        if se < sl - 2 or (se - ss + 1) < 0.6 * sl:
            continue
        if ntq[-3:] not in STOPS or nta[-3:] not in STOPS:
            continue
        if key not in chosen or bits > chosen[key][0]:
            chosen[key] = (bits, s, pid, qe, sl, se, ql, nta[-3:])

    ev_out = open(f"{outdir}/runon_events.tsv", "w")
    ev_out.write("B\tB_role\tR\tgene\tpartner\tpident\tql\te_star\tclass\tA_stop\tw_m1\tw_0\tw_p1\tw_p2\tw_p3\n")
    summ = defaultdict(lambda: {"n": 0, "runon": 0, "same": 0, "short": 0, "p1": defaultdict(int),
                                 "win": defaultdict(int), "down": defaultdict(int), "ndown": 0, "offs": defaultdict(int),
                                 "same_stop": defaultdict(int)})
    for (q, grp), (bits, s, pid, qe, sl, se, ql, astop) in chosen.items():
        gq = q.split("|")[0]
        nt = genes[q][1]
        cods = [nt[k:k + 3] for k in range(0, len(nt), 3)]  # index i = residue i+1; index ql = stop
        es = qe + (sl - se)
        d = summ[(gq, grp)]
        d["n"] += 1
        if ql - es >= 5:
            cls = "runon"
        elif ql == es:
            cls = "same"
        elif ql < es:
            cls = "short"
        else:
            cls = "near"
        d[cls] = d.get(cls, 0) + 1
        w = [cods[i] if 0 <= i < len(cods) else "NA" for i in range(es - 2, es + 3)]  # e*-1 .. e*+3
        ev_out.write(f"{gq}\t{role[gq]}\t{grp}\t{q.split('|', 1)[1]}\t{s}\t{pid}\t{ql}\t{es}\t{cls}\t{astop}\t"
                     + "\t".join(w) + "\n")
        if cls == "same":
            d["same_stop"][cods[ql]] += 1
        if cls == "runon":
            d["p1"][cods[es]] += 1
            for c in set(w[1:]):
                d["win"][c] += 1
            for off, c in zip((-1, 0, 1, 2, 3), w):
                if c == "TTA":
                    d["offs"][off] += 1
            for i in range(es + 3, min(es + 30, ql)):
                d["down"][cods[i]] += 1
                d["ndown"] += 1
    ev_out.close()

    # interior aligned codons
    inter = defaultdict(lambda: defaultdict(int))
    sites = open(f"{outdir}/interior_tta_sites.tsv", "w")
    sites.write("B\tR\tgene\tresidue\tql\tdist_to_end\tqs\tqe\tpartner\tpident\n")
    for (q, grp), (bits, qs, qe, ql, s, pid) in interior_best.items():
        gq = q.split("|")[0]
        nt = genes[q][1]
        for res in range(qs + 10, qe - 10 + 1):
            c = nt[3 * (res - 1):3 * res]
            inter[(gq, grp)][c] += 1
            inter[(gq, grp)]["_all"] += 1
            if c == "TTA":
                sites.write(f"{gq}\t{grp}\t{q.split('|', 1)[1]}\t{res}\t{ql}\t{ql - res}\t{qs}\t{qe}\t{s}\t{pid}\n")
    sites.close()

    with open(f"{outdir}/interior.tsv", "w") as f:
        f.write("B\tB_role\tR\taligned_codons\tLeu\tTTA\tCTA\tTTG\tCTT\tTTA_per1e5\tCTA_per1e5\tTTA_over_CTA\n")
        for (gq, grp), c in sorted(inter.items()):
            leu = sum(c[x] for x in ("TTA", "TTG", "CTT", "CTC", "CTA", "CTG"))
            n = c["_all"]
            f.write(f"{gq}\t{role[gq]}\t{grp}\t{n}\t{leu}\t{c['TTA']}\t{c['CTA']}\t{c['TTG']}\t{c['CTT']}\t"
                    f"{1e5 * c['TTA'] / n:.2f}\t{1e5 * c['CTA'] / n:.2f}\t"
                    f"{(c['TTA'] / c['CTA']) if c['CTA'] else float('nan'):.4f}\n")

    with open(f"{outdir}/end_codon_summary.tsv", "w") as f:
        f.write("B\tB_role\tR\tn\tsame\trunon\tshort\tTTA_p1\tTTA_win\tTTA_offsets\tf_down\tE_down_p1\t"
                "log10p_down\tf_int\tE_int_p1\tlog10p_int\tTTA_rank_by_ratio_down\ttop_codons_p1\n")
        for (gq, grp), d in sorted(summ.items()):
            n_run = d.get("runon", 0)
            o = d["p1"]["TTA"]
            fd = d["down"]["TTA"] / d["ndown"] if d["ndown"] else 0
            ic = inter[(gq, grp)]
            fi = ic["TTA"] / ic["_all"] if ic["_all"] else 0
            ed, ei = n_run * fd, n_run * fi
            ed_f, ei_f = max(ed, 0.5), max(ei, 0.5)
            # rank of TTA among sense codons by O/E_down (E floored at 0.5)
            ratios = {}
            for c in SENSE:
                fc = d["down"][c] / d["ndown"] if d["ndown"] else 0
                ratios[c] = d["p1"][c] / max(n_run * fc, 0.5)
            rank = 1 + sum(1 for c in SENSE if ratios[c] > ratios["TTA"])
            top = sorted(d["p1"].items(), key=lambda x: -x[1])[:6]
            offs = ",".join(f"{k}:{v}" for k, v in sorted(d["offs"].items()))
            f.write(f"{gq}\t{role[gq]}\t{grp}\t{d['n']}\t{d.get('same', 0)}\t{n_run}\t{d.get('short', 0)}\t{o}\t"
                    f"{d['win']['TTA']}\t{offs}\t{fd:.5f}\t{ed:.3f}\t{pois_sf_log10(o, ed_f):.2f}\t{fi:.6f}\t{ei:.3f}\t"
                    f"{pois_sf_log10(o, ei_f):.2f}\t{rank}\t{' '.join(f'{c}:{v}' for c, v in top)}\n")


if __name__ == "__main__":
    main(*sys.argv[1:5])
