"""Which residue sits at a CGG that is followed by a frameshift, and what is the local nucleotide context?

For every frameshift event (in an alignment that is the best one for its reference at that place), look for a CGG
in the old (upstream) reading frame within 3 codons before / 2 codons from the event. For a +1 event with CGG at
oriented offset p, build the hypothesis translation  old-frame ... | X (= the 4 nt p..p+3) | new frame from p+4 ...
and align it locally (BLOSUM62, gap 11/1) to the reference protein around the event. The reference residue facing X
is the residue the reference has where the genome has the CGG+1nt unit. For -1 events the unit is CGG with the
next codon starting at p+2 (X = 2 nt overlap); reported separately.
Also tabulates, for every event, the codon just before it in the old frame (null comparison).

usage: 04_fs_detail.py LABEL [LABEL...]  -> out/LABEL.fs_cgg_sites.tsv, out/fs_detail_summary.json (appended dict)
"""
import sys, json
from collections import Counter, defaultdict
sys.path.insert(0, "src/recheck")
from mp import read_fasta, parse_paf, gcoord, best_per_ref, STD
from Bio import Align
from Bio.Align import substitution_matrices

W = "recheck_work"
prots = read_fasta(f"{W}/refs/refs_all.faa")
al = Align.PairwiseAligner(mode="local", substitution_matrix=substitution_matrices.load("BLOSUM62"),
                           open_gap_score=-11, extend_gap_score=-1)

def tr(nt):
    return "".join(STD.get(nt[i:i + 3], "X") for i in range(0, len(nt) - 2, 3))

def residue_at(hyp, xi, ref):
    a = al.align(hyp, ref)[0]
    # map hyp index -> ref index via aligned blocks
    blocks_h, blocks_r = a.aligned
    for (h0, h1), (r0, r1) in zip(blocks_h, blocks_r):
        if h0 <= xi < h1:
            ri = r0 + (xi - h0)
            # clean if X is >= 3 residues inside an ungapped block
            clean = (xi - h0 >= 3) and (h1 - xi >= 4)
            lo, hi = max(h0, xi - 10), min(h1, xi + 11)
            idn = sum(1 for k in range(lo, hi) if k != xi and hyp[k] == ref[r0 + k - h0]) / max(1, hi - lo - 1)
            return ref[ri], clean, idn, a.score
    return "-", False, 0.0, a.score

summ = {}
for lab in sys.argv[1:]:
    genome = read_fasta(f"{W}/genomes/{lab}.fna")
    alns = [a for a in parse_paf(f"{W}/out/{lab}.paf", genome, prots) if a["nM"] >= 30 and a["AS"] >= 60]
    best = best_per_ref(alns)
    sites = defaultdict(list)   # (contig, strand, gpos of CGG, dir) -> [(ref, res, clean, idn)]
    ctx = {}; ctxl = {}
    prev_codon = Counter(); n_ev = Counter(); n_ev_cgg = Counter(); seen_ev = set()
    for ai, a in enumerate(alns):
        seg = a["seg"]
        # qpos lookup from M columns: oriented offset -> query position
        mcols = {o: k for k, (o, *_r) in enumerate(a["recs"])}
        for (e, d, n, op) in a["events"]:
            if d == 0: continue
            b = best.get((a["contig"], a["strand"], gcoord(a, e - 3), a["ref"]))
            if b is None or b[1] != ai: continue
            evkey = (a["contig"], a["strand"], gcoord(a, e, 1) // 12, d)
            first = evkey not in seen_ev
            seen_ev.add(evkey)
            cand = [e + 3 * j for j in (-1, 0, -2, 1, -3) if e + 3 * j >= 0 and seg[e + 3 * j:e + 3 * j + 3] == "CGG"]
            if first:
                n_ev[d] += 1
                prev_codon[(d, seg[e - 3:e])] += 1
                if cand: n_ev_cgg[d] += 1
            if not cand: continue
            p = cand[0]
            if d == 1:
                up = seg[(p % 3):p]; dn = seg[p + 4:]
            else:
                up = seg[(p % 3):p]; dn = seg[p + 2:]
            up = up[max(0, len(up) - 120):]; dn = dn[:120]
            hyp = (tr(up) + "X" + tr(dn)).replace("*", "*")
            xi = len(tr(up))
            # reference window around the event
            near = [r for r in a["recs"] if abs(r[0] - e) <= 6]
            # query positions are not stored in recs; recover via protein index from cols walk
            ref = prots[a["qn"]]
            # approximate reference position: alignment qs + (number of M cols before e)
            k_before = sum(1 for r in a["recs"] if r[0] < e)
            qguess = a["qs"] + k_before
            rwin = ref[max(0, qguess - 60):qguess + 60]
            res, clean, idn, sc = residue_at(hyp, xi, rwin)
            gp = gcoord(a, p, 3)
            sites[(a["contig"], a["strand"], gp, d)].append((a["ref"], res, clean, idn))
            ctx[(a["contig"], a["strand"], gp, d)] = seg[max(0, p - 6):p + 10]
            ctxl[(a["contig"], a["strand"], gp, d)] = seg[max(0, p - 30):p + 33]
    # aggregate
    out = open(f"{W}/out/{lab}.fs_cgg_sites.tsv", "w")
    out.write("contig\tstrand\tcgg_start0\tdir\tcontext(-6..+10, CGG at 7-9)\tnext_nt\tn_refs\tclean_residues\tconsensus\tper_ref\n")
    cons = {1: Counter(), 2: Counter()}; nxt = {1: Counter(), 2: Counter()}; allpairs = {1: Counter(), 2: Counter()}
    seenctx = set(); ndup = 0
    for key, v in sorted(sites.items()):
        if ctxl[key] in seenctx: ndup += 1; continue
        seenctx.add(ctxl[key])
        byref = {}
        for (r, res, clean, idn) in v:
            if r not in byref or idn > byref[r][2]: byref[r] = (res, clean, idn)
        cl = [x[0] for x in byref.values() if x[1] and x[2] >= 0.3 and x[0] not in "-X*"]
        for x in cl: allpairs[key[3]][x] += 1
        c = ""
        if len(cl) >= 2:
            r, k = Counter(cl).most_common(1)[0]
            if k / len(cl) >= 0.6: c = r; cons[key[3]][r] += 1
        cx = ctx[key]
        nxt[key[3]][cx[9:10]] += 1
        out.write(f"{key[0]}\t{key[1]}\t{key[2]}\t{'+1' if key[3]==1 else '-1'}\t{cx[:6].lower()}{cx[6:9]}{cx[9:].lower()}\t{cx[9:10]}\t{len(byref)}\t{''.join(sorted(cl))}\t{c}\t"
                  + ";".join(f"{r}:{x[0]}:{int(x[1])}:{x[2]:.2f}" for r, x in sorted(byref.items())) + "\n")
    out.close()
    s = dict(collapsed_duplicate_sites=ndup, events_p1=n_ev[1], events_m1=n_ev[2], events_p1_with_CGG=n_ev_cgg[1], events_m1_with_CGG=n_ev_cgg[2],
             cgg_sites_p1=sum(1 for k in sites if k[3] == 1), cgg_sites_m1=sum(1 for k in sites if k[3] == 2), unique_sites_p1=sum(cons[1].values()) if False else None,
             consensus_p1=dict(cons[1].most_common()), consensus_m1=dict(cons[2].most_common()),
             pairs_p1=dict(allpairs[1].most_common()), next_nt_p1=dict(nxt[1].most_common()),
             prev_codon_top_p1=[(k[1], v) for k, v in prev_codon.most_common() if k[0] == 1][:8],
             prev_codon_top_m1=[(k[1], v) for k, v in prev_codon.most_common() if k[0] == 2][:5])
    summ[lab] = s
    print(lab, json.dumps(s))
try:
    old = json.load(open(f"{W}/out/fs_detail_summary.json"))
except Exception:
    old = {}
old.update(summ)
json.dump(old, open(f"{W}/out/fs_detail_summary.json", "w"), indent=1)
