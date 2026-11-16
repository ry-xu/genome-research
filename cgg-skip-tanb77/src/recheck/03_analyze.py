"""Per-genome codon statistics from miniprot alignments of standard-code reference proteomes.

usage: 03_analyze.py LABEL [PAF] [GENOME]   (defaults out/LABEL.paf, genomes/LABEL.fna)
Writes out/LABEL.summary.json and out/LABEL.cgg_loci.tsv.

Statistics per focal codon (in-frame = the codon sits in an M column of a reference alignment):
  n_loci       distinct in-frame loci (union over references; local window identity >= 0.3)
  cons         loci at conserved positions: >= 2 references with a clean column (window identity >= 0.4, no
               indel/frameshift within 3 columns, reference residue not X) and >= 60% of them agree;
               counts of the consensus residue
  pairs        all clean (locus, reference) residue pairs
  fs           loci with a frameshift event within 9 nt (3 codons) of the codon, supported by >= half
               (and >= 1) of the references covering the locus; by direction (+1 / -1)
  D            genome codons facing a protein gap (codon skipped by every aligned protein)
  inframe/outframe  occurrences in frame 0 vs frames +1/+2 at the same aligned loci (stop-codon test)
  junction     frameshift events (merged across references) with the codon in the old-frame window
               (5 codons: 3 before the event, 2 from it) or in the new-frame window (5 codons)
"""
import sys, json
from collections import Counter, defaultdict
sys.path.insert(0, "src/recheck")
from mp import read_fasta, parse_paf, gcoord, best_per_ref, rc, STD, CODONS

W = "recheck_work"
lab = sys.argv[1]
paf = sys.argv[2] if len(sys.argv) > 2 else f"{W}/out/{lab}.paf"
gfa = sys.argv[3] if len(sys.argv) > 3 else f"{W}/genomes/{lab}.fna"
outp = sys.argv[4] if len(sys.argv) > 4 else f"{W}/out/{lab}"
genome = read_fasta(gfa)
prots = read_fasta(f"{W}/refs/refs_all.faa")
alns = [a for a in parse_paf(paf, genome, prots) if a["nM"] >= 30 and a["AS"] >= 60]
best = best_per_ref(alns)

# ---- loci
loci = defaultdict(list)  # (contig, strand, start) -> [(ref, raa, wid, neargap, ai, o)]
for (ctg, st, s, ref), (AS, ai, o, codon, raa, wid, neargap, kf, kt) in best.items():
    loci[(ctg, st, s)].append((ref, raa, wid, neargap, ai, o, codon))
loci = {k: v for k, v in loci.items() if max(x[2] for x in v) >= 0.3}

def codon_at(ctg, st, s, shift=0):
    g = genome[ctg]
    if st == "+":
        a = s + shift
        return g[a:a + 3] if 0 <= a and a + 3 <= len(g) else ""
    a = s - shift
    return rc(g[a:a + 3]) if 0 <= a and a + 3 <= len(g) else ""

# ---- frameshift events: only from an alignment that is the best one (for its reference) at the codon before it
ev_raw = []
for ai, a in enumerate(alns):
    for (e, d, n, op) in a["events"]:
        if d == 0: continue
        prev = e - 3
        key = (a["contig"], a["strand"], gcoord(a, prev), a["ref"])
        b = best.get(key)
        if b is None or b[1] != ai: continue
        g = gcoord(a, e, max(n, 1))
        seg = a["seg"]
        old = [seg[e + 3 * j:e + 3 * j + 3] for j in range(-3, 2) if e + 3 * j >= 0]
        new = [seg[e + n + 3 * j:e + n + 3 * j + 3] for j in range(-2, 3) if e + n + 3 * j >= 0]
        ev_raw.append((a["contig"], a["strand"], g, d, a["ref"], old, new, a["qn"], e))
ev_raw.sort()
events = []  # merged: same contig/strand/direction within 9 nt
for x in ev_raw:
    if events and events[-1]["ctg"] == x[0] and events[-1]["st"] == x[1] and events[-1]["d"] == x[3] and x[2] - events[-1]["g1"] <= 9:
        m = events[-1]; m["g1"] = x[2]; m["refs"].add(x[4]); m["old"].update(x[5]); m["new"].update(x[6])
    else:
        events.append(dict(ctg=x[0], st=x[1], g0=x[2], g1=x[2], d=x[3], refs={x[4]}, old=set(x[5]), new=set(x[6]), qn=x[7]))
ev_index = defaultdict(list)
for m in events:
    ev_index[(m["ctg"], m["st"])].append(m)

# ---- collapse repeated elements: a locus (or event) whose oriented +-30 nt context was already seen is a copy
def ctx_at(ctg, st, s, half=30):
    g = genome[ctg]
    w = g[max(0, s - half):s + 3 + half]
    return w if st == "+" else rc(w)
seen = set(); dup_loci = 0
for k in sorted(loci):
    c = ctx_at(*k)
    if c in seen: del loci[k]; dup_loci += 1
    else: seen.add(c)
seen = set(); ev_unique = []
for m in events:
    c = ctx_at(m["ctg"], m["st"], m["g0"])
    if c not in seen: seen.add(c); ev_unique.append(m)
print(lab, "collapsed duplicate-context loci:", dup_loci, "events:", len(events) - len(ev_unique), file=sys.stderr)

def fs_near(ctg, st, s, ncover):
    need = max(1, (ncover + 1) // 2)
    out = Counter()
    for m in ev_index.get((ctg, st), []):
        if m["g1"] + 9 < s - 1 or m["g0"] - 9 > s + 3: continue
        dist = 0 if (m["g0"] <= s + 2 and m["g1"] >= s) else min(abs(m["g0"] - (s + 2)), abs(s - m["g1"]))
        if dist <= 9 and len(m["refs"]) >= need:
            out["+1" if m["d"] == 1 else "-1"] += 1
    return out

# ---- per-codon stats
stats = {c: dict(n_loci=0, cons=Counter(), pairs=Counter(), fs=Counter(), D=0, inframe=0, outframe=0,
                 junction_old=0, junction_new=0, cons_n=0) for c in CODONS}
cgg_rows = []
for (ctg, st, s), v in loci.items():
    codon = v[0][6]
    if codon not in stats: continue
    S = stats[codon]
    S["n_loci"] += 1
    S["inframe"] += 1
    o1, o2 = codon_at(ctg, st, s, 1), codon_at(ctg, st, s, 2)
    for oc in (o1, o2):
        if oc in stats: stats[oc]["outframe"] += 0.5
    clean = [x for x in v if x[2] >= 0.4 and not x[3] and x[1] != "X"]
    for x in clean: S["pairs"][x[1]] += 1
    cons = None
    if len(clean) >= 2:
        r, k = Counter(x[1] for x in clean).most_common(1)[0]
        if k / len(clean) >= 0.6:
            cons = r; S["cons"][r] += 1; S["cons_n"] += 1
    fs = fs_near(ctg, st, s, len(v))
    for k2, n2 in fs.items(): S["fs"][k2] += 1
    if codon == "CGG":
        cgg_rows.append((ctg, st, s, len(v), len(clean), cons or "", "".join(sorted(x[1] for x in clean)),
                         ";".join(f"{x[0]}:{x[1]}:{x[2]:.2f}" for x in v), dict(fs), v[0][4]))
# D columns (genome codon facing protein gap) from best alignments only
seenD = set()
for ai, a in enumerate(alns):
    for (o, cod) in a["Dcols"]:
        key = (a["contig"], a["strand"], gcoord(a, o))
        if key in seenD or key in loci: continue
        seenD.add(key)
        if cod in stats: stats[cod]["D"] += 1
for m in ev_unique:
    for c in m["old"]:
        if c in stats: stats[c]["junction_old"] += 1
    for c in m["new"]:
        if c in stats: stats[c]["junction_new"] += 1

tot = sum(S["n_loci"] for S in stats.values())
summary = dict(label=lab, n_alignments=len(alns), n_loci=tot, n_events=len(ev_unique), n_events_raw=len(events), dup_loci=dup_loci,
               n_events_p1=sum(1 for m in ev_unique if m["d"] == 1), n_events_m1=sum(1 for m in ev_unique if m["d"] == 2),
               n_events_ge2refs=sum(1 for m in ev_unique if len(m["refs"]) >= 2),
               codons={c: dict(n_loci=S["n_loci"], cons_n=S["cons_n"], cons=dict(S["cons"].most_common()),
                               pairs=dict(S["pairs"].most_common()), fs=dict(S["fs"]), D=S["D"],
                               inframe=S["inframe"], outframe=S["outframe"],
                               junction_old=S["junction_old"], junction_new=S["junction_new"])
                       for c, S in stats.items()})
json.dump(summary, open(f"{outp}.summary.json", "w"), indent=1)
with open(f"{outp}.cgg_loci.tsv", "w") as f:
    f.write("contig\tstrand\tstart0\tn_refs\tn_clean\tconsensus\tclean_residues\tper_ref(ref:res:window_id)\tfs_near\tfirst_ref_aa\n")
    for r in sorted(cgg_rows):
        f.write("\t".join(map(str, r)) + "\n")
# events list for inspection
with open(f"{outp}.events.tsv", "w") as f:
    f.write("contig\tstrand\tg0\tg1\tdir\tn_refs\trefs\texample_query\told_frame_codons\tnew_frame_codons\n")
    for m in events:
        f.write(f"{m['ctg']}\t{m['st']}\t{m['g0']}\t{m['g1']}\t{'+1' if m['d']==1 else '-1'}\t{len(m['refs'])}\t{','.join(sorted(m['refs']))}\t{m['qn']}\t{','.join(sorted(m['old']))}\t{','.join(sorted(m['new']))}\n")
c = summary["codons"]["CGG"]
print(lab, "alns", len(alns), "loci", tot, "events", len(events), "| CGG loci", c["n_loci"], "cons", c["cons"], "fs", c["fs"],
      "D", c["D"], "in/out", c["inframe"], c["outframe"], "junc", c["junction_old"], c["junction_new"])
