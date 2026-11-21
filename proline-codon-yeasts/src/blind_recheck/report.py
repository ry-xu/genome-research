"""Assemble report.md from the Codetta outputs (method 1) and the ortholog-column tallies (method 2).

usage: report.py work_dir refsum_dir src_dir out.md
  work_dir   : blind_cun/work (codetta/, genes/, genomes/, refs/, tally*/)
  refsum_dir : download dir holding NCBI esummary JSON (for reference organism names)
  src_dir    : where the scripts live (only printed into the report)
"""
import os
import re
import sys
import csv
import glob
import json
import datetime
from collections import Counter, defaultdict

from scipy.stats import fisher_exact

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import CODONS, T4, CONTESTED, read_fasta  # noqa: E402

PANEL = ["BK074526.1", "BK074532.1", "BK074539.1", "BK074262.1",
         "BK074538.1", "BK074252.1", "BK074542.1", "BK074222.1"]
GENES = ["cox1", "cox2", "cox3", "cob", "atp6", "atp8", "atp9",
         "nad1", "nad2", "nad3", "nad4", "nad4L", "nad5", "nad6"]
CTN = ["CTT", "CTC", "CTA", "CTG"]
CCN = ["CCT", "CCC", "CCA", "CCG"]
TTR = ["TTA", "TTG"]
MIN_N = 5        # conserved columns needed to call a codon by method 2
MIN_TOP = 0.8    # top-residue share needed to call it
N_EX = 8         # example positions per deviating codon


# ---------------------------------------------------------------- io helpers
def md_table(header, rows):
    out = ["| " + " | ".join(header) + " |", "|" + "|".join("---" for _ in header) + "|"]
    for r in rows:
        out.append("| " + " | ".join(str(x) for x in r) + " |")
    return "\n".join(out)


def short_org(title):
    t = re.sub(r"^TPA_asm:\s*", "", title.strip())
    w = t.split()
    return " ".join(w[:2]) if len(w) >= 2 else t


def parse_codetta(path):
    inf, logp, final, cols, sec = {}, {}, None, None, None
    for line in open(path):
        if line.startswith("# Codon inferences"):
            sec = "inf"
            continue
        if line.startswith("# Log decoding"):
            sec = "logp"
            continue
        if line.startswith("# Final genetic code"):
            sec = "final"
            continue
        parts = line.split()
        if sec == "logp" and line.startswith("# codon"):
            cols = [c[5:-1] for c in parts[2:]]
            continue
        if line.startswith("#") or not parts:
            continue
        if sec == "inf" and len(parts) >= 6 and len(parts[0]) == 3:
            inf[parts[0]] = dict(inf=parts[1], N_al=int(parts[4]), N_used=int(parts[5]))
        elif sec == "logp" and len(parts[0]) == 3 and cols:
            logp[parts[0]] = dict(zip(cols, map(float, parts[1:])))
        elif sec == "final":
            final = parts[0]
    return inf, logp, final


def read_codon_table(path):
    d = {}
    for r in csv.DictReader(open(path), delimiter="\t"):
        dist = Counter()
        for x in (r["dist"] or "").split(","):
            if x:
                k, v = x.split(":")
                dist[k] = int(v)
        d[(r["acc"], r["filter"], r["codon"])] = dict(usage=int(r["usage"]), dist=dist)
    return d


def dist_of(tab, acc, codons, filt="mask"):
    c = Counter()
    for cod in codons:
        c.update(tab.get((acc, filt, cod), {}).get("dist", Counter()))
    return c


def fmt_dist(c, k=4):
    if not c:
        return "-"
    return ", ".join(f"{a}{n}" for a, n in c.most_common(k)) + (" ..." if len(c) > k else "")


def fmt_p(p):
    return "1" if p >= 0.9995 else (f"{p:.3f}" if p >= 0.001 else f"{p:.1e}")


# ---------------------------------------------------------------- method-2 helpers
def call_m2(c, t4aa):
    N = sum(c.values())
    if N == 0:
        return "no data", N, "-", 0.0
    top, n = c.most_common(1)[0]
    f = n / N
    if N >= MIN_N and f >= MIN_TOP:
        return ("DEVIATES" if top != t4aa else "t4"), N, top, f
    if N >= 2 and top != t4aa and f == 1.0:
        return "tentative-deviation", N, top, f
    if top == t4aa and f >= 0.5:
        return "t4 (weak)" if N < MIN_N else "t4 (mixed)", N, top, f
    return "unresolved", N, top, f


def aln_pos(row, col):
    """residue index (0-based) of column col in an aligned row, and whether col is a gap."""
    k = sum(1 for ch in row[:col] if ch != "-")
    return k, row[col] == "-"


def peptide(seq, k, gap, w=5, lower=None):
    lower = lower or set()
    def r(j):
        a = seq[j]
        return a.lower() if j in lower else a
    left = "".join(r(j) for j in range(max(0, k - w), k))
    if gap:
        right = "".join(r(j) for j in range(k, min(len(seq), k + w)))
        return f"{left}[-]{right}"
    right = "".join(r(j) for j in range(k + 1, min(len(seq), k + w + 1)))
    return f"{left}[{r(k)}]{right}"


def main():
    work, refsum_dir, src_dir, out_md = sys.argv[1:5]
    tally_dir = os.path.join(work, "tally")
    tab = read_codon_table(os.path.join(tally_dir, "codon_table.tsv"))
    sens = {name: read_codon_table(os.path.join(work, name, "codon_table.tsv"))
            for name in ("tally_n7_f0.9", "tally_n5_f0.6", "tally_n3_f0.8")}
    cons = list(csv.DictReader(open(os.path.join(tally_dir, "conserved_records.tsv")), delimiter="\t"))

    # organism names
    panel_org = {}
    for acc in PANEL:
        hdr = read_fasta(os.path.join(work, "genomes", f"{acc}.fna"))[0][0]
        panel_org[acc] = short_org(hdr.split(None, 1)[1] if " " in hdr else hdr)
    ref_org = {}
    for f in glob.glob(os.path.join(refsum_dir, "*.json")):
        try:
            d = json.load(open(f))
        except Exception:
            continue
        for v in (d.get("result", {}) if isinstance(d, dict) else {}).values():
            if isinstance(v, dict) and v.get("caption"):
                ref_org[v["caption"]] = short_org(v.get("title", ""))
    refsel = list(csv.DictReader(open(os.path.join(work, "refs", "refs_selected.tsv")), delimiter="\t"))
    ref_genes = defaultdict(list)
    for r in refsel:
        ref_genes[r["acc"]].append(r["gene"])

    # gene models
    models = {acc: json.load(open(os.path.join(work, "genes", f"{acc}.genes.json"))) for acc in PANEL}

    # codetta
    cdt, cdt_r = {}, {}
    for acc in PANEL:
        p = glob.glob(os.path.join(work, "codetta", f"{acc}.fna.Pfam-A_enone.hmm.1e-10_0.9999_0.01_*.genetic_code.out"))
        cdt[acc] = parse_codetta(p[0])
        pr = os.path.join(work, "codetta", f"{acc}.relaxed_r0.99.genetic_code.out")
        cdt_r[acc] = parse_codetta(pr) if os.path.exists(pr) else None

    def cdt_dev(res):
        inf = res[0]
        return [(c, inf[c]["inf"]) for c in CODONS if inf[c]["inf"] not in ("?", T4[c])]

    # method-2 calls
    calls = {}
    for acc in PANEL:
        for c in CODONS:
            calls[(acc, c)] = call_m2(tab.get((acc, "mask", c), {}).get("dist", Counter()), T4[c])

    # 2x2 tests
    def two_by_two(acc, test, ctrl, a, b):
        ct, cc = dist_of(tab, acc, test), dist_of(tab, acc, ctrl)
        ta, tb = ct.get(a, 0), (sum(v for k, v in ct.items() if k != a) if b is None else ct.get(b, 0))
        ca, cb = cc.get(a, 0), (sum(v for k, v in cc.items() if k != a) if b is None else cc.get(b, 0))
        orr, p = fisher_exact([[ta, tb], [ca, cb]]) if (ta + tb) and (ca + cb) else (float("nan"), 1.0)
        return ta, tb, ca, cb, orr, p

    L = []
    w = L.append
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M %Z").strip()
    w("# Blinded codon-meaning inference: 8 Dipodascomycetes mitogenomes (NCBI TPA)\n")
    w(f"Generated {now} by `{os.path.join(src_dir, 'report.py')}`. Null hypothesis: NCBI table 4 "
      "(the candidate alternative was table 3). Each codon's meaning was inferred from the data with two "
      "independent methods: (1) Codetta, which uses Pfam profiles and no relatives; (2) an ortholog-column "
      "tally against curated relatives. Organism names come from the FASTA headers and are listed for "
      "orientation only; they were not used in any inference.\n")

    # ------------------------------------------------------------ summary table
    w("## 1. Per-genome summary\n")
    pg = [a for a in PANEL if any(calls[(a, c)][0] == "DEVIATES" for c in CODONS)]
    w(f"**Bottom line.** Four genomes ({', '.join(pg)}) read CTT and CTA as **proline**, not leucine. Both methods "
      "support this independently, and no other codon departs from table 4 in any genome. In these four genomes, CTC "
      "and CTG are nearly unused (0-1 codons each), so their meaning cannot be determined. The other four genomes fit "
      "table 4. BK074539.1 is a strong CUN = Leu control. BK074526.1, BK074252.1 and BK074222.1 have few CUN codons "
      "but read them as Leu wherever they can be scored, never as Pro. Every genome reads TGA as Trp and ATA as Ile, "
      "which is table 4 and rules out table 3.\n")
    rows = []
    for acc in PANEL:
        devs = [c for c in CODONS if calls[(acc, c)][0] == "DEVIATES"]
        tent = [c for c in CODONS if calls[(acc, c)][0] == "tentative-deviation"]
        cd = cdt_dev(cdt[acc])
        cdr = cdt_dev(cdt_r[acc]) if cdt_r[acc] else []
        ncod = sum(len(m["codons"]) for m in models[acc].values())
        if devs:
            dev_txt = ", ".join(f"{c}→{calls[(acc, c)][2]}" for c in devs)
            m2 = "; ".join(f"{c}: {calls[(acc, c)][2]} at {calls[(acc, c)][1] and round(calls[(acc, c)][3] * calls[(acc, c)][1])}/{calls[(acc, c)][1]} conserved columns" for c in devs)
            m1 = ", ".join(f"{c}={a}" for c, a in cd if c != "TGA") or "none"
            fx = [two_by_two(acc, [c], TTR, "P", "L") for c in devs]
            m2 += " (Fisher vs TTR, P vs L: " + ", ".join(f"{c} p={fmt_p(f[5])}" for c, f in zip(devs, fx)) + ")"
            unres = [c for c in CTN if c not in devs]
            un_txt = ", ".join(f"{c} used {tab[(acc, 'mask', c)]['usage']}x / {calls[(acc, c)][1]} cons. cols" for c in unres)
            inferred = f"table 4 except {' and '.join(devs)} = Pro; {'/'.join(unres)} not callable ({un_txt})"
        else:
            dev_txt = "none detected"
            m1 = ", ".join(f"{c}={a}" for c, a in cd if c != "TGA") or "no deviation (only TGA=W, which is table 4)"
            cun = dist_of(tab, acc, CTN)
            m2 = f"no codon deviates; CUN at conserved columns: {fmt_dist(cun)} (n={sum(cun.values())})"
            if tent:
                m2 += "; single-column oddities: " + ", ".join(f"{c}→{calls[(acc, c)][2]} (n={calls[(acc, c)][1]})" for c in tent)
            inferred = "table 4"
        cun_m1 = "/".join(cdt[acc][0][c]["inf"] for c in CTN)
        cun_m1r = "/".join(cdt_r[acc][0][c]["inf"] for c in CTN) if cdt_r[acc] else "-"
        rows.append((f"{acc} ({panel_org[acc]})", f"{len(models[acc])} genes, {ncod} codons", dev_txt, inferred,
                     f"default: CTT/CTC/CTA/CTG = {cun_m1}; r=0.99: {cun_m1r}", m2))
    w(md_table(["genome", "modelled", "deviations from table 4", "inferred code",
                "method 1 (Codetta, CUN shown)", "method 2 (ortholog columns)"], rows))
    w("")
    w("**Codons checked in all 8 genomes.** TGA = Trp: method 2 gives "
      + ", ".join(f"{acc[:8]} {dist_of(tab, acc, ['TGA']).get('W', 0)}/{sum(dist_of(tab, acc, ['TGA']).values())}" for acc in PANEL)
      + "; Codetta calls TGA = W in every genome. That is a table-4 meaning, so it is not a deviation. "
      "ATA = Ile, not Met (table 3 would read Met): pooled conserved columns give I "
      f"{dist_of(tab, 'ALL', ['ATA']).get('I', 0)}, M {dist_of(tab, 'ALL', ['ATA']).get('M', 0)} "
      f"of {sum(dist_of(tab, 'ALL', ['ATA']).values())}, and Codetta calls ATA = I in all 8. "
      "TAA and TAG never occur inside the modelled genes. The remaining codons are either table 4 or too rare to call "
      "(see section 6). Table 3's CUN = Thr is not supported anywhere: pooled over all genomes, CUN codons sit at "
      f"{dist_of(tab, 'ALL', CTN).get('T', 0)} conserved Thr columns out of {sum(dist_of(tab, 'ALL', CTN).values())}.\n")

    # ------------------------------------------------------------ 2x2 tables
    w("## 2. Method 2 internal controls: 2×2 tests against table-4 synonyms\n")
    w("Each cell counts conserved ortholog columns (≥5 relatives with a residue, ≥80% of them identical; "
      "relatives' CUN/AUA/CGA/CGC codons and stops masked) whose consensus is the named residue. "
      "Fisher's exact test is two-sided. The control row holds the codon's unambiguous table-4 synonyms, "
      "so a codon that means what table 4 says should look like its control.\n")
    tests = [
        ("CTT vs TTA+TTG", ["CTT"], TTR, "P", "L"),
        ("CTA vs TTA+TTG", ["CTA"], TTR, "P", "L"),
        ("CTN (all 4) vs TTA+TTG", CTN, TTR, "P", "L"),
        ("CTT vs all other table-4 Leu codons (TTR, CTC, CTA, CTG)", ["CTT"], TTR + ["CTC", "CTA", "CTG"], "P", "L"),
        ("CTA vs all other table-4 Leu codons (TTR, CTT, CTC, CTG)", ["CTA"], TTR + ["CTT", "CTC", "CTG"], "P", "L"),
        ("comparator: CCN (Pro) vs TTA+TTG", CCN, TTR, "P", "L"),
        ("ATA vs ATT+ATC (table 3 would make ATA Met)", ["ATA"], ["ATT", "ATC"], "M", "I"),
        ("TGA vs TGG (Trp vs any other residue)", ["TGA"], ["TGG"], "W", None),
    ]
    for title, test, ctrl, a, b in tests:
        bl = b or f"non-{a}"
        w(f"**{title}**: columns {a} vs {bl}\n")
        rows = []
        for acc in PANEL:
            ta, tb, ca, cb, orr, p = two_by_two(acc, test, ctrl, a, b)
            rows.append((acc, f"{ta} / {tb}", f"{ca} / {cb}", "inf" if orr == float("inf") else (f"{orr:.2g}" if orr == orr else "-"), fmt_p(p),
                         fmt_dist(dist_of(tab, acc, test))))
        w(md_table(["genome", f"test {a} / {bl}", f"control {a} / {bl}", "odds ratio", "Fisher p",
                    "test codon(s): all consensus residues"], rows))
        w("")

    # ------------------------------------------------------------ Pro/Leu encoding
    w("## 3. Which codons fill conserved Pro and Leu columns (method 2)\n")
    w("In the four CUN→Pro genomes, conserved Pro columns are filled mainly by CUN and conserved Leu columns "
      "almost only by TTR. Leu is recoded onto TTA/TTG rather than lost, and Pro is carried by CUN.\n")
    rows = []
    for acc in PANEL:
        pc = Counter(r["codon"] for r in cons if r["acc"] == acc and r["cons_aa"] == "P")
        lc = Counter(r["codon"] for r in cons if r["acc"] == acc and r["cons_aa"] == "L")
        rows.append((acc, sum(pc[c] for c in CCN), sum(pc[c] for c in CTN), sum(pc.values()) - sum(pc[c] for c in CCN + CTN),
                     sum(lc[c] for c in TTR), sum(lc[c] for c in CTN), sum(lc.values()) - sum(lc[c] for c in TTR + CTN),
                     " ".join(f"{c}:{tab[(acc, 'mask', c)]['usage']}" for c in CCN),
                     " ".join(f"{c}:{tab[(acc, 'mask', c)]['usage']}" for c in CTN)))
    w(md_table(["genome", "Pro cols: CCN", "Pro cols: CUN", "Pro cols: other", "Leu cols: TTR", "Leu cols: CUN",
                "Leu cols: other", "CCN usage (all modelled codons)", "CUN usage"], rows))
    w("")

    # ------------------------------------------------------------ examples
    w("## 4. Example positions for every deviating codon\n")
    w("`codon #` is the codon's 1-based index in the panel gene model. The model is built from tblastn HSPs, so its "
      "first codon may not be the annotated start. `pos` is the 1-based genome coordinate span of the codon. "
      "Peptides are ±5 residues, focal residue in brackets. The panel peptide is translated with table 4, so the "
      "focal CUN prints as `l`. Lowercase marks a contested codon (CUN/AUA/CGA/CGC) shown with its table-4 letter. "
      "The relatives are the two most identical to the panel protein in that gene's alignment that have a residue in "
      "the column. `column` is the masked consensus of all relatives at that column.\n")
    for acc in PANEL:
        devs = [c for c in CODONS if calls[(acc, c)][0] == "DEVIATES"]
        for cod in devs:
            top = calls[(acc, cod)][2]
            recs = [r for r in cons if r["acc"] == acc and r["codon"] == cod and r["cons_aa"] == top]
            recs.sort(key=lambda r: (-float(r["frac"]), -int(r["n_ref"])))
            by_gene = defaultdict(list)
            for r in recs:
                by_gene[r["gene"]].append(r)
            picked = []
            while len(picked) < N_EX and any(by_gene.values()):
                for g in sorted(by_gene, key=lambda g: -len(by_gene[g])):
                    if by_gene[g] and len(picked) < N_EX:
                        picked.append(by_gene[g].pop(0))
            w(f"### {acc} ({panel_org[acc]}) {cod} → {top}: {len(recs)} of {calls[(acc, cod)][1]} conserved columns "
              f"(used {tab[(acc, 'mask', cod)]['usage']}× in the modelled genes; genes: "
              + ", ".join(f"{g} {n}" for g, n in sorted(Counter(r['gene'] for r in recs).items())) + ")\n")
            rows = []
            for r in picked:
                g, idx, col = r["gene"], int(r["idx"]), int(r["col"])
                m = models[acc][g]
                cods = [c[1] for c in m["codons"]]
                pseq = "".join(T4.get(c, "X") for c in cods)
                plow = {j for j, c in enumerate(cods) if c in CONTESTED}
                ppep = peptide(pseq, idx - 1, False, lower=plow)
                gp = int(r["gpos0"])
                pos = f"{gp + 1}-{gp + 3} (+)" if m["strand"] == "+" else f"{gp - 1}-{gp + 1} (−)"
                aln = dict(read_fasta(os.path.join(tally_dir, "aln", f"{acc}.{g}.mask.aln")))
                rt4 = dict(read_fasta(os.path.join(work, "refs", f"ref_{g}.t4.faa")))
                rmk = dict(read_fasta(os.path.join(work, "refs", f"ref_{g}.mask.faa")))
                prow = aln[f"PANEL|{acc}"]
                sims = []
                for name, row in aln.items():
                    if name.startswith("PANEL") or row[col] == "-":
                        continue
                    same = tot = 0
                    for x, y in zip(prow, row):
                        if x in "-X" or y in "-Xx":
                            continue
                        tot += 1
                        same += x.upper() == y.upper()
                    sims.append((same / tot if tot else 0.0, name))
                sims.sort(reverse=True)
                rel = []
                for ident, name in sims[:2]:
                    k, gap = aln_pos(aln[name], col)
                    low = {j for j, ch in enumerate(rmk[name]) if ch == "X" and rt4[name][j] != "X"}
                    ra = name.split("|")[0]
                    rel.append(f"`{peptide(rt4[name], k, gap, lower=low)}` {ra} {ref_org.get(ra, '')} ({ident:.0%} id)")
                while len(rel) < 2:
                    rel.append("-")
                colc = Counter(ch for ch in r["ref_res"] if ch != "-")
                rows.append((g, idx, pos, f"`{ppep}`", rel[0], rel[1], " ".join(f"{a}{n}" for a, n in colc.most_common())))
            w(md_table(["gene", "codon #", "pos", "panel ±5", "relative 1", "relative 2", "column"], rows))
            w("")

    # ------------------------------------------------------------ codetta detail
    w("## 5. Method 1: Codetta per-codon results\n")
    w("Codetta v2.0 run as cloned (`codetta.py <genome>.fna -m`, so mitochondrial Pfams are included; "
      "defaults e-value 1e-10, probability threshold 0.9999, max_fraction 0.01). A second pass "
      "(`codetta_infer.py -m -r 0.99`) reuses the same alignments with a relaxed probability threshold. "
      "`N` = Pfam consensus columns where the codon aligned. logP values are Codetta's log decoding probabilities; "
      "0 marks the most probable meaning, and the more negative a value, the less probable that meaning. "
      "Codetta compares against the standard code, so its TGA = W 'difference' is table 4's own meaning.\n")
    rows = []
    for acc in PANEL:
        inf, lp, final = cdt[acc]
        infr = cdt_r[acc][0] if cdt_r[acc] else {}
        dv = cdt_dev(cdt[acc])
        dvr = cdt_dev(cdt_r[acc]) if cdt_r[acc] else []
        nq = sum(1 for c in CODONS if inf[c]["inf"] == "?" and c not in ("TAA", "TAG"))
        tot = sum(inf[c]["N_al"] for c in CODONS)
        rows.append((acc, tot, nq, ", ".join(f"{c}={a}" for c, a in dv) or "none",
                     ", ".join(f"{c}={a}" for c, a in dvr) or "none", f"`{final}`"))
    w(md_table(["genome", "total aligned codon-columns", "'?' sense codons (default)", "deviations from table 4 (default)",
                "deviations from table 4 (r=0.99)", "default inferred code (TCAG order)"], rows))
    w("")
    w("Key codons. Each cell gives `inference [relaxed inference] N; logP(L) / logP(P) / logP(T)` for CUN, "
      "`logP(I) / logP(M)` for ATA, and `logP(W)` for TGA.\n")
    rows = []
    for acc in PANEL:
        inf, lp, _ = cdt[acc]
        infr = cdt_r[acc][0] if cdt_r[acc] else {}
        cells = [acc]
        for c in CTN:
            cells.append(f"{inf[c]['inf']} [{infr.get(c, {}).get('inf', '-')}] N={inf[c]['N_al']}; "
                         f"{lp[c]['L']:.1f} / {lp[c]['P']:.1f} / {lp[c]['T']:.1f}")
        c = "ATA"
        cells.append(f"{inf[c]['inf']} [{infr.get(c, {}).get('inf', '-')}] N={inf[c]['N_al']}; {lp[c]['I']:.1f} / {lp[c]['M']:.1f}")
        c = "TGA"
        cells.append(f"{inf[c]['inf']} [{infr.get(c, {}).get('inf', '-')}] N={inf[c]['N_al']}; {lp[c]['W']:.1f}")
        rows.append(cells)
    w(md_table(["genome", "CTT", "CTC", "CTA", "CTG", "ATA", "TGA"], rows))
    w("")
    w("Per-codon Codetta inference with N aligned for all 64 codons (default run):\n")
    rows = []
    for c in CODONS:
        rows.append([c, T4[c]] + [f"{cdt[acc][0][c]['inf']} ({cdt[acc][0][c]['N_al']})" for acc in PANEL])
    w(md_table(["codon", "t4"] + [a[:8] for a in PANEL], rows))
    w("")

    # ------------------------------------------------------------ method-2 per codon
    w("## 6. Method 2 per-codon consensus at conserved columns (all 64 codons)\n")
    w("Each cell gives `usage / N conserved columns: top residues`. A codon is called when N ≥ 5 and its top residue "
      "holds ≥ 80% of the columns. A **bold** cell is a codon whose call differs from table 4.\n")
    rows = []
    for c in CODONS:
        cells = [c, T4[c]]
        for acc in PANEL:
            d = tab.get((acc, "mask", c), {"usage": 0, "dist": Counter()})
            N = sum(d["dist"].values())
            txt = f"{d['usage']} / {N}: {fmt_dist(d['dist'], 3)}" if N else f"{d['usage']} / 0"
            if calls[(acc, c)][0] == "DEVIATES":
                txt = f"**{txt}**"
            cells.append(txt)
        rows.append(cells)
    w(md_table(["codon", "t4"] + [a[:8] for a in PANEL], rows))
    w("")

    # ------------------------------------------------------------ sensitivity
    w("## 7. Method 2 sensitivity: CTT and CTA under other settings\n")
    w("Each cell gives `Pro / Leu / other (N)` at conserved columns. The settings are: primary (masked panel protein in "
      "MAFFT, ≥5 relatives, ≥80% identity); `t4` (panel protein aligned unmasked, as table 4); `hsp-agree` (only codons "
      "where the MAFFT column equals the tblastn HSP column); and alternative conservation thresholds.\n")
    settings = [("primary", tab, "mask"), ("t4 panel", tab, "t4"), ("hsp-agree", tab, "mask_hspagree"),
                ("n≥7, ≥90%", sens["tally_n7_f0.9"], "mask"), ("n≥5, ≥60%", sens["tally_n5_f0.6"], "mask"),
                ("n≥3, ≥80%", sens["tally_n3_f0.8"], "mask")]
    rows = []
    for acc in PANEL:
        for c in ("CTT", "CTA"):
            cells = [acc, c]
            for _, t, f in settings:
                d = t.get((acc, f, c), {"dist": Counter()})["dist"]
                N = sum(d.values())
                cells.append(f"{d.get('P', 0)}/{d.get('L', 0)}/{N - d.get('P', 0) - d.get('L', 0)} ({N})")
            rows.append(cells)
    w(md_table(["genome", "codon"] + [s[0] for s in settings], rows))
    w("")

    # ------------------------------------------------------------ usage
    w("## 8. Codon usage per genome (all 64 codons, modelled genes only)\n")
    w("These are counts over the codons in the gene models (core genes cox1-3, cob, atp6/8/9, nad1-6, nad4L; "
      "intron ORFs excluded; segment junction codons trimmed).\n")
    rows = []
    for c in CODONS:
        rows.append([c, T4[c]] + [tab[(acc, "mask", c)]["usage"] for acc in PANEL])
    rows.append(["total", ""] + [sum(tab[(acc, "mask", c)]["usage"] for c in CODONS) for acc in PANEL])
    w(md_table(["codon", "t4"] + [a[:8] for a in PANEL], rows))
    w("")
    w("Gene models (codons per gene):\n")
    rows = []
    for acc in PANEL:
        rows.append([acc] + [len(models[acc][g]["codons"]) if g in models[acc] else "-" for g in GENES])
    w(md_table(["genome"] + GENES, rows))
    w("")

    # ------------------------------------------------------------ references
    w("## 9. Reference set (method 2)\n")
    w("These are curated relatives with GenBank CDS annotation. Core-gene CDSs were selected by gene name; intron ORFs, "
      "maturases, endonucleases and partial CDSs were dropped, and the longest copy was kept. They were translated with "
      "table 4 to build the MAFFT L-INS-i alignments. For the tally, each relative's CUN/AUA/CGA/CGC codons and stops "
      "were masked to X, so a relative's own code (table 3 vs 4, or a CUN reassignment) cannot leak into the consensus.\n")
    rows = [(a, ref_org.get(a, ""), len(gs), ",".join(gs)) for a, gs in sorted(ref_genes.items())]
    w(md_table(["accession", "organism", "genes", "which"], rows))
    w("")

    # ------------------------------------------------------------ methods, commands, caveats
    recs = list(csv.DictReader(open(os.path.join(tally_dir, "records.tsv")), delimiter="\t"))
    mk = [r for r in recs if r["mode"] == "mask"]
    nh = sum(1 for r in mk if r["hsp_col"] != "-1")
    na = sum(1 for r in mk if r["hsp_col"] == r["col"])
    w("## 10. Scripts, data and commands\n")
    w(f"All scripts are in `{src_dir}/`. They run as `uv run --no-project [--with scipy] python -I <script> ...`, and "
      "every data path is passed as an argument. Working data live in `" + work + "/`.\n")
    w(md_table(["step", "script / command", "output"], [
        ("download panel and references", "NCBI E-utilities `efetch` (rettype fasta / fasta_cds_na) and `esearch`/`esummary`, each into its own `mktemp -d blind_cun/dl.XXXX` directory, ≤3 requests/s", "`dl.*/panel.fna`, `dl.*/refs_cds_na.fna`, `dl.*/*.json`"),
        ("select reference core-gene CDSs", "`prep_refs.py refs_cds_na.fna work/refs ACC,...`", "`work/refs/ref_<gene>.{t4,mask}.faa`, `refs_selected.tsv`"),
        ("reference MSAs", "`build_ref_msa.py work/refs work/msa 4` (MAFFT L-INS-i)", "`work/msa/ref_<gene>.msa`"),
        ("panel gene models", "`find_genes.py work/refs work/msa work/genomes/<acc>.fna work/genes` (tblastn -db_gencode 4 -seg no, masked queries)", "`work/genes/<acc>.{genes.json,models.tsv,hsps.tsv}`"),
        ("align panel to MSAs, record columns", "`align_tally.py work/refs work/msa work/genes work/tally ACC,... 4` (mafft --add --localpair --maxiterate 1000; masked and t4 modes)", "`work/tally/records.tsv`, `work/tally/aln/`"),
        ("per-codon tallies", "`tally_codons.py records.tsv ref_order.json outdir [min_nref] [min_frac]` (5 0.8 primary; 7 0.9, 5 0.6, 3 0.8 sensitivity)", "`work/tally*/codon_table.tsv`, `conserved_records.tsv`"),
        ("Codetta (method 1)", "`HMMER_NCPU=4 uv run --no-project --python 3.12 --with 'numpy<2' --with 'scipy<1.14' python -E -s dl.bllj/codetta/codetta.py work/codetta/<acc>.fna -m`, then `codetta_infer.py <acc>.fna -m -r 0.99 --inference_output <acc>.relaxed_r0.99.genetic_code.out`", "`work/codetta/<acc>.fna.Pfam-A_enone.hmm.1e-10_0.9999_0.01_excl-tvuy.genetic_code.out`, `*.relaxed_r0.99.genetic_code.out`"),
        ("this report", "`report.py work refsum_dir src report.md` (needs scipy)", "`report.md`"),
    ]))
    w("")
    w("## 11. Caveats\n")
    cav = [
        "**Codetta setup differed from stock.** The bundled HMMER 3.3.2 does not configure on Apple-silicon arm64 (it needs SSE). "
        "Codetta therefore ran with Homebrew HMMER 3.4 (`hmmscan`/`hmmpress`/`hmmsearch`) symlinked into its expected `hmmer-3.3.2/bin`. "
        "The easel miniapps (`esl-sfetch`, `esl-translate`) were built from the bundled 3.3.2 easel source after running `autoconf`. "
        "Codetta needed numpy<2 and ran on Python 3.12. Pfam-A_enone was downloaded from eddylab.org/publications/Shulgina21/Pfam-A_enone.tar.gz. "
        "Codetta ran as cloned with `-m`, which includes the mitochondrial Pfams that it excludes by default; that suits mitogenomes. "
        "Setup took well under the 45-minute budget.",
        "**Codetta is underpowered for rare codons in mitogenomes.** Each genome has only about 2,200-2,650 aligned codon-columns. "
        "Rare codons (CTC, CTG, CCC/CCG, CGN, AGG and others) therefore come out '?', and the per-codon logP and N are given in section 5. "
        "Where CUN is rare (BK074526.1, BK074252.1, BK074222.1), Codetta leaves CUN as '?'. For CTT, logP(L) is near 0 in all "
        "three while logP(P) is -19 to -37. CTA has logP(P) of -11 to -19, except BK074526.1, where CTA aligned in only 1 "
        "column and is uninformative. The relaxed run calls CTT = L in BK074222.1 only.",
        "**Low CUN counts in three genomes.** BK074526.1 has 12 CUN codons and 4 scorable columns, BK074252.1 has 12 and 7, and "
        "BK074222.1 has 29 and 14. All scorable columns are Leu (12/14 in BK074222.1, with 0 Pro). This rules out a wholesale "
        "CUN→Pro reading in these genomes. It does not prove every CUN codon is Leu, and CTC/CTG are essentially unobserved there too.",
        "**CTC and CTG are uncalled in the four CUN→Pro genomes**, where each is used 0-1 times. The data show only that CTT and CTA "
        "are Pro. Whether CTC/CTG would also be Pro, or are unassigned, cannot be inferred from these sequences.",
        "**How Pro is encoded varies.** BK074538.1 uses no CCN codons at all in the modelled genes, BK074532.1 uses 3 and "
        "BK074542.1 uses 2, so CUN is effectively their only Pro codon family. BK074262.1 keeps CCN as well: its conserved Pro "
        "columns are 49 CCN and 56 CUN. Leu in all four is encoded almost only by TTR. Genome-wide CC/GG dinucleotides are not "
        "depleted in these genomes (checked separately), so a sequence-assembly artifact that turns CC into CT is unlikely.",
        "**Reference choice and masking.** There are 11 relatives with GenBank CDS annotation (section 9). The closest is 80-87% identical to "
        "the panel proteins in the example genes, and some lack genes (section 9). Relatives were translated with table 4 only to build the "
        "alignments. Their own CUN/AUA/CGA/CGC codons and stops were masked out of every consensus. A Pro consensus therefore "
        "comes only from the relatives' CCN codons, and the relatives' own code cannot create the signal. The panel protein was "
        "also masked at those codons before alignment (primary). An unmasked table-4 alignment ('t4 panel'), a check against "
        "tblastn HSP columns, and three conservation thresholds give the same calls (section 7).",
        f"**HSP column cross-check.** MAFFT columns agree with the tblastn HSP columns for {na}/{nh} ({na / nh:.1%}) of panel "
        "codons in masked mode. A key-mismatch bug at first left this cross-check empty. It was fixed, the alignments were "
        "reused unchanged, and the calls did not change.",
        "**Gene models are approximate.** They are chained tblastn HSPs, not annotated CDSs. Starts and stops are not modelled, "
        "intron ORFs are excluded, and a few codons at exon junctions are dropped where segment ranges overlap. All "
        f"{sum(len(m['codons']) for a in PANEL for m in models[a].values())} model codons were checked against the genome "
        "sequence at their reported coordinates and all match. Example flanking peptides follow the spliced model. For 4 of "
        "the 64 example rows, the ±5 window crosses an exon junction, so it differs at its edge from the contiguous genome "
        "translation; the focal codon is correct in all 64.",
        "**ATA.** At conserved columns the Ile share for ATA (pooled ~0.55) is lower than for ATT/ATC (~0.72). ATA sits at "
        "variable hydrophobic positions more often, where the consensus is L or V. Met appears at only 4 of 536 pooled columns, "
        "and Codetta gives logP(M) ≤ -89 in every genome. ATA is Ile (table 4), not Met (table 3).",
        "**Other low-consensus cells**, such as TCT in BK074526.1, have their top residue equal to the table-4 meaning and reflect "
        "variable columns. In section 6, no codon other than CTT/CTA has a non-table-4 top residue at 2 or more conserved columns.",
        "**No tRNA analysis.** aragorn, tRNAscan-SE and Infernal are not installed and were not installed for this run, so anticodon evidence (for example a "
        "tRNA-Pro(UAG)) was not examined. The inference rests on protein-level evidence only.",
        "**Blinding.** Organism names were visible in the panel FASTA headers (genus names Wickerhamiella, Sugiyamaella, Starmerella) "
        "and in the reference esummary records. No external information on the genetic codes of these taxa was consulted. Nothing under the "
        "excluded lab directories was read. Early on, one `ls` of the scratchpad's top level showed other files' names; none "
        "were opened.",
    ]
    for c in cav:
        w(f"- {c}")
    w("")

    with open(out_md, "w") as fh:
        fh.write("\n".join(L) + "\n")
    print("wrote", out_md, len(L), "blocks")


if __name__ == "__main__":
    main()
