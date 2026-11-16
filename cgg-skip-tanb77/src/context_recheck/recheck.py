"""Independent re-check of the skip_context_chicken count (own code, written from the rule text only).

Pass 1 (this file, `extract`): one row per site (Set A and Set C) and per random in-frame CDS position,
with the upstream codon, recomputed class features, and per-genome background tables.
Pass 2 (`analyse.py`) reads those tables.

Coordinates (found by probing, then checked for every site below): sites.tsv `cpos` is the 0-based forward
index of the C; strand + -> codon = s[c:c+3]; strand - -> codon = revcomp(s[c-2:c+1]).
Everything is done on the CDS strand: for - strand sites/CDSs the contig is reverse-complemented and
p = L-1-c, so codon = r[p:p+3] and the upstream codon = r[p-3:p].
"""
import csv, sys, random, collections, json, os

BASE = ""
GEN = BASE + "data/rna_chicken/genomes/"
SITES = BASE + "results/rna_chicken/sites.tsv"
CLASS = BASE + "results/rna_chicken/site_class.tsv"
EXCL = BASE + "scratch/cgg_test_T.tsv"
OUT = BASE + "results/skip_context_chicken/recheck/out/"
os.makedirs(OUT, exist_ok=True)

STOP = {"TAA", "TAG", "TGA"}
COMP = str.maketrans("ACGTNRYKMSWBDHV", "TGCANYRMKSWVHDB")
CAP = 60
N_RANDOM = 2000  # random in-frame CDS positions per genome (break test b)
AA = {}
_b = "TCAG"
_aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for i, a in enumerate(_b):
    for j, b in enumerate(_b):
        for k, c in enumerate(_b):
            AA[a + b + c] = _aa[16 * i + 4 * j + k]


def rc(s):
    return s.translate(COMP)[::-1]


def readfa(p):
    d, name, buf = {}, None, []
    with open(p) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    d[name] = "".join(buf).upper()
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        d[name] = "".join(buf).upper()
    return d


def old_run(seq, i):
    """codons read in the old frame from i (no skipping) before a stop; capped; None if the contig ends first."""
    n = 0
    while n < CAP:
        cod = seq[i:i + 3]
        if len(cod) < 3:
            return None
        if cod in STOP:
            return n
        n += 1
        i += 3
    return n


def plus1_run(seq, i):
    """codons read from i on a skip walk (every CGG met: drop the C, read on) before a stop; capped."""
    n = 0
    while n < CAP:
        cod = seq[i:i + 3]
        if len(cod) < 3:
            return None
        if cod == "CGG":
            i += 1
            continue
        if cod in STOP:
            return n
        n += 1
        i += 3
    return n


def classify(o, q):
    if o is None or q is None:
        return "unclear"
    if o <= 30 and q >= 60:
        return "skip"
    if o >= 60 and q <= 30:
        return "inframe"
    return "unclear"


def up_runs(seq, p, cds_s0):
    """upstream structure: old-frame codons back to the called CDS start (or stop), and the +1 frame's open
    length backwards from p+1 (codons ending at or before p+1, no skipping), capped at CAP."""
    # +1 frame, backwards: codons at p-2, p-5, ...
    n, i = 0, p - 2
    while n < CAP and i >= 0:
        cod = seq[i:i + 3]
        if cod in STOP:
            break
        n += 1
        i -= 3
    return n


def translate(s):
    return "".join(AA.get(s[i:i + 3], "X") for i in range(0, len(s) - 2, 3))


def main():
    excl = set()
    with open(EXCL) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            excl.add(r["gca"].strip())
    cls = {}
    with open(CLASS) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            cls[r["site_id"]] = r
    by_gen = collections.defaultdict(list)
    with open(SITES) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["set"] in ("A", "C"):
                by_gen[r["genome"]].append(r)
    genomes = sorted(g for g, rows in by_gen.items() if any(x["set"] == "A" for x in rows))
    print("genomes with Set A sites:", len(genomes), "of which on exclusion list:",
          len([g for g in genomes if g in excl]), file=sys.stderr)

    rng = random.Random(7)
    site_out = open(OUT + "sites_up.tsv", "w")
    site_out.write("site_id\tset\tgenome\tcds\tcodon\tup\tup2\tgiven_class\tmy_class\told\tplus1\t"
                   "in_cds\tcodon_idx\tcodons_to_cds_end\tup_plus1_open\tpep_up\n")
    rand_out = open(OUT + "random_up.tsv", "w")
    rand_out.write("genome\tcds\tcodon\tup\tmy_class\told\tplus1\tcodon_idx\tcodons_to_cds_end\n")
    bg = {}  # genome -> {name: Counter}
    cds_usage = {}  # cds id -> Counter of in-frame codons (for CDSs that hold a Set A site)
    checks = collections.Counter()
    for g in genomes:
        fa = readfa(GEN + g + ".fna")
        rcs = {}

        def strand_seq(contig, strand):
            if strand == "+":
                return fa[contig]
            if contig not in rcs:
                rcs[contig] = rc(fa[contig])
            return rcs[contig]

        cds = {}
        with open(GEN + g + ".gff") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 9 or f[2] != "CDS":
                    continue  # ncRNA features ignored
                cid = f[8].split("ID=")[1].split(";")[0]
                s, e, st = int(f[3]), int(f[4]), f[6]
                L = len(fa[f[0]])
                s0 = s - 1 if st == "+" else L - e
                e0 = e if st == "+" else L - s + 1
                cds[cid] = (f[0], st, s0, e0)
        need_usage = {r["cds"] for r in by_gen[g] if r["set"] == "A"}
        b_ggn, b_cgg, b_cgn, b_all = (collections.Counter() for _ in range(4))
        n_internal_stop = 0
        weights = []
        for cid, (contig, st, s0, e0) in cds.items():
            seq = strand_seq(contig, st)
            cods = [seq[i:i + 3] for i in range(s0, e0, 3)]
            if cods and cods[-1] in STOP:
                cods = cods[:-1]
            n_internal_stop += sum(1 for c in cods if c in STOP)
            for k in range(1, len(cods)):
                prev, cur = cods[k - 1], cods[k]
                if prev in STOP:
                    continue
                if cur.startswith("GG"):
                    b_ggn[prev] += 1
                if cur == "CGG":
                    b_cgg[prev] += 1
                elif cur in ("CGA", "CGC", "CGT"):
                    b_cgn[prev] += 1
                b_all[prev] += 1
            if cid in need_usage:
                cds_usage[cid] = collections.Counter(c for c in cods if c not in STOP)
            weights.append((cid, len(cods)))
        checks["internal_stops"] += n_internal_stop
        bg[g] = {"ggn": b_ggn, "cgg": b_cgg, "cgn": b_cgn, "all": b_all}

        for r in by_gen[g]:
            c = int(r["cpos"])
            contig = r["contig"]
            seq = strand_seq(contig, r["strand"])
            L = len(seq)
            p = c if r["strand"] == "+" else L - 1 - c
            cod = seq[p:p + 3]
            checks["codon_ok" if cod == r["codon"] else "codon_BAD"] += 1
            if cod != r["codon"]:
                continue
            up = seq[p - 3:p] if p >= 3 else ""
            up2 = seq[p - 6:p - 3] if p >= 6 else ""
            o = old_run(seq, p + 3)
            q = plus1_run(seq, p + 1)
            mine = classify(o, q)
            gc = ""
            if r["set"] == "A":
                k = cls.get(r["site_id"])
                if k is None:
                    checks["A_no_class"] += 1
                else:
                    gc = k["class"]
                    same = (gc == mine and k["old_frame_codons"] == ("" if o is None else str(o))
                            and k["plus1_frame_codons"] == ("" if q is None else str(q)))
                    checks["A_class_match" if same else "A_class_MISMATCH"] += 1
            ccontig, cst, s0, e0 = cds[r["cds"]]
            in_cds = int(ccontig == contig and cst == r["strand"] and s0 <= p < e0 and (p - s0) % 3 == 0)
            cidx = (p - s0) // 3 if in_cds else -1
            to_end = (e0 - p) // 3 - 1 if in_cds else -1
            upo = up_runs(seq, p, s0)
            pep = translate(seq[p - 63:p - 3]) if (r["set"] == "A" and p >= 63) else ""
            site_out.write("\t".join(map(str, [r["site_id"], r["set"], g, r["cds"], cod, up, up2, gc, mine, o, q,
                                                in_cds, cidx, to_end, upo, pep])) + "\n")

        # random in-frame CDS positions (codon index >= 1, not the terminal codon), classed by the same rule
        tot = sum(w for _, w in weights)
        cum, acc = [], 0
        for cid, w in weights:
            acc += max(w - 1, 0)
            cum.append(acc)
        import bisect
        for _ in range(N_RANDOM if acc > 0 else 0):
            x = rng.randrange(acc)
            j = bisect.bisect_right(cum, x)
            prior = cum[j - 1] if j else 0
            k = x - prior + 1  # codon index 1..ncod-1
            cid = weights[j][0]
            contig, st, s0, e0 = cds[cid]
            seq = strand_seq(contig, st)
            p = s0 + 3 * k
            cod = seq[p:p + 3]
            up = seq[p - 3:p]
            o = old_run(seq, p + 3)
            q = plus1_run(seq, p + 1)
            to_end = (e0 - p) // 3 - 1
            rand_out.write("\t".join(map(str, [g, cid, cod, up, classify(o, q), o, q, k, to_end])) + "\n")
        del fa, rcs
        print(g, "done", file=sys.stderr)
    site_out.close()
    rand_out.close()
    with open(OUT + "background.json", "w") as fh:
        json.dump({g: {k: dict(v) for k, v in d.items()} for g, d in bg.items()}, fh)
    with open(OUT + "cds_usage.json", "w") as fh:
        json.dump({k: dict(v) for k, v in cds_usage.items()}, fh)
    with open(OUT + "checks.json", "w") as fh:
        json.dump({"checks": dict(checks), "excluded": sorted(excl & set(genomes)), "genomes": genomes}, fh, indent=1)
    print(dict(checks), file=sys.stderr)


if __name__ == "__main__":
    main()
