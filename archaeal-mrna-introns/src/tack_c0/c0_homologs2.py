"""(v2, 14:0x: + Sabik 13:51 other-genus copies; inner-helix MFE pairing added) Comparative BHB check of the c0 insert in all 9 Methanosuratincola species reps (insert list: outbox/2026-11-20-tack-c0-intron/
sequences.txt). Altair's motif/fold functions (agents/altair/kills/bhb_fold/bhb_check.py, unchanged). Per species: motif class,
MFE helix pairing at the exact boundaries, and the base pairs at each helix position (covariation = pair changes, still paired)."""
import gzip
src = open("/Users/ryan/lab/agents/altair/kills/bhb_fold/bhb_check.py").read().split("# minus strand")[0]
exec("\n".join(l for l in src.split("\n") if not l.startswith("G = ")))
A = "/Users/ryan/lab/shared/codon-atlas/genomes/"; E = "/Users/ryan/lab/work/tack-rna/ext/msur5/"; X = "/Users/ryan/lab/agents/vega/c0_fold/ext/"  # X: own NCBI downloads
INS = [("GCA_049570475.1", A + "GCA_0495/GCA_049570475.1.fna.gz", "ACGGATCTAGAGTCCGGCCGAA"),
       ("GCA_038801655.1", A + "GCA_0388/GCA_038801655.1.fna.gz", "ACGGATCTAAAGTCCGGCCGAA"),
       ("GCA_030522375.1", A + "GCA_0305/GCA_030522375.1.fna.gz", "ACGGGTGAGAGAATCCGGCCGAA"),
       ("GCA_038801705.1", A + "GCA_0388/GCA_038801705.1.fna.gz", "ACGGATGAGATCCGGCCGAA"),
       ("GCA_029907755.1", E + "GCA_029907755.1.fna.gz", "GCGGATGAAAGAGTCCGGCCGAA"),
       ("GCA_030055545.1", E + "GCA_030055545.1.fna.gz", "ACGGATGAGAAAGTCCGGCCGAA"),
       ("GCA_037478935.1", E + "GCA_037478935.1.fna.gz", "GCGGATAAAAAAAGTCCGGCCGAA"),
       ("GCA_937930465.1", E + "GCA_937930465.1.fna.gz", "GCGGATTCTAAAGTCCGGCCGAA"),
       ("GCA_030055575.1", E + "GCA_030055575.1.fna.gz", "GCGGATCTTAGGGATCCGGCCGAG"),
       ("GCA_038801785.1", A + 'GCA_0388/GCA_038801785.1.fna.gz', "GCGGCCGCCACTGACGGATGGCTGCCGGCCGAG"),  # JAULTD01 33 nt JAULTD010000001.1 -
       ("GCA_964576015.1", X + 'GCA_964576015.1.fna.gz', "ACGACCGGAGGTCGGCCGAA"),  # VGLU01 20 nt CAZNTO010000002.1 -
       ("GCA_014361445.1", X + 'GCA_014361445.1.fna.gz', "GCGACCAAAGGTCGGCCGAG"),  # WYZ-LMO10 20 nt JACIVQ010000030.1 +
       ("GCA_016882195.1", A + 'GCA_0168/GCA_016882195.1.fna.gz', "ACGACCCTCGTATACGCACGTAGGTCGGCCGAG"),  # VGLU01 33 nt VGLU01000005.1 +
       ("GCA_004347955.1", X + 'GCA_004347955.1.fna.gz', "GCGACCGAAAGTCGGCCGAG"),  # WYZ-LMO10 20 nt QNVH01000001.1 +
       ("GCA_020832905.1", X + 'GCA_020832905.1.fna.gz', "ACGACTTTTAGTCGGCCGAA"),  # WYZ-LMO11 20 nt JAJHRG010000001.1 +
       ("GCA_964654895.1", X + 'GCA_964654895.1.fna.gz', "GCGACTAAAAGTCGGCCGAG"),  # WYZ-LMO10 20 nt CAZZJY010000064.1 -
       ("GCA_026414165.1", X + 'GCA_026414165.1.fna.gz', "ACGACCAAAAGTCGGCCGAA"),  # WYZ-LMO10 20 nt JAOAJN010000007.1 -
       ("GCA_964477685.1", X + 'GCA_964477685.1.fna.gz', "GCGGCCTCTGAGCCGGCCGAA"),  # JACTUC01 21 nt CAYZHK010000002.1 -
       ("GCA_029887795.1", A + 'GCA_0298/GCA_029887795.1.fna.gz', "ACGACCAAAAGTCGGCCGAA"),  # WYZ-LMO10 20 nt JARXPE010000004.1 -
       ("GCA_049970435.1", X + 'GCA_049970435.1.fna.gz', "ACGACAGAAAGTCGGCCGAA"),  # JBMFMY01 20 nt JBMFMY010000047.1 -
       ("GCA_004348015.1", A + 'GCA_0043/GCA_004348015.1.fna.gz', "ACGACTTTAAGTCGGCCGAA"),  # WYZ-LMO11 20 nt QNVI01000024.1 -
       ("GCA_964595415.1", X + 'GCA_964595415.1.fna.gz', "GCGACTCTCAGTCGGCCGAA"),  # CAZQQB01 20 nt CAZQQB010000003.1 -
       ("GCA_029887765.1", A + 'GCA_0298/GCA_029887765.1.fna.gz', "ACGGCTATCATAGCCGGCCGAA"),  # JANXER01 22 nt JARXPF010000001.1 -
       ("GCA_964617245.1", X + 'GCA_964617245.1.fna.gz', "ACGACTTTTAGTCGGCCGAA"),  # WYZ-LMO11 20 nt CAZTWF010000085.1 +
       ("GCA_025057955.1", X + 'GCA_025057955.1.fna.gz', "GCGGTCGAGTAGACCGGCCGAA"),  # JANXER01 22 nt JANXER010000010.1 -
       ]
comp = str.maketrans("ACGTN", "TGCAN"); rc = lambda s: s.translate(comp)[::-1]
def contigs(p):
    name, buf = None, []
    for l in gzip.open(p, "rt"):
        if l[0] == ">":
            if name: yield name, "".join(buf)
            name, buf = l[1:].split()[0], []
        else: buf.append(l.strip().upper())
    yield name, "".join(buf)
pairs = {}; INNER = {}
for acc, path, ins in INS:
    hits = []
    for name, s in contigs(path):
        for strand, t in (("+", s), ("-", rc(s))):
            k = t.find("GGTAT" + ins + "TTAG")
            if k < 0: k = t.find("GTAT" + ins + "TTA")
            if k >= 0: hits.append((name, strand, t, t.find(ins, k)))
    if len(hits) != 1:
        print(acc, "hits", len(hits)); continue
    name, strand, t, i = hits[0]; j = i + len(ins)
    # microhomology: can the same removal slide?
    left = 0
    while t[i - 1 - left] == t[j - 1 - left]: left += 1
    right = 0
    while t[i + right] == t[j + right]: right += 1
    m = motif(t, i, j); dg, st, hel, b5, b3, w = fold_check(t, i, j)
    hp = [t[i - 6 + q] + t[j - 3 - q] for q in range(4)]
    import RNA
    a = 40; pt = RNA.ptable(st)
    inner = [t[i + 1 + q] + t[j - 7 - q] for q in range(5)]
    inmfe = sum(pt[a + 1 + q + 1] - 1 == a + (j - i) - 7 - q for q in range(m[2]))
    INNER[acc] = (inner, m[2], inmfe)
    pairs[acc] = hp
    print(f"{acc} {name}{strand} len {len(ins)} slide L{left}/R{right} motif {m} {cls(m)} | MFE dG {dg:.1f} helix {hel}/4 "
          f"bulge5 {b5} bulge3 {b3} | helix pairs {' '.join(hp)} | 5' {t[i-8:i].lower()}|{t[i:i+6]} 3' {t[j-8:j]}|{t[j:j+6].lower()}")
print("\nhelix position: pairs seen across species (covariation if >1 canonical pair at a position)")
for q in range(4):
    c = sorted({p[q] for p in pairs.values()})
    print(q, c, "all paired" if all(x in PAIR for x in (p[q] for p in pairs.values())) else "NOT all paired")

print("\ninner helix: first 5 pairs (insert nt 2-6 with nt L-6..L-10), motif length n, pairs of n present in MFE")
for acc, (inner, n, inm) in INNER.items():
    print(acc, " ".join(inner), "n", n, "MFE", inm, "all canonical(5)", all(x in PAIR for x in inner))
