"""Recount the README's numbers from the tables in data/.
usage (from the folder root): uv run python -I src/check_numbers.py

- in-frame stops in the 76-protein tblastn panel (data/review/bil.sites.tsv, caf.sites.tsv; the reviewer's rules in
  src/review/analyse_sites.py: aligned = opposite a query residue; clean = aligned, >= 10 columns from the HSP ends,
  no query gap within 6 columns, >= 8 of 25 identical in the window)
- gene ends (data/ends/*.ends.tsv, deduplicated by locus as in src/review/terminal_codons.py)
- RNA: codon at the clean TAG sites in the three pooled 1KP runs (data/rna/ERR*.sites.tsv, as src/rna/site_table.py),
  and between-strain C->T rates over the panel (data/rna/ALL3.codons.tsv.gz, as src/rna/codon_pairs.py)
"""
import gzip
from collections import Counter, defaultdict

STOPS = ("TAA", "TAG", "TGA")


def rows(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            yield dict(zip(head, line.rstrip("\n").split("\t")))


def stop_sites(prefix):
    aligned, clean, clean_res = Counter(), Counter(), Counter()
    out = []
    for r in rows(f"data/review/{prefix}.sites.tsv"):
        if r["qres"] in ("-", "None"):
            continue
        aligned[r["codon"]] += 1
        ok = int(r["edge_dist"]) >= 10 and r["qgaps6"] == "0" and int(r["win_ident"]) >= 8
        if ok:
            clean[r["codon"]] += 1
            if r["codon"] == "TAG":
                clean_res[r["qres"]] += 1
        out.append((r, ok))
    return aligned, clean, clean_res, out


def gene_ends(path):
    seen, c = set(), Counter()
    for r in list(open(path))[1:]:
        contig, sid, slen, send, pid, pre, post = r.rstrip("\n").split("\t")[:7]
        key = (contig.split(":")[0], post[:30])
        if key in seen or int(send) != int(slen):
            continue
        seen.add(key)
        c[post[:3] if post[:3] in STOPS else "sense"] += 1
    return c


def main():
    al, cl, res, bil = stop_sites("bil")
    print(f"Bilabrum, aligned in-frame stops: TAG {al['TAG']}, TAA {al['TAA']}, TGA {al['TGA']}"
          "   (README: 166, 1, 4)")
    print(f"  clean TAG sites {cl['TAG']}; human residue there: " + ", ".join(f"{k} {v}" for k, v in res.most_common(4))
          + "   (README: 91; Q 44, K 11)")
    al_c, _, _, _ = stop_sites("caf")
    print(f"Cafeteria, aligned in-frame TAG: {al_c['TAG']}   (README: none)")

    names = {"GCA_965178125.1": "Bilabrum", "GCA_008330645.1": "Cafeteria", "GCA_965644315.1": "Caecitellus",
             "GCA_057816405.1": "Halocafeteria", "GCA_964026335.1": "uncultured Cafeteriaceae"}
    for acc, name in names.items():
        c = gene_ends(f"data/ends/{acc}.ends.tsv")
        n = c["TAA"] + c["TAG"] + c["TGA"]
        print(f"gene ends, {name} {acc}: TAG {c['TAG']} of {n} (TGA {c['TGA']}, TAA {c['TAA']})")
    print("   (README: Bilabrum 0 of 149, TGA 136, TAA 13; Cafeteria 11 of 145; the other three 35-39 of 149-172)")

    # RNA at clean TAG sites, the three runs pooled
    clean_keys = {(r["contig"], r["pos1"], r["strand"]) for r, ok in bil if ok and r["codon"] == "TAG"}
    pooled = defaultdict(Counter)
    for run in ("ERR2041020", "ERR2041128", "ERR3487419"):
        for r in rows(f"data/rna/{run}.sites.tsv"):
            key = (r["contig"], r["pos1"], r["strand"])
            if key in clean_keys:
                pooled[key]["TAG"] += int(r["TAG"])
                pooled[key]["CAG"] += int(r["CAG"])
                for x in filter(None, r["other"].split(",")):
                    cod, v = x.split(":")
                    if cod != "indel":
                        pooled[key][cod] += int(v)
    maj = Counter(c.most_common(1)[0][0] for c in pooled.values() if sum(c.values()) >= 2)
    print(f"RNA, clean TAG sites covered by >= 2 reads: {sum(maj.values())}; majority read codon "
          + ", ".join(f"{k} {v}" for k, v in maj.most_common()) + "   (README: 60; TAG 41, CAG 18)")

    # between-strain changes over the panel
    gcount, pairs = Counter(), Counter()
    for r in rows("data/rna/ALL3.codons.tsv.gz"):
        if r["call"]:
            gcount[r["genome"]] += 1
            pairs[(r["genome"], r["call"])] += 1

    def rate(g, t):
        n = sum(gcount[x] for x in g)
        k = sum(pairs[(x, y)] for x in g for y in t)
        return k, n, 1000 * k / n

    pro = ("CCA", "CCC", "CCG", "CCT")
    k1, n1, r1 = rate(("CAG",), ("TAG",))
    k2, n2, r2 = rate(("CTG",), ("TTG",))
    ks = [rate(("CAC", "CAT"), ("TAC", "TAT")), rate(("CGG",), ("TGG",)), rate(pro, tuple("T" + p[1:] for p in pro))]
    k3, n3 = sum(x[0] for x in ks), sum(x[1] for x in ks)
    print(f"genome CAG, 1KP strain TAG: {k1} of {n1} ({r1:.1f} per 1,000)   (README: 12.3)")
    print(f"genome CTG, 1KP strain TTG (silent): {k2} of {n2} ({r2:.1f} per 1,000)   (README: 14.4)")
    print(f"His->Tyr, Arg->Trp, Pro->Ser together: {k3} of {n3} ({1000 * k3 / n3:.1f} per 1,000)   (README: 0.8)")


if __name__ == "__main__":
    main()
