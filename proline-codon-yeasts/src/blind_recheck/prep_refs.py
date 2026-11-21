"""Select canonical core-gene CDSs from reference fasta_cds_na, translate with table 4.

usage: prep_refs.py refs_cds_na.fna outdir ACC1,ACC2,...
Writes per gene: ref_<gene>.t4.faa (table-4 translation, stops '*'),
ref_<gene>.mask.faa (contested codons -> X), ref_<gene>.cds.fna, and refs_selected.tsv.
"""
import os
import re
import sys
import json
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import read_fasta, write_fasta, codons_of, translate_codons, T4, CONTESTED  # noqa: E402

GENES = ["cox1", "cox2", "cox3", "cob", "atp6", "atp8", "atp9",
         "nad1", "nad2", "nad3", "nad4", "nad4L", "nad5", "nad6"]
ALIASES = {
    "cox1": ["cox1", "coxi", "co1", "coi"], "cox2": ["cox2", "coxii", "co2"], "cox3": ["cox3", "coxiii", "co3"],
    "cob": ["cob", "cytb", "cyb"], "atp6": ["atp6"], "atp8": ["atp8"], "atp9": ["atp9", "oli1"],
    "nad1": ["nad1", "nd1"], "nad2": ["nad2", "nd2"], "nad3": ["nad3", "nd3"], "nad4": ["nad4", "nd4"],
    "nad4L": ["nad4l", "nd4l"], "nad5": ["nad5", "nd5"], "nad6": ["nad6", "nd6"],
}
PROT_ALIASES = {  # fall back on protein name when gene tag missing
    "nad4": "nadh dehydrogenase subunit 4$", "nad5": "nadh dehydrogenase subunit 5$",
    "nad3": "nadh dehydrogenase subunit 3$", "nad1": "nadh dehydrogenase subunit 1$",
    "nad2": "nadh dehydrogenase subunit 2$", "nad6": "nadh dehydrogenase subunit 6$",
    "nad4L": "nadh dehydrogenase subunit 4l$",
}
BAD = re.compile(r"endonuclease|maturase|intron|reverse transcriptase|hypothetical|orf|giy|laglidadg", re.I)


def norm_gene(g):
    g = g.lower()
    for k, al in ALIASES.items():
        if g in al:
            return k
    return None


def main():
    src, outdir, accs = sys.argv[1], sys.argv[2], sys.argv[3].split(",")
    os.makedirs(outdir, exist_ok=True)
    recs = read_fasta(src)
    cand = defaultdict(list)
    for hdr, seq in recs:
        acc = re.match(r"lcl\|(.+?)_cds_", hdr).group(1).split(".")[0]
        if acc not in accs:
            continue
        gene = re.search(r"\[gene=([^\]]+)\]", hdr)
        prot = re.search(r"\[protein=([^\]]+)\]", hdr)
        prot = prot.group(1) if prot else ""
        g = norm_gene(gene.group(1)) if gene else None
        if g is None and not gene:
            for k, pat in PROT_ALIASES.items():
                if re.search(pat, prot.lower()):
                    g = k
        if g is None or BAD.search(prot) or "partial=" in hdr:
            continue
        cand[(acc, g)].append((hdr, seq))
    sel = []
    per_gene = defaultdict(list)
    info = {}
    for (acc, g), lst in sorted(cand.items()):
        # if several, take the longest non-partial
        hdr, seq = max(lst, key=lambda x: len(x[1]))
        cods = codons_of(seq)
        t4 = translate_codons(cods, T4)
        if t4.endswith("*"):
            t4, cods = t4[:-1], cods[:-1]
        nstop = t4.count("*")
        mask = translate_codons(cods, T4, mask=CONTESTED | {"TAA", "TAG"})
        name = f"{acc}|{g}"
        per_gene[g].append((name, t4, mask, "".join(cods)))
        info[name] = {"len": len(t4), "internal_stops": nstop, "n_alts": len(lst),
                      "loc": (re.search(r"\[location=([^\]]+)\]", hdr) or [None, ""])[1][:80]}
        sel.append((acc, g, len(t4), nstop, len(lst)))
    for g in GENES:
        lst = per_gene.get(g, [])
        write_fasta(os.path.join(outdir, f"ref_{g}.t4.faa"), [(n, t) for n, t, m, c in lst])
        write_fasta(os.path.join(outdir, f"ref_{g}.mask.faa"), [(n, m) for n, t, m, c in lst])
        write_fasta(os.path.join(outdir, f"ref_{g}.cds.fna"), [(n, c) for n, t, m, c in lst])
    with open(os.path.join(outdir, "refs_selected.tsv"), "w") as fh:
        fh.write("acc\tgene\tlen_aa\tinternal_stops_t4\tn_candidates\n")
        for r in sel:
            fh.write("\t".join(map(str, r)) + "\n")
    json.dump(info, open(os.path.join(outdir, "refs_info.json"), "w"), indent=1)
    # matrix summary
    accs_present = sorted({a for a, *_ in sel})
    print("acc\t" + "\t".join(GENES))
    have = {(a, g): (L, s) for a, g, L, s, n in sel}
    for a in accs_present:
        print(a + "\t" + "\t".join(f"{have[(a, g)][0]}{'!' + str(have[(a, g)][1]) if have.get((a, g), (0, 0))[1] else ''}" if (a, g) in have else "-" for g in GENES))


if __name__ == "__main__":
    main()
