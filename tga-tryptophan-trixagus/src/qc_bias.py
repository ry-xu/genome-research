"""Is GTDB QC biased against recoded genomes? QC-fail fraction among good assemblies (same 'good' definition as
tiers.py) for NCBI genera with a known stop->sense code vs standard-code relatives with similar genomes.
Passed = GTDB R232 bac120 metadata (NCBI genus from ncbi_taxonomy); failed = qc_failed_r232 (genus from NCBI species).
usage: uv run --no-project --with pandas python src/qc_bias.py
"""
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
REF = ROOT.parent.parent / "shared" / "codon-atlas" / "ref"
GROUPS = {
    "TGA-recoded Mollicutes": ["Mycoplasma", "Mycoplasmoides", "Mycoplasmopsis", "Metamycoplasma", "Mesomycoplasma",
                               "Malacoplasma", "Ureaplasma", "Spiroplasma", "Entomoplasma", "Mesoplasma",
                               "Candidatus Hepatoplasma", "Candidatus Bacilliplasma"],
    "standard-code Mollicutes": ["Acholeplasma", "Candidatus Phytoplasma", "Haploplasma", "Alteracholeplasma",
                                 "Paracholeplasma", "Candidatus Phytoplasma"],
    "TGA-recoded insect symbionts": ["Candidatus Hodgkinia", "Candidatus Nasuia", "Candidatus Zinderia",
                                     "Candidatus Stammera"],
    "standard-code insect symbionts": ["Buchnera", "Candidatus Carsonella", "Candidatus Sulcia", "Candidatus Tremblaya",
                                       "Candidatus Blochmannia", "Wigglesworthia", "Candidatus Portiera",
                                       "Candidatus Vidania", "Candidatus Uzinura", "Candidatus Walczuchella",
                                       "Candidatus Baumannia"],
}


def genus_from_species(s):
    w = str(s).replace("s__", "", 1).replace("[", "").replace("]", "").split()
    if not w:
        return ""
    return " ".join(w[:2]) if w[0] == "Candidatus" and len(w) > 1 else w[0]


def good(contigs, n50):
    return ((contigs <= 150) & (n50 >= 20000)) | ((contigs <= 30) & (n50 >= 5000))


def main():
    p = pd.read_csv(REF / "bac120_metadata.tsv.gz", sep="\t", usecols=["accession", "ncbi_taxonomy", "contig_count",
                    "n50_contigs", "ncbi_organism_name"], low_memory=False)
    p["genus"] = p.ncbi_taxonomy.str.extract(r"g__([^;]*)")[0].fillna("")
    alt = p.ncbi_organism_name.map(genus_from_species)
    p.loc[alt.str.startswith("Candidatus"), "genus"] = alt[alt.str.startswith("Candidatus")]
    p["good"] = good(p.contig_count, p.n50_contigs)
    f = pd.read_csv(ROOT / "data" / "qc_failed_r232.tsv.gz", sep="\t", low_memory=False)
    f["genus"] = f["NCBI species"].map(genus_from_species)
    f["good"] = good(f["No. contigs"], f["N50 contigs"])
    rows = []
    for g, names in list(GROUPS.items()) + [("all bacteria", None)]:
        pp = p if names is None else p[p.genus.isin(names)]
        ff = f if names is None else f[f.genus.isin(names)]
        for lab, a, b in (("all", pp, ff), ("good assembly", pp[pp.good], ff[ff.good])):
            n = len(a) + len(b)
            rows.append(dict(group=g, subset=lab, passed=len(a), failed=len(b), fail_frac=round(len(b) / max(n, 1), 3)))
        if names is not None:
            gb = pd.concat([pp[pp.good].genus.value_counts().rename("passed_good"),
                            ff[ff.good].genus.value_counts().rename("failed_good")], axis=1).fillna(0).astype(int)
            print(g, gb.to_dict("index"))
    out = pd.DataFrame(rows)
    out.to_csv(ROOT / "results" / "qc_bias.tsv", sep="\t", index=False)
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
