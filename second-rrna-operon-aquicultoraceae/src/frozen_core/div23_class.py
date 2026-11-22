"""Classify flagged divergent 23S copies from div23_check.py output.

usage: python3 -I src/div23_class.py CHECK.tsv SCANDIR > CLASS.tsv   (SCANDIR holds copies.*.fna)
Classes (id_other = identity of the divergent copy to its nearest census 23S in another genome):
  foreign   id_other >= 0.97 and the nearest genome is in another genus: the copy belongs to another organism (contamination/co-binning)
  shared    id_other >= 0.97 and the nearest genome is in the host's genus: the copy recurs in relatives (best leads)
  singleton id_other < 0.97: seen once; a MAG singleton needs the co-assembly checks (results/lsu/bac/cca/arist/NOTES.md)
  ref-foreign  the genome's best (reference) copy is >= 0.97 to a genome in another genus: the reference itself is a
            contaminant, so id_to_best and the flag are meaningless (Glycomyces GCF_003265355 and Actinoplanes GCF_051386375
            in the 22k partial scan: a Bacillus-like 23S taken as the reference). Checked before the other classes.
  rrna-pile more than 10 23S loci (any length) in the genome: rRNA contigs from many organisms binned together (the
            "99.95% complete" Collinsella MAG GCA_938023335 has 27). Checked first; the n_loci column is added for all rows.
Genus comes from the GTDB taxonomy strings in other_tax and host_tax.
"""
import collections
import glob
import sys


def genus(tax):
    return next((x for x in tax.split(";") if x.startswith("g__")), "?")


def main(path, scandir):
    nloci = collections.Counter()
    for f in glob.glob(scandir + "/copies.*.fna"):
        for line in open(f):
            if line[0] == ">":
                nloci[line[1:].split("|")[0]] += 1
    rows = [l.rstrip("\n").split("\t") for l in open(path)]
    hdr, rows = rows[0], rows[1:]
    ix = {h: i for i, h in enumerate(hdr)}
    print("\t".join(["class", "n_loci"] + hdr))
    for r in rows:
        try:
            idn = float(r[ix["id_other"]])
        except ValueError:
            idn = 0.0
        same = genus(r[ix["other_tax"]]) == genus(r[ix["host_tax"]])
        try:
            bid = float(r[ix["best_id"]])
        except ValueError:
            bid = 0.0
        ref_foreign = bid >= 0.97 and genus(r[ix["best_tax"]]) != genus(r[ix["host_tax"]])
        n = nloci[r[ix["acc"]]]
        c = "rrna-pile" if n > 10 else "ref-foreign" if ref_foreign else "singleton" if idn < 0.97 else "shared" if same else "foreign"
        print("\t".join([c, str(n)] + r))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
