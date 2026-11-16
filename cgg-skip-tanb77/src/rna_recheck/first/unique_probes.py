"""Re-check, step 2: drop sites whose keep or lack probe occurs at more than one place in the 24 genomes.
  .venv/bin/python -I unique_probes.py
Every keep and lack probe of the Set A and Set C candidates is searched (Aho-Corasick, exact) in all contigs of the 24
genomes, on both strands. A place is (genome, contig, plus-strand start, length), so a probe that is its own reverse
complement counts once at its own place. A site is kept only if its keep probe occurs at exactly 1 place (its own)
and its lack probe at 0 places, and neither probe holds a non-ACGT base.
Amendment 3: a site is dropped (status ncrna, checked before uniqueness) if its observation window
[run_start0 - 20, run_end0 + 20) (0-based half-open, plus strand) overlaps any GFF feature of type ncRNA (rRNA, tRNA,
tmRNA, RNase P, riboswitches...) in the catalogue GFFs; GFF start..end (1-based inclusive) is [start - 1, end).
Output: sites_recheck.tsv (all candidates with status; only status == ok are counted)
"""
import csv
import glob
from collections import Counter, defaultdict

import ahocorasick

R = "results/rna_cgg/recheck/"
MDIR = "data/proteomics/mgnify_pig/"
COMP = str.maketrans("ACGTN", "TGCAN")


def rc(s):
    return s.translate(COMP)[::-1]


def fasta(path):
    name, buf = None, []
    for line in open(path):
        if line.startswith(">"):
            if name:
                yield name, "".join(buf).upper()
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name:
        yield name, "".join(buf).upper()


def ncrna_features():
    """(genome, contig) -> list of 0-based half-open [start - 1, end) of every GFF feature of type ncRNA"""
    feat, n = defaultdict(list), 0
    for f in sorted(glob.glob(MDIR + "*.gff")):
        g = f.rsplit("/", 1)[1][:-4]
        for line in open(f):
            if line.startswith("##FASTA"):
                break
            x = line.rstrip("\n").split("\t")
            if line[0] == "#" or len(x) < 9 or x[2] != "ncRNA":
                continue
            feat[(g, x[0])].append((int(x[3]) - 1, int(x[4])))
            n += 1
    print("ncRNA features", n)
    return feat


def main():
    nc = ncrna_features()
    cand = [r for r in csv.DictReader(open(R + "candidates.tsv"), delimiter="\t") if r["set"] in ("A", "C")]
    A = ahocorasick.Automaton()
    for r in cand:
        for p in (r["keep_probe"], r["lack_probe"]):
            if not A.exists(p):
                A.add_word(p, p)
    A.make_automaton()
    places = defaultdict(set)
    fnas = sorted(glob.glob(MDIR + "*.fna"))
    assert len(fnas) == 24, len(fnas)
    for f in fnas:
        g = f.rsplit("/", 1)[1][:-4]
        for c, x in fasta(f):
            n = len(x)
            for e, p in A.iter(x):
                places[p].add((g, c, e - len(p) + 1, len(p)))
            for e, p in A.iter(rc(x)):
                places[p].add((g, c, n - 1 - e, len(p)))
    out, why = [], Counter()
    for r in cand:
        k, l = r["keep_probe"], r["lack_probe"]
        own = (r["genome"], r["contig"], int(r["run_start0"]) - 20, len(k))
        ws, we = int(r["run_start0"]) - 20, int(r["run_end0"]) + 20
        r["ncrna"] = int(any(fs < we and ws < fe for fs, fe in nc[(r["genome"], r["contig"])]))
        if set(k) - set("ACGT"):
            st = "non_acgt"
        elif r["ncrna"]:
            st = "ncrna"
        elif own not in places[k]:
            st = "own_place_missing"  # should never happen
        elif len(places[k]) > 1:
            st = "keep_multi"
        elif places[l]:
            st = "lack_present"
        else:
            st = "ok"
        why[(r["set"], st)] += 1
        r["keep_places"], r["lack_places"], r["status"] = len(places[k]), len(places[l]), st
        r["site_id"] = f'{r["genome"]}:{r["contig"]}:{r["cpos1"]}:{r["strand"]}'
        out.append(r)
    cols = ["site_id", "set", "genome", "contig", "cpos1", "strand", "codon", "run_start0", "run_end0", "run_len", "cds",
            "walk_has_cgg", "best_walk", "ncrna", "keep_places", "lack_places", "status", "keep_probe", "lack_probe"]
    with open(R + "sites_recheck.tsv", "w") as fh:
        w = csv.DictWriter(fh, cols, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(out)
    for k in sorted(why):
        print(*k, why[k])


if __name__ == "__main__":
    main()
