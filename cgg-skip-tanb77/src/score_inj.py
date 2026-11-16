"""Recall of planted frame breaks: a truth event is found if a break of the same type on the same strand lies within
+-TOL nt. Also reports how many planted events fell inside an aligned locus at all (the ceiling).
usage: python -I src/score_inj.py <truth.tsv> <name.fs.tsv.gz> <name.hits.tsv.gz>
Prints one line per level: all breaks, clean, prog (summarize.flag's definitions)."""
import csv
import gzip
import sys
from bisect import bisect_left
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from summarize import flag  # noqa: E402

TOL = 9


def main():
    truth = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
    fs = list(csv.DictReader(gzip.open(sys.argv[2], "rt"), delimiter="\t"))
    hits = list(csv.DictReader(gzip.open(sys.argv[3], "rt"), delimiter="\t"))
    flag(fs)
    cover = {}
    for h in hits:
        cover.setdefault((h["contig"], h["strand"]), []).append((int(h["gstart"]), int(h["gend"])))
    for level in ("all", "clean", "prog"):
        by = {}
        for r in fs:
            if level == "all" or r[level]:
                by.setdefault((r["contig"], r["strand"], r["type"]), []).append(int(r["gpos"]))
        for v in by.values():
            v.sort()
        found = inside = 0
        for t in truth:
            p = int(t["pos"])
            strands = ["+", "-"] if t["strand"] == "." else [t["strand"]]
            hit = False
            for st in strands:
                v = by.get((t["contig"], st, t["type"]), [])
                i = bisect_left(v, p - TOL)
                if i < len(v) and v[i] <= p + TOL:
                    hit = True
            found += hit
            inside += any(a <= p <= b for st in strands for a, b in cover.get((t["contig"], st), []))
        n = len(truth)
        print(f"{Path(sys.argv[1]).name[:-10]}\t{level}\tplanted {n}\tin_locus {inside}\tfound {found}\t"
              f"recall {found / max(n, 1):.3f}\trecall_in_locus {found / max(inside, 1):.3f}\t"
              f"breaks {sum(len(v) for v in by.values())}")


if __name__ == "__main__":
    main()
