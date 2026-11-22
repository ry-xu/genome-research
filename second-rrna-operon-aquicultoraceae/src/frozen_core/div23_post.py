"""Flag 23S copies that score far below the best copy in the same genome (candidate divergent paralogs).

usage: python -I src/div23_post.py SCANDIR [MIN_SPAN] [MAX_RATIO]
Reads SCANDIR/hits.*.tsv (from div23_scan.py). Density = nhmmer score / model span covered. A copy is flagged when it
spans >= MIN_SPAN (default 2000) model positions and its density is <= MAX_RATIO (default 0.85) x the density of the
best copy of span >= MIN_SPAN in the same genome. Prints acc, contig, lo, hi, strand, hmm span, density, best density,
ratio, and the best copy's coordinates; then a summary line.
"""
import collections
import glob
import sys


def main(scandir, min_span=2000, max_ratio=0.85):
    by_acc = collections.defaultdict(list)
    for p in glob.glob(f"{scandir}/hits.*.tsv"):
        for line in open(p):
            acc, ctg, lo, hi, strand, hf, ht, score = line.rstrip("\n").split("\t")
            span = int(ht) - int(hf) + 1
            if span >= min_span:
                by_acc[acc].append((float(score) / span, ctg, int(lo), int(hi), strand, int(hf), int(ht)))
    n_multi = n_flag = 0
    for acc in sorted(by_acc):
        cps = sorted(by_acc[acc], reverse=True)
        if len(cps) < 2:
            continue
        n_multi += 1
        best = cps[0]
        for c in cps[1:]:
            if c[0] <= max_ratio * best[0]:
                n_flag += 1
                print(acc, c[1], c[2], c[3], c[4], f"{c[5]}-{c[6]}", f"{c[0]:.3f}", f"{best[0]:.3f}",
                      f"{c[0] / best[0]:.3f}", f"{best[1]}:{best[2]}-{best[3]}{best[4]}", sep="\t")
    print(f"# genomes with >=2 copies spanning >={min_span}: {n_multi}; flagged copies: {n_flag}; "
          f"genomes with any copy >={min_span}: {len(by_acc)}", file=sys.stderr)


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 2000,
         float(sys.argv[3]) if len(sys.argv) > 3 else 0.85)
