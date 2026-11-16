"""Amendment 5 of the RNA test (results/rna_cgg/README.md; post hoc, the fact-checker's definition adopted unchanged):
split Set A sites by what the two frames do after the CGG, then count keeps/lacks per class in both pipelines.
  skip-like:     the old frame meets a stop within 30 codons after the CGG, and the +1 frame (from GGN) runs at least
                 60 codons without a stop;
  in-frame-like: the reverse (old frame runs >= 60 codons, +1 frame stops within 30);
  unclear:       anything else (including a contig end before either bound is settled).
Stops are TAA/TAG/TGA (table 11). Written independently of the fact-checker's script; its conventions, which are the
definition adopted, are the default: "within 30" counts up to 30 codons before the stop (<= 30), and the +1 frame
skips again at any later CGG, as the skip walk does. "strict" (my first reading): < 30, and the +1 frame read plainly.
usage: python -I src/rna_skiplike.py <out_dir> [strict]   (run with the atlas Python, from the project root)
"""
import csv
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from proteomics_db import rc, read_fasta  # noqa: E402

# RNA_MDIR / RNA_RES override the pig inputs (chicken replication; main() reads the pig follow-up files)
MGN = Path(os.environ.get("RNA_MDIR", ROOT / "data" / "proteomics" / "mgnify_pig"))
RNA = Path(os.environ.get("RNA_RES", ROOT / "results" / "rna_cgg"))
STOPS = {"TAA", "TAG", "TGA"}
NEAR, FAR = 30, 60
STRICT = False


def run_len(s, start, cap, skip=False):
    """codons from s[start] before the first stop (cap if none within cap codons); None if the sequence ends first.
    skip: a CGG met on the way loses its C (the walk moves on by one base, no codon counted)"""
    k, i = 0, start
    while k < cap:
        cod = s[i:i + 3]
        if len(cod) < 3:
            return None
        if cod in STOPS:
            return k
        if skip and cod == "CGG":
            i += 1
            continue
        k += 1
        i += 3
    return cap


def classify(seqs, contig, cpos, strand):
    """cpos: 0-based plus-strand position of the C (a G on the plus strand for '-' genes)"""
    x = seqs[contig]
    s = x[cpos:] if strand == "+" else rc(x[:cpos + 1])
    assert s[:3] == "CGG", (contig, cpos, strand, s[:3])
    old = run_len(s, 3, FAR)    # codons after the CGG in the old frame
    new = run_len(s, 1, FAR, skip=not STRICT)    # codons from GGN in the +1 frame
    near = (lambda n: n < NEAR) if STRICT else (lambda n: n <= NEAR)
    if old is not None and new is not None:
        if near(old) and new >= FAR:
            return "skip", old, new
        if old >= FAR and near(new):
            return "inframe", old, new
    return "unclear", old, new


def obs_table(path, fmt):
    """per-site keeps/lacks of Set A (and C, kept for decide's power rule); key (contig, 0-based cpos)"""
    out = []
    for r in csv.DictReader(open(path), delimiter="\t"):
        if r["set"] not in ("A", "C"):
            continue
        if fmt == "mine":
            out.append((r["site_id"], r["set"], r["genome"], None, None, int(r["keeps"]), int(r["lacks"])))
        else:  # re-check ids: genome:contig:pos(1-based):strand
            g, c, p, st = r["site_id"].split(":")
            out.append((r["site_id"], r["set"], g, (c, int(p) - 1), st, int(r["keeps"]), int(r["lacks"])))
    return out


def main(outdir):
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    sites = {r["site_id"]: r for r in csv.DictReader(open(RNA / "sites.tsv"), delimiter="\t")}
    seqs = {}
    for g in sorted({r["genome"] for r in sites.values()}):
        seqs.update(read_fasta(MGN / f"{g}.fna"))
    cls = {}  # (contig, cpos) -> class
    with open(out / "site_class.tsv", "w") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["site_id", "set", "genome", "contig", "cpos", "strand", "cds", "class", "old_frame_codons",
                    "plus1_frame_codons"])
        for sid, r in sites.items():
            if r["set"] not in ("A", "B"):
                continue
            c, p = r["contig"], int(r["cpos"])
            k, old, new = classify(seqs, c, p, r["strand"])
            cls[(c, p)] = k
            w.writerow([sid, r["set"], r["genome"], c, p, r["strand"], r["cds"], k, old, new])
    tab = Counter((r["set"], cls[(r["contig"], int(r["cpos"]))]) for r in sites.values() if r["set"] in ("A", "B"))
    print("site classes (all sites):", dict(sorted(tab.items())))
    mine = {r["site_id"]: (r["contig"], int(r["cpos"])) for r in sites.values()}
    runs = [("S8", "mine", RNA / "followS" / "S8.obs.tsv"), ("S9", "mine", RNA / "followS" / "S9.obs.tsv"),
            ("S8", "recheck", RNA / "recheck_S" / "per_site_counts.tsv"),
            ("S9", "recheck", RNA / "recheck_S" / "per_site_counts_plus_SRR37726778.tsv")]
    summary = []
    for tag, fmt, path in runs:
        rows = obs_table(path, fmt)
        per = defaultdict(lambda: [0, 0, 0, Counter()])  # obs, lacks, sites, genome -> obs
        unknown = 0
        files = {k: open(out / f"{tag}_{fmt}_{k}.obs.tsv", "w") for k in ("skip", "inframe", "unclear")}
        for f in files.values():
            f.write("site_id\tset\tkeeps\tlacks\n")
        for sid, st, g, key, strand, k, l in rows:
            if st == "C":
                for f in files.values():
                    f.write(f"{sid}\tC\t{k}\t{l}\n")
                continue
            key = mine[sid] if fmt == "mine" else key
            if key in cls:
                c = cls[key]
            else:  # a re-check site absent from my Set A: classify it the same way
                unknown += 1
                c = classify(seqs, key[0], key[1], strand)[0]
            files[c].write(f"{sid}\tA\t{k}\t{l}\n")
            if k + l:
                p = per[c]
                p[0] += k + l
                p[1] += l
                p[2] += 1
                p[3][g] += k + l
        for f in files.values():
            f.close()
        for c in ("skip", "inframe", "unclear"):
            o, l, s, gs = per[c]
            top = max(gs.values()) if gs else 0
            summary.append([tag, fmt, c, o, l, s, len(gs), top])
        print(tag, fmt, "re-check sites not in my Set A:", unknown)
    with open(out / "summary.tsv", "w") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["followup", "pipeline", "class", "observations", "lacks", "sites_observed", "genomes",
                    "top_genome_obs"])
        w.writerows(summary)
    for r in summary:
        print("\t".join(map(str, r)))


if __name__ == "__main__":
    STRICT = sys.argv[2:3] == ["strict"]
    main(sys.argv[1])
