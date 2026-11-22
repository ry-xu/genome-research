"""Check flagged divergent 23S copies: identity to the genome's best copy, nearest census 23S in other genomes, host QC.

usage: python -I src/div23_check.py FLAGGED.tsv SCANDIR CENSUS_DB CENSUS_TAX OUT.tsv
  FLAGGED.tsv  output of div23_post.py
  CENSUS_DB    blastn db of census 23S, one copy per genome, ids = GTDB accession (GB_/RS_ prefix)
  CENSUS_TAX   acc<TAB>completeness<TAB>contamination<TAB>gtdb taxonomy (GB_/RS_ prefix), bacteria + archaea
dc-megablast (megablast misses ~75%-identity pairs). Identities are summed over all HSPs of a subject (matches / aligned length; aligned length also reported), so a
divergent copy is not judged on its single best local block. Writes one row per flagged copy; the
best copy's nearest other genome is reported too, to tell which copy is the foreign one.
"""
import collections
import glob
import os
import subprocess
import sys
import tempfile


def blast(q, db_args):
    out = subprocess.run(["blastn", "-task", "dc-megablast", "-query", q, *db_args, "-outfmt", "6 sseqid pident length nident qstart qend",
                          "-max_target_seqs", "50", "-max_hsps", "20", "-evalue", "1e-10"],
                         capture_output=True, text=True, check=True).stdout
    agg = collections.defaultdict(lambda: [0, 0])
    for line in out.splitlines():
        s, _, ln, ni, _, _ = line.split("\t")
        agg[s][0] += int(ni)
        agg[s][1] += int(ln)
    return sorted(((v[0] / v[1], v[1], s) for s, v in agg.items() if v[1] >= 300), key=lambda x: (-x[1] * x[0]))


def main(flagged, scandir, db, taxp, outp):
    tax = {}
    for line in open(taxp):
        a, comp, cont, t = line.rstrip("\n").split("\t")
        tax[a[3:]] = (comp, cont, t)
    seqs = {}
    for p in glob.glob(f"{scandir}/copies.*.fna"):
        name = None
        for line in open(p):
            if line.startswith(">"):
                name = line[1:].strip()
            else:
                seqs[name] = line.strip()
    tmp = tempfile.mkdtemp()
    with open(outp, "w") as fo:
        fo.write("acc\tcopy\tratio\tid_to_best\taln_to_best\tnearest_other\tid_other\taln_other\tother_tax\t"
                 "best_nearest\tbest_id\tbest_aln\tbest_tax\t"
                 "host_comp\thost_cont\thost_tax\n")
        for line in open(flagged):
            if line.startswith("#"):
                continue
            acc, ctg, lo, hi, strand, hspan, dens, bdens, ratio, best = line.rstrip("\n").split("\t")
            key = f"{acc}|{ctg}|{lo}|{hi}|{strand}|{hspan}"
            bctg, brest = best.split(":")
            blo, bhi = brest[:-1].split("-")
            bkey = f"{acc}|{bctg}|{blo}|{bhi}|{brest[-1]}|"
            bname = next((k for k in seqs if k.startswith(bkey)), None)
            if key not in seqs or bname is None:
                print(f"{key}: sequence missing", file=sys.stderr)
                continue
            q, s = f"{tmp}/q.fa", f"{tmp}/s.fa"
            open(q, "w").write(f">q\n{seqs[key]}\n")
            open(s, "w").write(f">b\n{seqs[bname]}\n")
            tb = blast(q, ["-subject", s])
            idb, alb = (f"{tb[0][0]:.3f}", tb[0][1]) if tb else ("na", 0)
            near = []
            for path in (q, s):
                others = [h for h in blast(path, ["-db", db]) if h[2][3:] != acc]
                if others:
                    o = max(others, key=lambda h: h[0] if h[1] >= 1000 else h[0] - 1)
                    near += [o[2], f"{o[0]:.3f}", o[1], tax.get(o[2][3:], ("", "", ""))[2]]
                else:
                    near += ["none", "na", 0, ""]
            hc, hx, ht = tax.get(acc, ("", "", ""))
            fo.write("\t".join(map(str, [acc, key, ratio, idb, alb, *near, hc, hx, ht])) + "\n")
            fo.flush()


if __name__ == "__main__":
    main(*sys.argv[1:6])
