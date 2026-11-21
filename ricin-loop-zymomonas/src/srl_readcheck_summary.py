"""Summary of read-level SRL checks (tally v2): one row per carrier from results/lsu/bac/reads/runs/*/tally.tsv.
usage: python -I src/srl_readcheck_summary.py OUT.tsv
Identity is outside the 12 site bases (src/srl_readcheck.py v2). hi_* = >= 0.98: reads carrying the assembled dodecamer
(reads / molecules = distinct both-mate placements, see srl_readcheck.py), reads with another full dodecamer, reads with an indel
near the site. mid_assembled / mid_other = reads at 0.90-0.98 with the assembled / the most common other dodecamers
(relatives; reads dropped by the span and coverage filters are not in these). dropped = touching the site without
spanning it +/- 10 nt."""
import collections, csv, glob, os, re, sys
W = os.environ.get("FROZEN_CORE", "frozen-core")  # lab working tree of the census
pick = {r["acc"]: r for r in csv.DictReader(open(f"{W}/results/lsu/bac/reads/pick.tsv"), delimiter="\t")}
out = []
for t in sorted(glob.glob(f"{W}/results/lsu/bac/reads/runs/*/tally.tsv")):
    d = os.path.dirname(t); acc = os.path.basename(d)
    lines = open(t).readlines()
    if "molecules_by_run" not in lines[1]:
        out.append([acc, pick[acc]["set"], pick[acc]["phylum"][3:], pick[acc]["order"][3:], "", "", "PENDING (old tally format)"] + [""] * 8)
        continue
    asm = lines[0].split()[3]
    dropped = re.search(r"not spanning \+/-10: (\d+)", lines[0]).group(1)
    runs = sorted(os.path.basename(p)[:-9] for p in glob.glob(f"{d}/*.paf.done"))
    hi_agree = hi_place = hi_other = hi_indel = mid_asm = 0; byrun = ""; oth = []; mid = collections.Counter()
    for r in csv.DictReader(lines[1:], delimiter="\t"):
        n = int(r["reads"])
        if r["identity_outside_site"] == "≥0.98":
            if r["dodecamer"] == asm:
                hi_agree, hi_place = n, int(r["molecules"]); byrun = r["molecules_by_run"]
            elif r["dodecamer"] == "indel_or_gap":
                hi_indel += n
            else:
                hi_other += n; oth.append(f"{r['dodecamer']}:{n}")
        elif r["identity_outside_site"] == "0.90-0.98" and r["dodecamer"] == asm:
            mid_asm += n
        elif r["identity_outside_site"] == "0.90-0.98" and r["dodecamer"] != "indel_or_gap":
            mid[r["dodecamer"]] += n
    out.append([acc, pick[acc]["set"], pick[acc]["phylum"][3:], pick[acc]["order"][3:], asm, ",".join(runs),
                hi_agree, hi_place, byrun, hi_other, hi_indel, ";".join(oth), mid_asm, ";".join(f"{k}:{v}" for k, v in mid.most_common(2)), dropped])
with open(sys.argv[1], "w") as f:
    f.write("acc\tset\tphylum\torder\tassembled\truns\thi_agree_reads\thi_agree_molecules\thi_agree_molecules_by_run\thi_other\thi_indel\thi_other_detail\tmid_assembled\tmid_other\tdropped_partial\n")
    for r in out:
        f.write("\t".join(map(str, r)) + "\n")
for r in out:
    print(*r, sep="\t")
