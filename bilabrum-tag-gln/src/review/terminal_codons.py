"""Terminal codon after DIAMOND hits with send == slen, from new-stops-euk stored ends tables (deduplicated by locus).
usage: python3 -I terminal_codons.py ENDS.tsv [...]"""
import sys
from collections import Counter
for f in sys.argv[1:]:
    seen, c = set(), Counter()
    for line in list(open(f))[1:]:
        contig, sid, slen, send, pid, pre, post = line.rstrip("\n").split("\t")[:7]
        key = (contig.split(":")[0], post[:30])
        if key in seen or int(send) != int(slen):
            continue
        seen.add(key)
        c[post[:3] if post[:3] in ("TAA", "TAG", "TGA") else "sense"] += 1
    n = c["TAA"] + c["TAG"] + c["TGA"]
    print(f"{f.split('/')[-1]}\tstops {n}\tTGA {c['TGA']}\tTAA {c['TAA']}\tTAG {c['TAG']}\tTAG share {c['TAG'] / max(n, 1):.3f}\tnext codon sense {c['sense']}")
