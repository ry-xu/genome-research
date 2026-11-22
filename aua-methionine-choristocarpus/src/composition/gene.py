"""Per gene: Choristocarpus ATG, ATA, ATT+ATC counts vs Phaeophyceae relatives' median Met and Ile counts per same gene
(relatives: ATG = Met, ATA+ATT+ATC = Ile; all scaled to Choristocarpus gene length). Run after comp.py logic."""
import re, gzip, glob, sys, collections, statistics as st
import polars as pl
exec(open("comp.py").read().split("out = []")[0])  # reuse taxonomy + cds + parse
phae = {a for a, t in at.items() if grp[t][0] == "Phaeophyceae"}
rows = collections.defaultdict(dict)
for r in cds.filter(pl.col("acc").is_in(list(phae))).iter_rows(named=True):
    p = parse(r["loc"]); g = (r["gene"] or "").lower()
    if not p or not g: continue
    fn = f"{D}/genomes/mito/{r['acc']}.fna.gz"
    if r["acc"] not in rows or "_seq" not in rows[r["acc"]]:
        try: rows[r["acc"]]["_seq"] = "".join(l.strip() for l in gzip.open(fn, "rt") if not l.startswith(">")).upper()
        except FileNotFoundError: continue
    seq = rows[r["acc"]]["_seq"]; st_, parts = p
    nt = "".join(seq[a - 1:b] for a, b in parts)
    if st_ < 0: nt = rc(nt)
    if len(nt) % 3: continue
    c = collections.Counter(nt[i:i + 3] for i in range(3, len(nt) - 3, 3))
    n = sum(c.values())
    rows[r["acc"]][g] = (n, c["ATG"], c["ATA"], c["ATT"] + c["ATC"])
CH = "NC_085326.1"
print("gene\tn_cod\tATG\tATA\tATT+C\trel_n\trelMet_scaled_med\trelIle_scaled_med\tATG+ATA\tATT+C+ATA")
T = collections.Counter()
for g, v in sorted(rows[CH].items()):
    if g == "_seq": continue
    n, atg, ata, att = v
    rel = [v[g] for a, v in rows.items() if a != CH and g in v and v[g][0] > 0]
    if len(rel) < 10: continue
    m = st.median(r[1] / r[0] for r in rel) * n; i = st.median((r[2] + r[3]) / r[0] for r in rel) * n
    print(f"{g}\t{n}\t{atg}\t{ata}\t{att}\t{len(rel)}\t{m:.1f}\t{i:.1f}\t{atg+ata}\t{att+ata}")
    T.update(dict(n=n, atg=atg, ata=ata, att=att)); T["m"] += m; T["i"] += i
print("TOTAL", dict(T))
