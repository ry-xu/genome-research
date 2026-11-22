"""Homolog-free check of Lyra's Choristocarpus mito ATA=Met: per Phaeophyceae (+ other stramenopile) mito genome,
codon shares over annotated CDS (first and last codon dropped). If ATA=Met, Choristocarpus (ATG+ATA) should sit
with the relatives' Met share (ATG) and its ATT+ATC with their Ile share (ATT+ATC+ATA)."""
import re, gzip, glob, sys, collections
import polars as pl
D = "lab/work/organelle-reading/data"
# taxids with Stramenopiles lineage; keep group label
grp = {}
for f in sorted(glob.glob(D + "/raw_tax/*.xml")):
    s = open(f).read()
    for blk in re.split(r"(?m)^<Taxon>", s)[1:]:
        m = re.search(r"<TaxId>(\d+)</TaxId>\s*<ScientificName>([^<]*)</ScientificName>.*?<Lineage>([^<]*)</Lineage>", blk, re.S)
        if not m: continue
        tid, name, lin = m.groups()
        if "Stramenopiles" not in lin: continue
        g = next((x for x in ["Phaeophyceae", "Xanthophyceae", "Schizocladiophyceae", "Bacillariophyta", "Eustigmatophyceae",
                               "Pelagophyceae", "Chrysophyceae", "Raphidophyceae", "Oomycota"] if x in lin), "other")
        grp[tid] = (g, name, lin)
at = {}
for line in open(D + "/acc_taxid.tsv"):
    a, t = line.split()
    if t in grp: at[a] = t
print("stramenopile accs:", len(at), file=sys.stderr)
cds = pl.read_parquet(D + "/cds_mito.parquet").filter(pl.col("acc").is_in(list(at)) & ~pl.col("pseudo") & ~pl.col("trans_splicing"))
def parse(loc):
    strand = -1 if loc.startswith("complement(") else 1
    parts = [(int(a), int(b)) for a, b in re.findall(r"<?(\d+)\.\.>?(\d+)", loc)]
    if "complement(join" in loc or loc.startswith("complement("):
        pass
    if re.search(r"join\(complement|order", loc): return None
    return strand, parts
comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
rc = lambda s: "".join(comp.get(c, "N") for c in reversed(s))
out = []
for acc, sub in cds.group_by("acc"):
    acc = acc[0]
    fn = f"{D}/genomes/mito/{acc}.fna.gz"
    try: seq = "".join(l.strip() for l in gzip.open(fn, "rt") if not l.startswith(">")).upper()
    except FileNotFoundError: continue
    cnt = collections.Counter(); ngene = 0
    for loc in sub["loc"]:
        p = parse(loc)
        if not p: continue
        st, parts = p
        nt = "".join(seq[a - 1:b] for a, b in parts)
        if st < 0: nt = rc(nt)
        if len(nt) % 3: continue
        cods = [nt[i:i + 3] for i in range(3, len(nt) - 3, 3)]
        cnt.update(cods); ngene += 1
    tot = sum(v for k, v in cnt.items() if set(k) <= set("ACGT"))
    if tot < 3000: continue
    g, name, _ = grp[at[acc]]
    out.append((g, name, acc, ngene, tot, cnt["ATG"], cnt["ATA"], cnt["ATT"] + cnt["ATC"]))
out.sort()
print("group\tname\tacc\tgenes\tcodons\tATG/1k\tATA/1k\tATT+ATC/1k\tATA_share_of_ATA+ATT+ATC")
for g, name, acc, ng, tot, atg, ata, att in out:
    print(f"{g}\t{name}\t{acc}\t{ng}\t{tot}\t{1000*atg/tot:.1f}\t{1000*ata/tot:.1f}\t{1000*att/tot:.1f}\t{ata/(ata+att):.3f}")
