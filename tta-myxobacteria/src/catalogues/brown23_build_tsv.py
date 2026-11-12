"""Brown et al. 2023 activated-sludge bioreactor MAGs (Zenodo 10028566, BROWN23-MAGs.tar.gz, 2023-10-20).
Reads GTDB-Tk summary + checkm lineage output from the extracted dir; writes Polyangiaceae TSV and
gzips target/(b) FASTAs into outdir. Usage: extracted_dir out_tsv fasta_outdir"""
import csv, sys, os, re, glob, gzip, shutil
src, out, fdir = sys.argv[1], sys.argv[2], sys.argv[3]
TARGET = {"JADLFU01", "JANWXX01", "DASTPP01", "JAXLYQ01", "CAYNEK01"}
R232 = set(l.strip() for l in open("data/outside_gtdb/r232_polyangiaceae_genera.txt") if l.strip())
tax = {}
sums = [p for p in glob.glob(os.path.join(src, "**", "*"), recursive=True) if os.path.isfile(p) and re.search(r"(gtdb|summary|classif)", os.path.basename(p), re.I)]
print("taxonomy files:", sums)
for p in sums:
    with open(p, errors="replace") as fh:
        rd = csv.reader(fh, delimiter="\t")
        hdr = next(rd, [])
        if "user_genome" in hdr and "classification" in hdr:
            iu, ic = hdr.index("user_genome"), hdr.index("classification")
            for r in rd:
                if len(r) > ic: tax[r[iu]] = r[ic]
cm = {}
for p in glob.glob(os.path.join(src, "**", "checkm*"), recursive=True):
    for line in open(p, errors="replace"):
        f = line.split()
        if len(f) >= 14 and re.match(r"^\d", f[-1]) and not line.startswith("-") and f[0] != "Bin":
            try: cm[f[0]] = (f[-3], f[-2])
            except Exception: pass
fastas = {re.sub(r"\.(fna|fa|fasta)$", "", os.path.basename(p)): p for p in glob.glob(os.path.join(src, "**", "*.f*a"), recursive=True)}
print("genomes with taxonomy", len(tax), "fasta files", len(fastas), "checkm rows", len(cm))
n = t = dl = 0
with open(out, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["genome_id", "taxonomy", "completeness", "contamination", "environment", "ncbi_ena_accession", "fasta_url", "genus", "target_genus", "download_set", "gtdb_version"])
    for gid, tx in sorted(tax.items()):
        if "f__Polyangiaceae" not in tx: continue
        n += 1
        g = (re.search(r";g__([^;]*)", tx) or [None, ""])[1]
        tg = g in TARGET; t += tg
        dset = "a" if tg else ("b" if (g == "" or g not in R232) else "")
        if dset and gid in fastas:
            shutil.copyfileobj(open(fastas[gid], "rb"), gzip.open(os.path.join(fdir, gid + ".fna.gz"), "wb")); dl += 1
        c = cm.get(gid, ("", ""))
        w.writerow([gid, tx, c[0], c[1], "activated sludge bioreactors (Brown et al. 2023)", "",
                    "https://zenodo.org/records/10028566/files/BROWN23-MAGs.tar.gz (member BROWN23-MAGs/%s.fna)" % gid,
                    g, "yes" if tg else "no", dset, "GTDB-Tk, release not stated (230105_classifications.tsv, dated 2023-01-05, so r207 or earlier)"])
print("polyangiaceae", n, "target", t, "downloaded", dl)
