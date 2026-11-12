"""AnMAGC global anammox genome catalogue (Wang et al. 2024 Water Research 10.1016/j.watres.2024.122356;
preprint Research Square rs-4257542 v1 Additional file 2: Table S3 = 1,768 strain-level MAGs with GTDB r207 + r214
(GTDB-Tk v2.1.1); FASTA in figshare 10.6084/m9.figshare.25476583 Strain_level_MAGs.zip (file 45271516)).
Environment = sample with max relative abundance (Table S4, species-level MAGs) joined to Table S1.
Usage: data_dir out_tsv"""
import csv, sys, os, re
d, out = sys.argv[1], sys.argv[2]
TARGET = {"JADLFU01", "JANWXX01", "DASTPP01", "JAXLYQ01", "CAYNEK01"}
R232 = set(l.strip() for l in open("data/outside_gtdb/r232_polyangiaceae_genera.txt") if l.strip())
rd = lambda f: list(csv.reader(open(os.path.join(d, f)), delimiter="\t"))
s3 = rd("AddF2.Table_S3.tsv"); h3 = s3[0]
s1 = rd("AddF2.Table_S1.tsv"); h1 = s1[0]; samples = {r[0]: dict(zip(h1, r)) for r in s1[1:]}
s4 = rd("AddF2.Table_S4.tsv"); h4 = s4[0]; ab = {r[0]: r[1:] for r in s4[1:]}
genus = lambda t: (re.search(r";g__([^;]*)", t) or [None, ""])[1]
print("total strain-level genomes", len(s3) - 1, "of which GCA refs", sum(r[0].startswith("GC") for r in s3[1:]))
n = t = dl = 0
with open(out, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["genome_id", "taxonomy", "completeness", "contamination", "environment", "ncbi_ena_accession", "fasta_url",
                "taxonomy_r207", "genus_r214", "target_genus", "download_set", "species_rep", "n_contigs", "size", "gtdb_version"])
    for r in s3[1:]:
        rec = dict(zip(h3, r))
        t214, t207 = rec["Taxonomy_GTDB214"], rec["Taxonomy_GTDB207"]
        if "f__Polyangiaceae" not in t214 and "f__Polyangiaceae" not in t207: continue
        n += 1
        gid = rec["Strain_level (99% identity)"]; g = genus(t214)
        tg = g in TARGET or genus(t207) in TARGET; t += tg
        isgca = gid.startswith("GC")
        dset = "a" if tg else ("b" if (g == "" or g not in R232) else "")
        if dset and not isgca: dl += 1
        rep = rec["Species_level (95% identity)"]
        env = ""
        a = ab.get(rep) or ab.get(gid)
        if a:
            vals = [float(x) if x not in ("", "NA") else 0.0 for x in a]
            i = max(range(len(vals)), key=vals.__getitem__); sname = h4[1 + i]; s = samples.get(sname, {})
            env = f"anammox-system metagenome; max-abundance sample {sname} ({vals[i]:.3g}%): {s.get('Process','')}, {s.get('Bioreactor','')}, {s.get('Accession_number','')}, {s.get('Reference','')}"
        w.writerow([gid, t214, rec["completeness"], rec["contamination"], env, gid if isgca else "",
                    "" if isgca else "https://ndownloader.figshare.com/files/45271516 (zip member Strain_level_MAGs/%s.fa)" % gid,
                    t207, g, "yes" if tg else "no", dset + ("_ncbi" if dset and isgca else ""), rep, rec["number_contigs"], rec["size"],
                    "GTDB r214 (primary) and r207, GTDB-Tk v2.1.1"])
print("polyangiaceae", n, "target", t, "to download", dl)
