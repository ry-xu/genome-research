"""GtRNAdb browse.html (all genomes): look for brown algae / stramenopiles."""
import re, html
D = "lab/work/organelle-reading/checks/chor_wgs/dbcheck/data/"
t = html.unescape(re.sub(r"<[^>]+>", " ", open(D + "gtrnadb_browse.html", encoding="utf-8", errors="replace").read()))
links = re.findall(r'href="([^"]+)"', open(D + "gtrnadb_browse.html", encoding="utf-8", errors="replace").read())
print("index text length:", len(t), "links:", len(links))
for term in ["Choristocarpus", "Discosporangium", "Ectocarpus", "Saccharina", "Phaeophyceae", "brown", "Nannochloropsis",
             "Phaeodactylum", "Thalassiosira", "Aureococcus", "Phytophthora", "Blastocystis", "Stramenopile", "mitochondri"]:
    hits = sorted(set(m.group(0) for m in re.finditer(term + r"[\w ._-]{0,30}", t, flags=re.I)))
    print(f"   {term:16s} {len(hits)} {hits[:4]}")
euk = [l for l in links if "genomes/eukaryota/" in l]
print("eukaryota genome links:", len(set(euk)))
for term in ["Ectocarpus", "Esili", "Saccharina", "Sjapo", "Choristocarpus", "Ctene", "Discosporangium", "Undaria"]:
    print(f"   link containing '{term}':", sorted(set(l for l in links if term.lower() in l.lower()))[:5])
print("sanity terms:", {k: len(re.findall(k, t)) for k in ["Plasmodium", "Saccharomyces", "Homo sapiens", "Chlamydomonas", "Phaeo", "Ochrophyta", "Oomycet"]})
print("'Phaeo' contexts:", [t[m.start()-40:m.start()+40].replace("\n"," ") for m in re.finditer("Phaeo", t)])
