"""tRNAdb / mitotRNAdb (tdb.bioinf.uni-leipzig.de): search the organism and clade autocomplete lists for target names."""
import re
D = "lab/work/organelle-reading/checks/chor_wgs/dbcheck/data/"
for fn in ["tdb_autocomplete_org_data.js", "tdb_autocomplete_clade_data.js"]:
    t = open(D + fn, encoding="utf-8", errors="replace").read()
    items = re.findall(r'"([^"]+)"', t)
    print(f"== {fn}: {len(items)} quoted items")
    for term in ["Choristocarpus", "Discosporangium", "Phaeophyceae", "Phaeophyta", "Sphacelari", "Saccharina", "Laminaria",
                 "Desmarestia", "Fucus", "Ectocarpus", "Pylaiella", "Dictyota", "Xanthophyceae", "Vaucheria", "Sphaerosorus",
                 "Thraustochytri", "Labyrinthul", "Stramenopiles", "Ochrophyta"]:
        hits = [x for x in items if term.lower() in x.lower()]
        print(f"   {term:16s} {len(hits):4d}  {hits[:6]}")
