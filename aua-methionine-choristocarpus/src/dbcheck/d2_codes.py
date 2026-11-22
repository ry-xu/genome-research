"""(a) NCBI Datasets taxonomy report: genetic codes for target + comparison taxa.
(b) NCBI gc.prt: which translation tables read AUA as Met; (c) wprintgc page: organisms listed for those tables."""
import json, re, html
D = "lab/work/organelle-reading/checks/chor_wgs/dbcheck/data/"
r = json.load(open(D + "datasets_taxreport.json"))
print("== (a) NCBI Datasets v2 taxonomy dataset_report (query 116065,2870,2833,35131,240563)")
for rep in r["reports"]:
    t = rep["taxonomy"]; gc = t.get("genetic_code", {})
    print(f"taxid {t['tax_id']} {t['current_scientific_name']['name']} ({t['rank']}): "
          f"primary={gc.get('primary')} mito={gc.get('mitochondrial')} plastid={gc.get('plastid')}")
    print("   keys in genetic_code:", sorted(gc.keys()))
print("\n== (b) gc.prt tables with AUA (ATA) = M")
txt = open(D + "gc.prt").read()
bases = {}
for m in re.finditer(r'\{\s*name "([^"]+)"\s*(?:,\s*name "([^"]+)"\s*)?,\s*id (\d+)\s*,\s*ncbieaa\s+"([^"]+)"\s*,\s*sncbieaa\s+"([^"]+)"', txt):
    name, alt, gid, aa, st = m.groups()
    # codon order TCAG, base1/base2/base3
    order = "TCAG"; codons = [a+b+c for a in order for b in order for c in order]
    t = dict(zip(codons, aa)); s = dict(zip(codons, st))
    print(f"table {gid:>2}: ATA={t['ATA']} ATG={t['ATG']} ATA-start={s['ATA']}  {name}{' / '+alt if alt else ''}"
          + ("   <== ATA=Met" if t['ATA']=='M' else ""))
print("\n== (c) wprintgc.html: text for tables with ATA=M, searching for stramenopile/brown-algal terms")
w = html.unescape(re.sub(r"<[^>]+>", " ", open(D + "wprintgc.html", encoding="utf-8", errors="replace").read()))
w = re.sub(r"[ \t]+", " ", w)
terms = ["Phaeophyceae", "brown alga", "Choristocarpus", "Discosporangi", "stramenopile", "Stramenopile", "Xanthophyceae",
         "yellow-green", "Labyrinthul", "heterokont", "Heterokont", "Ochrophyta"]
for term in terms:
    hits = [m.start() for m in re.finditer(term, w)]
    print(f"term '{term}': {len(hits)} hits")
    for h in hits[:3]:
        print("    ..." + w[max(0,h-150):h+150].replace("\n"," ") + "...")
print("\n== (d) wprintgc.html sanity: does it carry the table text?")
for term in ["Vertebrate Mitochondrial", "Thraustochytrium", "Ascidian", "AUA", "ATA", "Systematic Range", "alga"]:
    print(f"term '{term}': {len(re.findall(term, w))} hits")
i = w.find("Thraustochytrium Mitochondrial Code")
print("Thraustochytrium section excerpt:", w[i:i+600].replace("\n", " ") if i >= 0 else "not found")
hdrs = [(m.start(), m.group(1), m.group(2)) for m in re.finditer(r"(\d+)\. ([A-Z][A-Za-z ,()-]{3,90}? Code)", w)]
for sr in [m.start() for m in re.finditer("Systematic Range", w)]:
    prev = [h for h in hdrs if h[0] < sr]
    gid, nm = (prev[-1][1], prev[-1][2]) if prev else ("?", "?")
    if gid in ("2","3","5","13","21","23","4"):
        print(f"\n-- table {gid} ({nm}) systematic range:", w[sr:sr+400].replace("\n"," "))
print("\n== (e) for each ATA=Met table, every sense/stop difference from table 1 (is any table 'standard + AUA=Met' only?)")
order = "TCAG"; codons = [a+b+c for a in order for b in order for c in order]
tabs = {}
for m in re.finditer(r'\{\s*name "([^"]+)"\s*(?:,\s*name "([^"]+)"\s*)?,\s*id (\d+)\s*,\s*ncbieaa\s+"([^"]+)"', txt):
    tabs[int(m.group(3))] = dict(zip(codons, m.group(4)))
for gid, t in sorted(tabs.items()):
    if t["ATA"] != "M": continue
    diffs = [f"{c}:{tabs[1][c]}->{t[c]}" for c in codons if t[c] != tabs[1][c]]
    print(f"   table {gid}: {diffs}")
only = [gid for gid, t in tabs.items() if [c for c in codons if t[c] != tabs[1][c]] == ["ATA"]]
print("   tables differing from table 1 only at ATA:", only if only else "none")
