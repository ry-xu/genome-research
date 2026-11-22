"""NCBI Taxonomy: genetic codes (nuclear GC, mito MGC, plastid PGC) and lineage for target taxa."""
import sys, xml.etree.ElementTree as ET
D = "lab/work/organelle-reading/checks/chor_wgs/dbcheck/data/"
root = ET.parse(D + "taxonomy.xml").getroot()
for t in root.findall("Taxon"):
    g = lambda p: (t.findtext(p) or "").strip()
    print(f"taxid {g('TaxId')}  {g('ScientificName')}  rank={g('Rank')}")
    print(f"  GeneticCode   id={g('GeneticCode/GCId')} name={g('GeneticCode/GCName')}")
    print(f"  MitoGenCode   id={g('MitoGeneticCode/MGCId')} name={g('MitoGeneticCode/MGCName')}")
    pg = t.find('PlastidGeneticCode')
    if pg is not None:
        print(f"  PlastidGenCode id={g('PlastidGeneticCode/PGCId')} name={g('PlastidGeneticCode/PGCName')}")
    # inherited flags are in the full taxonomy dump, not in efetch XML; print any child tags mentioning Inherit
    other = [c.tag for c in t if 'nherit' in c.tag or 'Flag' in c.tag]
    print(f"  inherit/flag tags in XML: {other if other else 'none present'}")
    print(f"  Lineage: {g('Lineage')}")
    for k in ('CreateDate','UpdateDate','PubDate'):
        print(f"  {k}: {g(k)}")
    # lineage codes are not given per ancestor; list ancestor names with rank
    lx = t.find('LineageEx')
    if lx is not None:
        print("  LineageEx ranks: " + "; ".join(f"{a.findtext('ScientificName')}({a.findtext('Rank')},{a.findtext('TaxId')})" for a in lx.findall('Taxon')[-6:]))
