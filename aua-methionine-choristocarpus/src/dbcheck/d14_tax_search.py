"""Counts/IDs of NCBI taxonomy esearch results saved in data/tax_search_*.xml."""
import glob, re, os
for f in sorted(glob.glob("lab/work/organelle-reading/checks/chor_wgs/dbcheck/data/tax_search_*.xml")):
    t = open(f).read()
    c = re.search(r"<Count>(\d+)</Count>", t).group(1); ids = re.findall(r"<Id>(\d+)</Id>", t)
    q = re.search(r"<QueryTranslation>(.*?)</QueryTranslation>", t)
    print(f"{os.path.basename(f):45s} count={c} ids={ids} query={q.group(1) if q else ''}")
