"""tRNAdb 2026 (tdb.bioinf.uni-leipzig.de, successor of mitotRNAdb) query results: Choristocarpus tRNAs; Phaeophyceae CAT tRNAs."""
import json, collections
D = "lab/work/organelle-reading/checks/chor_wgs/dbcheck/data/"
r = json.load(open(D + "tdb_r1.json"))
print(f"== q1 organism='Choristocarpus tenellus', databases mttrna,mttrna_rna,trna,trna_rna: total_results={r['total_results']} returned={len(r['results'])}")
for x in sorted(r["results"], key=lambda x: x["ams"]):
    print(f"   {x['ams']:5s} {x['ac']}  db={x['database']}  ref={x['ref']}  version={x['tRNAdb_version']}")
aa = collections.Counter(x["ams"] for x in r["results"])
print("   amino acids:", dict(aa))
print("   CAT entries:", [(x['ams'], x['ac']) for x in r['results'] if x['ac'] == 'CAT'])
print("   Ile entries:", [(x['ams'], x['ac']) for x in r['results'] if x['ams'] == 'Ile'])
r3 = json.load(open(D + "tdb_r3.json"))
print(f"\n== q3 organism='Discosporangium mesarthrocarpum': total_results={r3['total_results']}")
r2 = json.load(open(D + "tdb_r2.json"))
print(f"\n== q2 clade='Phaeophyceae', anticodon=CAT: total_results={r2['total_results']} returned={len(r2['results'])}")
c = collections.Counter((x["ams"], x["database"]) for x in r2["results"])
print("   (amino acid, db) counts:", dict(c))
per = collections.defaultdict(list)
for x in r2["results"]:
    per[x["tax_name"]].append(x["ams"])
print(f"   organisms: {len(per)}")
for k in sorted(per):
    print(f"   {k:40s} {sorted(per[k])}")
