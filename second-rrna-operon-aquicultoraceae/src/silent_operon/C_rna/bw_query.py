"""POST one sourmash signature (k=21, scaled=1000; the first sig of a .sig file) to branchwater and save the JSON result.

usage: python3 -I src/bw_query.py SIG OUT.json [THRESHOLD]
Branchwater (branchwater.sourmash.bio, formerly branchwater.jgi.doe.gov) searches ~1.16 M SRA metagenome runs by
containment of the query's k-mers. The payload is {"signatures": <one signature as a JSON string>, "threshold": t}.
"""
import json, sys, urllib.request

sig = json.dumps(json.load(open(sys.argv[1]))[0])
thr = float(sys.argv[3]) if len(sys.argv) > 3 else 0.1
req = urllib.request.Request('https://branchwater.sourmash.bio/', data=json.dumps({'signatures': sig, 'threshold': thr}).encode(),
                             headers={'Content-Type': 'application/json'}, method='POST')
txt = urllib.request.urlopen(req, timeout=900).read().decode()
open(sys.argv[2], 'w').write(txt)
print(sys.argv[1], len(json.loads(txt)), 'hits')
