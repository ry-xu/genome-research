import sys, time, urllib.request, urllib.parse, json
# usage: esearch.py "<term1>" "<term2>" ...  -> prints accession, length, title of refseq/genbank mito hits
base="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
for term in sys.argv[1:]:
    q=urllib.parse.urlencode({"db":"nuccore","term":term,"retmax":"8","retmode":"json"})
    r=json.load(urllib.request.urlopen(base+"esearch.fcgi?"+q)); time.sleep(0.4)
    ids=r["esearchresult"]["idlist"]
    print("##",term,r["esearchresult"]["count"])
    if not ids: continue
    q=urllib.parse.urlencode({"db":"nuccore","id":",".join(ids),"retmode":"json"})
    s=json.load(urllib.request.urlopen(base+"esummary.fcgi?"+q)); time.sleep(0.4)
    for i in ids:
        d=s["result"][i]; print(d["accessionversion"],d["slen"],d["title"])
