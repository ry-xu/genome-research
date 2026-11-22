import sys, time, urllib.request, urllib.parse, os
# usage: efetch.py outdir rettype acc1 acc2 ...
base="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
out, rt = sys.argv[1], sys.argv[2]
for acc in sys.argv[3:]:
    q=urllib.parse.urlencode({"db":"nuccore","id":acc,"rettype":rt,"retmode":"text"})
    data=urllib.request.urlopen(base+q).read()
    open(os.path.join(out,f"{acc}.{rt}"),"wb").write(data)
    print(acc, rt, len(data)); time.sleep(0.4)
