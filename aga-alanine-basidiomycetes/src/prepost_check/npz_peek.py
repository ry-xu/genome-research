import sys, numpy as np, zipfile
for p in sys.argv[1:]:
    z = zipfile.ZipFile(p)
    for i in z.infolist():
        print(p.split('/')[-1], i.filename, i.file_size)
