import sys, numpy as np
for p in sys.argv[1:]:
    d = np.load(p)
    print(p, {k: d[k].shape for k in d.files})
    print(d['models'].tolist()[:200])
