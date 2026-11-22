"""Print TMT reporter-ion intensities (126-131 region, 10 ppm) for given scans. usage: reporters.py file.mzML scan [scan...]"""
import sys
from pyteomics import mzml
import numpy as np
REP = {'126': 126.127726, '127N': 127.124761, '127C': 127.131081, '128N': 128.128116, '128C': 128.134436,
       '129N': 129.131471, '129C': 129.137790, '130N': 130.134825, '130C': 130.141145, '131': 131.138180}
want = set(sys.argv[2:])
with mzml.MzML(sys.argv[1]) as r:
    for s in r:
        sc = s['id'].split('scan=')[-1]
        if sc not in want:
            continue
        mz, it = s['m/z array'], s['intensity array']
        out = []
        for k, m in REP.items():
            w = np.abs(mz - m) <= m * 10e-6
            out.append(f'{k}:{it[w].max() if w.any() else 0:.0f}')
        print(sc, ' '.join(out))
