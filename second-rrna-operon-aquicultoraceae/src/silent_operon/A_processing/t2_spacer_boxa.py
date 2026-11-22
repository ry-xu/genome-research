"""T2: boxA in the 16S-23S spacer with a 200-shuffle null (dinucleotide), and boxA core TCTTTAA positions."""
import random, sys
sys.path.insert(0, 'src')
from common import operons, region, dishuffle
from t2_upstream_lib import best_boxa
rng = random.Random(4); seen = set()
print('genome\tcontig\tcls\tspacer_len\tboxA_mm\tboxA\tpos_after_16S\tp_null')
for r in operons():
    if not (r['s16b'] and r['l23a']): continue
    sp = region(r, r['s16b'] + 1, r['l23a'] - 1)
    if sp in seen and r['genome'] not in ('ECOLI',): continue
    seen.add(sp)
    d, i, m = best_boxa(sp)
    nb = [best_boxa(dishuffle(sp, rng))[0] for _ in range(200)]
    p = (sum(x <= d for x in nb) + 1) / 201
    print(f"{r['genome']}\t{r['contig']}\t{r['cls']}\t{len(sp)}\t{d}\t{m}\t{i}\t{p:.3f}")
