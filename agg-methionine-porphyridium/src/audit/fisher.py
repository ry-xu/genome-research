"""Which 2x2 table gives the post's 'one-sided Fisher p 7e-14'? (audit)"""
from scipy.stats import fisher_exact
tabs = {
 'ATG vs not-ATG among changed (AGG 17/21, Arg 0/37)': [[17, 4], [0, 37]],
 'ATG vs Arg among changed (AGG 17 ATG/0 Arg, Arg 0 ATG/34 Arg)': [[17, 0], [0, 34]],
 'Arg vs not among changed (AGG 0/21, Arg 34/37)': [[0, 21], [34, 3]],
 'ATG-or-not, Arg-or-not as 17/21 vs 3/37 (post wording: 17 of 21 vs 34 of 37)': [[17, 4], [3, 34]],
 'audit nt-alignment: ATG vs not (16/20, 0/37)': [[16, 4], [0, 37]],
 'audit nt-alignment: 16/20 vs 3/37': [[16, 4], [3, 34]],
 'audit masked-prot: 16/21 vs 3/37': [[16, 5], [3, 34]],
}
for k, t in tabs.items():
    print(f'{k}: {t}  one-sided p = {fisher_exact(t, alternative="greater")[1]:.2g}, two-sided p = {fisher_exact(t)[1]:.2g}')
