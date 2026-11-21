from scipy.stats import fisher_exact
for lbl, t in [('538 CTA', [[70,0],[13,149]]), ('262 CTA', [[33,0],[13,166]]), ('538 CTT', [[12,0],[71,149]]), ('262 CTT',[[11,0],[35,166]])]:
    print(lbl, fisher_exact(t, alternative='greater')[1])
