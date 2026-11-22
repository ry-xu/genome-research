import numpy as np
S = 'scratch/'
C = np.load(S + 'mycache.npz')
acc = C['acc'].tolist(); tax = C['tax'].tolist()
for g in np.load(S + 'ala32.npy').tolist():
    t = tax[g].split(';')
    print(acc[g], ' | ', '; '.join(t[-4:]))
