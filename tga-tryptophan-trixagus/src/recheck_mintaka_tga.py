import glob
names={'GCA_040210545.2':'Tser','GCA_040211365.2':'Tder','GCA_040409785.2':'Ttur','GCA_040211755.2':'Tcar','GCA_040210955.2':'Tduv','GCA_037094695.1':'Tmey','GCA_040212085.2':'AUS03'}
def fa(f):
    s={};k=None
    for l in open(f):
        if l[0]=='>': k=l[1:].split()[0]; s[k]=[]
        else: s[k].append(l.strip())
    return {k:''.join(v).upper() for k,v in s.items()}
for a,n in names.items():
    G=fa(a+'.g4.ffn'); cod=tga=tgg=0
    for s in G.values():
        c=[s[i:i+3] for i in range(0,len(s)-2,3)][:-1]  # internal codons incl. start, excl. terminal stop
        cod+=len(c); tga+=c.count('TGA'); tgg+=c.count('TGG')
    print(f'{n} genes {len(G)} codons {cod} TGA/1000 {1000*tga/cod:.2f} TGA {tga} TGG {tgg} TGA share of Trp {100*tga/(tga+tgg):.1f}%')
