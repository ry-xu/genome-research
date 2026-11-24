"""Words, gate and verdict (RULES.md §3, §4, §6 with amendment A2). Pools counts over the .done runs of each study.
usage: words.py OUT_PREFIX [--dna DNACHECK.tsv] [--reads READCHECK.tsv] [--metag METAG_WORDS.tsv] STUDY=DIR:RUN,RUN,... [STUDY=...]
  -> OUT_PREFIX.events.tsv (one row per event group x study), prints diagnostics, gate, decoy bar, verdict and leads.
Events whose +-100 nt context is identical (same locus in several genomes) are pooled into one group, named by the lowest id."""
import sys, csv, os
from collections import defaultdict
import polars as pl

args = sys.argv[1:]
out = args.pop(0)
dna = metag = rck = None
while args and args[0].startswith('--'):
    o = args.pop(0)
    if o == '--dna':
        dna = args.pop(0)
    elif o == '--metag':
        metag = args.pop(0)
    elif o == '--reads':
        rck = args.pop(0)
TACK = ['Thermoproteota', 'Korarchaeota']
CREN = ['Sulfolobales', 'Desulfurococcales', 'Thermoproteales']

# context groups
seqs, name = {}, None
for l in open('/Users/ryan/lab/work/tack-rna/out/regions.fa'):
    if l[0] == '>':
        name = l[1:].strip()
    else:
        seqs[name] = l.strip().upper()
comp = str.maketrans('ACGTN', 'TGCAN')
E = pl.read_csv('/Users/ryan/lab/work/tack-rna/out/events.tsv', separator='\t', infer_schema_length=0)
ctx = {}
for r in E.iter_rows(named=True):
    s, x, y = seqs[r['region']], int(r['rx']), int(r['ry'])
    c = s[max(0, x - 100):y + 100]
    ctx[r['id']] = (r['kind'], min(c, c.translate(comp)[::-1]))
grp_of, members = {}, defaultdict(list)
for eid in sorted(ctx):
    members[ctx[eid]].append(eid)
for k, ids in members.items():
    for e in ids:
        grp_of[e] = ids[0]
E = E.with_columns(group=pl.col('id').replace_strict(grp_of), n_members=pl.col('id').replace_strict({e: len(members[ctx[e]]) for e in ctx}))

rows = []
for arg in args:
    st, rest = arg.split('=', 1)
    d, runs = rest.split(':', 1)
    tot = defaultdict(lambda: [0, 0, 0])
    used = [r for r in runs.split(',') if os.path.exists(f'{d}/{r}.done')]
    print(f'{st}: {len(used)} of {len(runs.split(","))} runs done', file=sys.stderr)
    for r in used:
        for x in csv.DictReader(open(f'{d}/{r}.counts.tsv'), delimiter='\t'):
            t = tot[grp_of[x['event']]]
            t[0] += int(x['joined']); t[1] += int(x['through']); t[2] += int(x['other'])
    for e, (j, th, o) in tot.items():
        rows.append(dict(study=st, id=e, joined=j, through=th, other=o))
C = pl.DataFrame(rows, schema=dict(study=pl.Utf8, id=pl.Utf8, joined=pl.Int64, through=pl.Int64, other=pl.Int64))
D = C.join(E, on='id', how='left')
D = D.with_columns(share=pl.col('joined') / (pl.col('joined') + pl.col('through')).clip(1))
D = D.with_columns(word=pl.when((pl.col('joined') >= 2) & (pl.col('share') >= 0.5)).then(pl.lit('removed'))
                   .when(pl.col('joined') >= 2).then(pl.lit('partly removed'))
                   .when((pl.col('through') >= 5) & (pl.col('joined') <= 1)).then(pl.lit('read through'))
                   .otherwise(pl.lit('not covered')))
D = D.with_columns(spliced=pl.col('joined') >= 2, covered=pl.col('word') != 'not covered',
                   tack=pl.col('phylum').is_in(TACK),
                   cls2=pl.when(pl.col('kind') == 'trna').then(pl.lit('tRNA intron'))
                   .when(pl.col('kind') == 'decoy').then(pl.lit('decoy'))
                   .when(pl.col('bhb') == 'true').then(pl.when(pl.col('fam') == 'TRUB').then(pl.lit('cut BHB TRUB')).otherwise(pl.lit('cut BHB non-TRUB')))
                   .otherwise(pl.lit('cut non-BHB')))
pc = E.filter(pl.col('kind') == 'cut').select(pl.col('id').alias('parent'), pl.col('bhb').alias('p_bhb'), pl.col('fam').alias('p_fam'))
D = D.join(pc, on='parent', how='left').with_columns(
    pcls=pl.when(pl.col('p_bhb') == 'true').then(pl.when(pl.col('p_fam') == 'TRUB').then(pl.lit('cut BHB TRUB')).otherwise(pl.lit('cut BHB non-TRUB'))).otherwise(pl.lit('cut non-BHB')))
# DNA check: census 40-mer and metaG
dv = set()
if dna:
    for x in csv.DictReader(open(dna), delimiter='\t'):
        if int(x['n_genomes']) > 0:
            dv.add(grp_of[x['event']])
rc = defaultdict(lambda: [0, 0])               # group -> [joined templates, explained templates] (A3)
if rck:
    for x in csv.DictReader(open(rck), delimiter='\t'):
        t = rc[grp_of[x['event']]]
        t[0] += int(x['joined_templates']); t[1] += int(x['explained_templates'])
mg = {}
if metag:
    for x in csv.DictReader(open(metag), delimiter='\t'):
        mg[grp_of[x['event']]] = (int(x['joined']), int(x['through']))
def dna_word(eid, spliced):
    if not spliced:
        return ''
    if eid in dv:
        return 'DNA variant (census genome)'
    if eid not in rc or not rc[eid][0]:
        return 'DNA unchecked'
    j, x = rc[eid]
    if x == j:
        return 'DNA variant (read)'
    if x:
        return 'DNA mixed'
    if eid in mg:
        mj, mt = mg[eid]
        if mj >= 2 and mj / max(1, mj + mt) >= 0.1:
            return 'DNA variant (metaG)'
        if mj + mt >= 5:
            return 'DNA clean (census reads, metaG)'
    return 'DNA clean (census reads)'
D = D.with_columns(dna=pl.struct('id', 'spliced').map_elements(lambda s: dna_word(s['id'], s['spliced']), return_dtype=pl.Utf8))
D.sort('study', 'kind', 'id').write_csv(out + '.events.tsv', separator='\t')

pl.Config.set_tbl_rows(100); pl.Config.set_tbl_width_chars(220); pl.Config.set_fmt_str_lengths(60)
print('== words by study x class (all phyla, event groups)')
print(D.group_by('study', 'cls2', 'word').agg(n=pl.len()).pivot(on='word', index=['study', 'cls2'], values='n').fill_null(0).sort('study', 'cls2'))
print('== G1 diagnostic (A2): tRNA introns with joined+through >= 5 (spliced_clean: read check done, none explained)')
T = D.filter((pl.col('kind') == 'trna') & ((pl.col('joined') + pl.col('through')) >= 5))
print(T.group_by('study').agg(covered=pl.len(), spliced=pl.col('spliced').sum(), removed=(pl.col('word') == 'removed').sum(),
                              spliced_clean=(pl.col('spliced') & pl.col('dna').str.starts_with('DNA clean')).sum()).sort('study'))
print('== Gate (A2 = G2): crenarchaeal cbf5, covered ones')
G2 = D.filter((pl.col('cls2') == 'cut BHB TRUB') & pl.col('order').is_in(CREN) & pl.col('covered'))
print(G2.select('study', 'id', 'n_members', 'species', 'L', 'joined', 'through', 'word'))
gate = G2.group_by('study').agg(n=pl.len(), sp=pl.col('spliced').sum())
gated = set(gate.filter(pl.col('sp') >= 0.8 * pl.col('n'))['study'].to_list())
for st in sorted(D['study'].unique().to_list()):
    g = gate.filter(pl.col('study') == st)
    print(f"  {st}: cbf5 covered {g['n'][0] if len(g) else 0}, spliced {g['sp'][0] if len(g) else 0} -> {'GATED' if st in gated else 'ungated'}")
print('== TACK cuts and their decoys by parent class (all studies)')
TK = D.filter(pl.col('tack') & (pl.col('kind') != 'trna'))
print(TK.group_by('study', 'pcls', 'kind', 'word').agg(n=pl.len()).pivot(on='word', index=['study', 'pcls', 'kind'], values='n').fill_null(0).sort('study', 'pcls', 'kind'))
tgt = D.filter(pl.col('tack') & (pl.col('cls2') == 'cut BHB non-TRUB'))
dec = D.filter(pl.col('tack') & (pl.col('kind') == 'decoy') & (pl.col('pcls') == 'cut BHB non-TRUB'))
ns_t, ns_d = tgt['spliced'].sum(), dec['spliced'].sum()
dec_ok = ns_d <= 0.1 * ns_t if ns_t else ns_d == 0
print(f'== Decoy bar (TACK non-TRUB BHB, all studies): spliced targets {ns_t}, spliced decoys {ns_d} (decoys covered {dec["covered"].sum()}) -> {"PASS" if dec_ok else "FAIL"}')
g = tgt.filter(pl.col('study').is_in(list(gated)))
gs = g.filter(pl.col('spliced') & pl.col('dna').str.starts_with('DNA clean'))
# an event group counts once even if seen in two gated studies
sp_ids = gs['id'].unique().to_list()
cov_ids = g.filter(pl.col('covered'))['id'].unique().to_list()
ncls = gs['class'].n_unique() if len(gs) else 0
print(f'== Verdict (gated studies {sorted(gated)}): TACK non-TRUB BHB groups covered {len(cov_ids)}, spliced and DNA-clean {len(sp_ids)} in {ncls} classes; not counted: {dict(g.filter(pl.col("spliced") & ~pl.col("dna").str.starts_with("DNA clean")).group_by("dna").agg(pl.col("id").n_unique()).iter_rows())}')
if not dec_ok:
    v = 'decoys fail: no verdict'
elif len(sp_ids) >= 10 and ncls >= 2:
    v = 'HEADLINE: archaeal protein-gene splicing is common'
elif 2 <= len(sp_ids) <= 9 or (len(sp_ids) >= 10 and ncls < 2):
    v = 'MEDIUM'
elif len(cov_ids) >= 10 and len(sp_ids) <= 1:
    v = 'NOT SPLICED (kill)'
else:
    v = 'NOT TESTABLE (coverage)'
print('   verdict:', v)
print('== Spliced TACK non-TRUB BHB events, all studies (ungated ones are leads only)')
L = tgt.filter(pl.col('spliced')).with_columns(gated=pl.col('study').is_in(list(gated)))
print(L.select('study', 'gated', 'id', 'n_members', 'class', 'genus', 'fam', 'L', 'k', 'strict', 'joined', 'through', 'word', 'dna').sort('study', 'id'))
print('== Other spliced cuts (everything except TACK non-TRUB BHB), reported only')
O = D.filter((pl.col('kind') == 'cut') & pl.col('spliced') & ~(pl.col('tack') & (pl.col('cls2') == 'cut BHB non-TRUB')))
print(O.select('study', 'id', 'phylum', 'class', 'order', 'fam', 'cls2', 'L', 'k', 'joined', 'through', 'dna').sort('study', 'id'))
