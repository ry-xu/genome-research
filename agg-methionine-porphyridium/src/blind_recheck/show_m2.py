import sys, json
d = json.load(open(sys.argv[1]))
print(d['name'], 'len', d['genome_len'], 'sites', d['n_sites'], 'moved_pass1', d['moved_pass1'])
print('genes', sorted(d['genes_hit']))
full = len(sys.argv) > 2
for c, r in d['pass2'].items():
    flag = ''
    if r['strong']: flag = 'STRONG'
    elif r['support']: flag = 'support'
    if r.get('absent_flag'): flag += ' ABSENT'
    if full or flag or r['std'] == '*' or r['best'] != r['std']:
        ls = '%.1f' % r['llr_std'] if r['llr_std'] is not None else 'n/a'
        print('%s std=%s n=%4d best=%s 2nd=%s LLRstd=%7s LLR2nd=%6.1f top=[%s] %s' % (c, r['std'], r['n'], r['best'], r['second'], ls, r['llr_2nd'], r['top_ref'], flag))
