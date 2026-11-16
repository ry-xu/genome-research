import sys, csv, subprocess
meta, runids, out = sys.argv[1:4]
ids = set(l.split('\t')[0] for l in open(runids))
p = subprocess.Popen(['zstd', '-dcq', meta], stdout=subprocess.PIPE, text=True)
with open(out, 'w') as fo:
    fo.write('id\tvotu\tvgenus\tvfam\tcheckv_trimmed\tcheckv_host_markers\tlifestyle\tquality\tlength\tgc\tcrispr_host\tcrispr_conn\tcrispr_agree\n')
    for r in csv.DictReader(p.stdout, delimiter='\t'):
        if r['uhgv_genome'] not in ids: continue
        tax = dict(x.split('-', 1) for x in r['uhgv_taxonomy'].split(';') if '-' in x)
        fo.write('\t'.join([r['uhgv_genome'], r['uhgv_votu'], 'vGENUS-' + tax.get('vGENUS', 'NA'), 'vFAM-' + tax.get('vFAM', 'NA'), r['checkv_trimmed'], r['checkv_host_markers'], r['lifestyle'], r['checkv_quality'], r['genome_length'], r['gc_percent'], r['crispr_host_lineage_gtdb_r207'], r['crispr_total_connections'], r['crispr_host_agreement']]) + '\n')
