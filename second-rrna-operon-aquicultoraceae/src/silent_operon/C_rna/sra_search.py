"""Query NCBI SRA via E-utilities and print runinfo rows (Run, BioProject, LibrarySource, LibraryStrategy, LibrarySelection, spots, bases, ScientificName, SampleName).

usage: python3 -I src/sra_search.py 'QUERY' [OUT.csv]
"""
import csv, io, sys, time, urllib.parse, urllib.request

E = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'


def get(url):
    for i in range(4):
        try:
            return urllib.request.urlopen(url, timeout=60).read().decode()
        except Exception:
            time.sleep(2 + 3 * i)
    raise SystemExit('fetch failed: ' + url)


def main(q, out=None):
    r = get(E + 'esearch.fcgi?db=sra&usehistory=y&retmax=0&term=' + urllib.parse.quote(q))
    import re
    cnt = int(re.search(r'<Count>(\d+)</Count>', r).group(1))
    if cnt == 0:
        print(f'# {q!r}: 0 records'); return
    we = re.search(r'<WebEnv>(\S+?)</WebEnv>', r).group(1); qk = re.search(r'<QueryKey>(\d+)</QueryKey>', r).group(1)
    txt = get(E + f'efetch.fcgi?db=sra&rettype=runinfo&retmode=text&query_key={qk}&WebEnv={we}&retmax=5000')
    rows = [x for x in csv.DictReader(io.StringIO(txt)) if x.get('Run') and x['Run'] != 'Run']
    print(f'# {q!r}: {cnt} records, {len(rows)} runs')
    for x in rows:
        print('\t'.join(x.get(k, '') for k in ('Run', 'BioProject', 'BioSample', 'LibrarySource', 'LibraryStrategy', 'LibrarySelection', 'spots', 'bases', 'ScientificName', 'SampleName', 'ReleaseDate')))
    if out:
        with open(out, 'w') as f: f.write(txt)


if __name__ == '__main__':
    main(*sys.argv[1:])
