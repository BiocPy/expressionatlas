"""Check if experiments with TSV data also have condensed-sdrf.tsv metadata."""
import urllib.request

tests = ['E-MTAB-7841', 'E-MTAB-5214', 'E-GEOD-30352', 'E-MTAB-4344', 'E-MTAB-1733']
base = 'ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments'

def check_exists(url):
    try:
        with urllib.request.urlopen(url, timeout=15) as f:
            f.read(100)
        return True
    except:
        return False

print('Experiment        | raw-counts.tsv | condensed-sdrf.tsv')
print('-' * 58)

for acc in tests:
    tsv_url = f"{base}/{acc}/{acc}-raw-counts.tsv"
    sdrf_url = f"{base}/{acc}/{acc}.condensed-sdrf.tsv"
    
    has_tsv = "YES" if check_exists(tsv_url) else "NO"
    has_sdrf = "YES" if check_exists(sdrf_url) else "NO"
    
    print(f'{acc:17} | {has_tsv:14} | {has_sdrf}')
