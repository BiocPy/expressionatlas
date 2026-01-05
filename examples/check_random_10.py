"""Find 10 random experiments via API and check TSV/SDRF availability."""
import urllib.request
import json
from expression_atlas.download import _try_download_sdrf

# Step 1: Get 10 experiments from BioStudies API
print("=" * 80)
print("STEP 1: Fetching 10 experiments from BioStudies API...")
print("=" * 80)

api_url = "http://www.ebi.ac.uk/biostudies/api/v1/search?query=&gxa=TRUE&pageSize=10"
with urllib.request.urlopen(api_url, timeout=30) as f:
    data = json.loads(f.read().decode('utf-8'))

accessions = [hit['accession'] for hit in data['hits']]
print(f"Found {len(accessions)} experiments: {accessions}\n")

# Step 2: Check file availability for each
print("=" * 80)
print("STEP 2: Checking file availability on FTP...")
print("=" * 80)

base = 'ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments'

def check_url(url):
    try:
        with urllib.request.urlopen(url, timeout=15) as f:
            f.read(100)
        return True
    except:
        return False

results = []
for acc in accessions:
    tsv_url = f"{base}/{acc}/{acc}-raw-counts.tsv"
    sdrf_url = f"{base}/{acc}/{acc}.condensed-sdrf.tsv"
    
    has_tsv = check_url(tsv_url)
    has_sdrf = check_url(sdrf_url)
    
    results.append({
        'accession': acc,
        'has_tsv': has_tsv,
        'has_sdrf': has_sdrf,
    })
    print(f"  {acc}: TSV={'YES' if has_tsv else 'NO':3}, SDRF={'YES' if has_sdrf else 'NO'}")

# Step 3: Try to download/parse for experiments that have files
print("\n" + "=" * 80)
print("STEP 3: Testing download/parse for experiments with files...")
print("=" * 80)

for r in results:
    acc = r['accession']
    
    if r['has_sdrf']:
        sdrf_url = f"{base}/{acc}/{acc}.condensed-sdrf.tsv"
        df = _try_download_sdrf(sdrf_url)
        r['sdrf_parsed'] = df is not None
        r['sdrf_shape'] = f"{df.shape[0]}x{df.shape[1]}" if df is not None else "N/A"
    else:
        r['sdrf_parsed'] = False
        r['sdrf_shape'] = "N/A"

# Step 4: Summary table
print("\n" + "=" * 80)
print("SUMMARY TABLE")
print("=" * 80)
print(f"{'Accession':<15} | {'TSV Data':<10} | {'SDRF Meta':<10} | {'SDRF Parsed':<12} | {'Shape':<12} | {'Status'}")
print("-" * 85)

for r in results:
    tsv = "YES" if r['has_tsv'] else "NO"
    sdrf = "YES" if r['has_sdrf'] else "NO"
    parsed = "YES" if r['sdrf_parsed'] else "NO"
    shape = r['sdrf_shape']
    
    if r['has_tsv'] and r['sdrf_parsed']:
        status = "FULL (data+meta)"
    elif r['has_tsv']:
        status = "DATA ONLY"
    elif r['sdrf_parsed']:
        status = "META ONLY (needs R)"
    else:
        status = "NONE"
    
    print(f"{r['accession']:<15} | {tsv:<10} | {sdrf:<10} | {parsed:<12} | {shape:<12} | {status}")

# Final counts
print("\n" + "=" * 80)
print("COUNTS")
print("=" * 80)
both = sum(1 for r in results if r['has_tsv'] and r['sdrf_parsed'])
tsv_only = sum(1 for r in results if r['has_tsv'] and not r['sdrf_parsed'])
sdrf_only = sum(1 for r in results if not r['has_tsv'] and r['sdrf_parsed'])
none = sum(1 for r in results if not r['has_tsv'] and not r['sdrf_parsed'])

print(f"  Both TSV + SDRF:  {both}/10  (can use TSV fallback with full metadata)")
print(f"  TSV only:         {tsv_only}/10  (can use TSV fallback, no metadata)")
print(f"  SDRF only:        {sdrf_only}/10  (needs R for data, has metadata)")
print(f"  Neither:          {none}/10  (needs R)")
