"""Test if pyreadr and rdata packages can read Expression Atlas .Rdata files."""
import urllib.request
import tempfile
import os

# Download a sample .Rdata file
acc = 'E-MTAB-7841'
url = f'ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments/{acc}/{acc}-atlasExperimentSummary.Rdata'

print('Downloading .Rdata file...')
with urllib.request.urlopen(url, timeout=60) as f:
    data = f.read()

tmp_file = tempfile.NamedTemporaryFile(suffix='.Rdata', delete=False)
tmp_file.write(data)
tmp_file.close()
print(f'Saved to: {tmp_file.name}')
print(f'Size: {len(data)} bytes')

# Try pyreadr
print('\n--- Testing pyreadr ---')
try:
    import pyreadr
    result = pyreadr.read_r(tmp_file.name)
    print(f'SUCCESS! Keys: {list(result.keys())}')
    for k, v in result.items():
        print(f'  {k}: {type(v)}')
except Exception as e:
    print(f'FAILED: {type(e).__name__}: {e}')

# Try rdata
print('\n--- Testing rdata ---')
try:
    import rdata
    parsed = rdata.parser.parse_file(tmp_file.name)
    print(f'Parsed OK. Type: {type(parsed)}')
    converted = rdata.conversion.convert(parsed)
    print(f'Converted OK. Type: {type(converted)}')
    if isinstance(converted, dict):
        print(f'Keys: {list(converted.keys())}')
        for k, v in converted.items():
            print(f'  {k}: {type(v)}')
except Exception as e:
    print(f'FAILED: {type(e).__name__}: {e}')

os.unlink(tmp_file.name)
