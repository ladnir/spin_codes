"""Check paired output hashes and summarize the unchanged 128-bit path."""
import json
from pathlib import Path
from statistics import median
import sys

root=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).parent/'measurements/generic'
for exponent in (16,18,20):
    pairs=[]
    for repeat in (1,2,3):
        old=json.loads((root/f'build-length-general-m{exponent}-r{repeat}.json').read_text())
        new=json.loads((root/f'build-generic-release-m{exponent}-r{repeat}.json').read_text())
        assert old['output_hash']==new['output_hash']
        pairs.append((old['median_ms'],new['median_ms']))
    a,b=(median(x[i] for x in pairs) for i in (0,1))
    print(f'K=2^{exponent}: {a:.6f} -> {b:.6f} ms ({100*(b/a-1):+.2f}%); hashes match')
