"""Summarize sequential length/regression runs; reject changed output hashes."""
import json
from pathlib import Path
from statistics import median
import sys

root=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).parent/'measurements/lengths'

def rows(pattern):
    return [json.loads(p.read_text()) for p in sorted(root.glob(pattern))]

print('Power-of-two regression: median of three process medians (ms)')
for exponent in (16,18,20):
    a=rows(f'reg-control-m{exponent}-r*.json')
    b=rows(f'reg-general-m{exponent}-r*.json')
    assert len(a)==len(b)==3
    assert all(x['output_hash']==y['output_hash'] for x,y in zip(a,b))
    old,new=median(x['median_ms'] for x in a),median(x['median_ms'] for x in b)
    print(f'2^{exponent}: {old:.6f} -> {new:.6f} ({100*(new/old-1):+.2f}%)')

if list(root.glob('range-range-k*.json')):
    print('\nRuntime-length direct routing versus tiled/exact-size baseline (ms)')
    ks=sorted({x['K'] for x in rows('range-range-k*.json')})
    for k in ks:
        a=rows(f'range-general-k{k}-r*.json')
        b=rows(f'range-range-k{k}-r*.json')
        assert len(a)==len(b)==3
        assert all(x['output_hash']==y['output_hash'] for x,y in zip(a,b))
        old,new=median(x['median_ms'] for x in a),median(x['median_ms'] for x in b)
        print(f'{k:8d}: {old:.6f} -> {new:.6f} ({100*(new/old-1):+.2f}%)')
