"""Summarize K18/K20 routing experiments, checking comparable output hashes."""
import argparse
import json
from pathlib import Path
import statistics
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('directory',type=Path)
a=p.parse_args()
rows={path.name:json.loads(path.read_text()) for path in a.directory.glob('*.json')}
hashes={}
for name,row in rows.items():
    key=(row['m'],row['route_seed'],len(row['samples_ms']),row['configuration'].replace('imt_','').replace('_weight5_seed0_r1','').replace('_weight5',''))
    hashes.setdefault(key,set()).add(row['output_hash'])
assert all(len(values)==1 for values in hashes.values()),hashes
if any('mode' in name for name in rows):
    for mode in range(8):
        for kind in ('spin_benchmark','spin_bidirectional_benchmark'):
            group=[r for name,r in rows.items() if name.startswith(f'{kind}-mode{mode}-')]
            if group: print(kind,mode,statistics.median(r['median_ms'] for r in group),[r['median_ms'] for r in group])
    for name,row in sorted(rows.items()):
        if name.startswith('tile'): print(name,row['median_ms'])
else:
    for seed in (1,17):
        for kind in ('spin_benchmark','spin_bidirectional_benchmark'):
            for m in (16,18,20):
                pair=[]
                for case in ('control','selected'):
                    group=[r for name,r in rows.items() if name.startswith(f'{case}-{kind}-m{m}-s{seed}-')]
                    if group:
                        median=statistics.median(r['median_ms'] for r in group)
                        pair.append(median)
                        print(case,kind,m,seed,median,'setup',group[0]['retained_setup_bytes'],'workspace',group[0]['workspace_bytes'])
                if len(pair)==2: print('latency_reduction_percent',100*(1-pair[1]/pair[0]))
print('Equivalent outputs across all matched variants: PASS')
