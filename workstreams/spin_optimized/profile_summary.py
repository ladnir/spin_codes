"""Summarize diagnostic stages without substituting them for production timing."""
from pathlib import Path
import json
import statistics
import sys
p=Path(sys.argv[1]);groups={};hashes={}
for file in sorted(p.glob('*.json')):
    if file.name.endswith('-stages.json'): continue
    r=json.loads(file.read_text())
    key=(r['route_seed'],len(r['samples_ms']))
    hashes.setdefault(key,set()).add(r['output_hash'])
    if file.name=='sampled.json':continue
    kind='bidirectional' if 'bidirectional' in file.name else 'transpose-only'
    if file.name.startswith('control-'):
        groups.setdefault((kind,'control'),[]).append(r['median_ms']);continue
    stage=json.loads(file.with_name(file.stem+'-stages.json').read_text())
    assert stage['profile_calls']==len(r['samples_ms'])
    for item in stage['stages']:
        assert item['counter_running_fraction']>.999,('multiplexed',file,item)
        groups.setdefault((kind,r['tile_rows'],item['name']),[]).append(item)
    groups.setdefault((kind,r['tile_rows'],'total'),[]).append(r['median_ms'])
assert all(len(v)==1 for v in hashes.values()),hashes
for key,rows in groups.items():
    if isinstance(rows[0],dict):
        ms=statistics.mean(r['mean_ms'] for r in rows)
        cycles=statistics.mean(r['cycles'] for r in rows)
        ins=statistics.mean(r['instructions'] for r in rows)
        misses=statistics.mean(r['cache_misses'] for r in rows)
        print(key,'mean_ms',round(ms,5),'cycles_M',round(cycles/1e6,3),
              'instructions_M',round(ins/1e6,3),'IPC',round(ins/cycles,3),
              'generic_cache_misses',round(misses))
    else: print(key,'median_ms',round(statistics.median(rows),6),'samples',rows)
print('Output equivalence and non-multiplexed counters: PASS')
