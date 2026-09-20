"""Summarize matched one-stream measurements and separate diagnostic profiles."""
import json
from collections import defaultdict
from pathlib import Path
from statistics import median
import re
import sys

root = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).parent/'measurements/single_forward'
groups = defaultdict(list)
for path in sorted(root.glob('m*-r*.json')):
    row = json.loads(path.read_text())
    # Retain the requested default-vs-explicit distinction even if tiles agree.
    key = (row['m'], row['direction'], path.stem.split('-')[2], row['tile_rows'])
    groups[key].append(row['median_ms'])
print('m,direction,requested_tile,actual_tile,processes,median_ms,min_process_ms,max_process_ms')
for key, values in sorted(groups.items()):
    print(','.join(map(str, (*key,len(values),median(values),min(values),max(values)))))
print('profile,m,tile,outer_ms,route_ms,inner_gather_ms')
for path in sorted(root.glob('profile-m*.log')):
    match = re.search(r'outer_ms=([\d.]+) route_ms=([\d.]+) inner_gather_ms=([\d.]+)',path.read_text())
    if not match:
        raise RuntimeError(f'missing profile in {path}')
    row = json.loads(path.with_suffix('.json').read_text())
    print(','.join(map(str, ('profile',row['m'],row['tile_rows'],*match.groups()))))
print('direct_comparison,m,variant,processes,median_ms,min_process_ms,max_process_ms')
direct=defaultdict(list)
for path in sorted((root.parent/'single_direct').glob('m*-r*.json')):
    row=json.loads(path.read_text())
    direct[(row['m'],path.stem.split('-')[1])].append(row['median_ms'])
for key,values in sorted(direct.items()):
    print(','.join(map(str, ('direct_comparison',*key,len(values),median(values),min(values),max(values)))))
