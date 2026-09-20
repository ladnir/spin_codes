"""Summarize serial wide-encoder trials; no third-party dependencies."""
import csv
from collections import defaultdict
from pathlib import Path
from statistics import median
import sys

root = Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).parent/'measurements/wide'
rows = defaultdict(list)
per_run = defaultdict(list)
for path in sorted(root.glob('m*-r*.csv')):
    local = defaultdict(list)
    with path.open() as f:
        for row in csv.DictReader(f):
            key = (int(row['m']), int(row['lanes']), row['input_layout'])
            rows[key].append(row)
            local[key].append(float(row['total_ms']))
    for key, times in local.items():
        per_run[key].append(median(times))
if not rows:
    raise SystemExit('No wide benchmark results')
print('| log2 K | bits | input | pack ms | encode ms | columns ms | total ms | speedup | run-median range ms | scratch MiB |')
print('|---:|---:|:---|---:|---:|---:|---:|---:|:---|---:|')
for (m,w,layout), data in sorted(rows.items()):
    fields = [median(float(r[f]) for r in data) for f in ('pack_ms','encode_ms','column_ms','total_ms')]
    baseline = median(float(r['total_ms']) for r in rows[m,1,'native'])
    times = per_run[m,w,layout]
    print(f'| {m} | {128*w} | {layout} | ' + ' | '.join(f'{x:.3f}' for x in fields) +
          f' | {baseline/fields[-1]:.2f}x | {min(times):.3f}–{max(times):.3f} | {int(data[0]["workspace_bytes"])/2**20:.1f} |')
