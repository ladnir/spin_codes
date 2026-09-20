"""Medians and process-to-process ranges for the serial K16/K20 wide sweep."""
from collections import defaultdict
import csv
from pathlib import Path
from statistics import median
import sys

root = Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).parent/'measurements/wide_tuning'
for pattern in ('k16-*.csv', 'screen-*.csv', 'confirm-*.csv'):
    samples, runs = defaultdict(list), defaultdict(list)
    for path in sorted(root.glob(pattern)):
        local = defaultdict(list)
        with path.open() as f:
            for row in csv.DictReader(f):
                key = (row['configuration'], int(row['tile_rows']), int(row['lanes']), row['input_layout'])
                samples[key].append(row)
                local[key].append(float(row['total_ms']))
        for key, times in local.items():
            runs[key].append(median(times))
    if not samples:
        continue
    print('\n' + pattern)
    print('| configuration | tile rows | lanes | input | encode ms | columns ms | total ms | run range ms | scratch MiB |')
    print('|:---|---:|---:|:---|---:|---:|---:|:---|---:|')
    for key, data in sorted(samples.items()):
        cfg, tile, lanes, layout = key
        fields = [median(float(r[f]) for r in data) for f in ('encode_ms','column_ms','total_ms')]
        times = runs[key]
        print(f'| {cfg} | {tile} | {lanes} | {layout} | ' + ' | '.join(f'{v:.3f}' for v in fields) +
              f' | {min(times):.3f}-{max(times):.3f} | {int(data[0]["workspace_bytes"])/2**20:.2f} |')
