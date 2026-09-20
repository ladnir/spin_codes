"""Summarize serial forward A/B measurements; all timings are milliseconds."""
import csv
from collections import defaultdict
from pathlib import Path
from statistics import median
import sys

root = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).parent / 'measurements'
for experiment in ('forward_four', 'forward_prefetch'):
    samples = defaultdict(list)
    for path in sorted((root / experiment).glob('*.csv')):
        per_process = defaultdict(list)
        for row in csv.DictReader(path.open()):
            key = (path.name.split('-')[0], int(row['m']), int(row['lanes']), row['input_layout'])
            per_process[key].append((float(row['encode_ms']), float(row['total_ms'])))
        for key, rows in per_process.items():
            samples[key].append((median(x[0] for x in rows), median(x[1] for x in rows)))
    print(experiment)
    print('variant,m,lanes,layout,processes,encode_16_planes_ms,total_16_planes_ms,encode_per_128_plane_ms,encode_process_min_ms,encode_process_max_ms')
    for key, rows in sorted(samples.items()):
        enc = median(x[0] for x in rows)
        print(','.join(map(str, (*key, len(rows), enc, median(x[1] for x in rows), enc / 16,
                                  min(x[0] for x in rows), max(x[0] for x in rows)))))
