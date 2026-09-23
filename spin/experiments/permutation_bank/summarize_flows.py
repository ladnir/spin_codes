"""Summarize five serial full-encoder runs; checksums must repeat per flow."""
import csv
from pathlib import Path
import statistics as st
import sys
from collections import defaultdict

groups = defaultdict(list)
for path in Path(sys.argv[1]).glob('*.csv'):
    with path.open() as source:
        for row in csv.DictReader(source):
            groups[row['flow'], row['bank_size']].append(row)
fields = ['bank_once_ms', 'code_and_workspace_once_ms', 'workspace_once_ms',
          'bank_bytes', 'plan_bytes', 'workspace_bytes', 'fresh_setup_ms', 'encode_ms', 'total_ms']
print('flow,bank_size,' + ','.join(fields) + ',min_total_ms,max_total_ms')
for key, rows in sorted(groups.items()):
    if len(rows) != 5 or len({r['checksum'] for r in rows}) != 1:
        raise ValueError(f'expected five matching runs: {key}')
    totals = [float(row['total_ms']) for row in rows]
    print(*key, *(f'{st.median(float(row[f]) for row in rows):.6f}' for f in fields),
          f'{min(totals):.6f}', f'{max(totals):.6f}', sep=',')
