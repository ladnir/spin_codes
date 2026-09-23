"""Summarize three serial runs and check scalar/batched address equality."""
import csv
from pathlib import Path
import statistics as st
import sys
from collections import defaultdict

groups, hashes = defaultdict(list), {}
for path in Path(sys.argv[1]).glob('*.csv'):
    with path.open() as source:
        for row in csv.DictReader(source):
            family = row['family']
            if family in hashes and hashes[family] != row['checksum']:
                raise ValueError(f'address checksum mismatch: {path}')
            hashes[family] = row['checksum']
            groups[family, row['batch']].append(row)
fields = ['bank_ms', 'bank_bytes', 'instance_bytes', 'setup_ms', 'first_route_ms', 'warm_route_ms']
print('family,batch,' + ','.join(fields))
for key, rows in sorted(groups.items()):
    if len(rows) != 3:
        raise ValueError(f'expected three runs: {key}')
    print(*key, *(f'{st.median(float(row[f]) for row in rows):.6f}' for f in fields), sep=',')
