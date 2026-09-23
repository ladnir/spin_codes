"""Summarize three serial process repetitions; check cross-path map equality."""
import csv
from pathlib import Path
import statistics
import sys

groups = {}
hashes = {}
for directory in sys.argv[1:]:
    for path in sorted(Path(directory).glob('*.csv')):
        label = path.stem.rsplit('-', 1)[0]
        for row in csv.DictReader(path.open()):
            k = row.get('log_k', '18')
            groups.setdefault((label, k), []).append(row)
            key = k, row['rounds']
            checksum = row['checksum']
            if key in hashes and hashes[key] != checksum:
                raise ValueError(f'encoding mismatch across paths: {path}, {key}')
            hashes[key] = checksum

print('variant,log_k,plan_ms,workspace_ms,first_ms,total_ms,warm_ms,setup_bytes,checksum')
for (name, k), rows in sorted(groups.items()):
    if len(rows) != 3:
        raise ValueError(f'expected three processes for {(name, k)}, got {len(rows)}')
    def med(key):
        if key == 'first_ms' and key not in rows[0]:
            key = 'first_encode_ms'
        return statistics.median(float(row[key]) for row in rows)
    values = ','.join(f'{med(key):.6f}' for key in ('plan_ms', 'workspace_ms', 'first_ms', 'total_ms', 'warm_ms'))
    print(f"{name},{k},{values},{rows[0]['setup_bytes']},{rows[0]['checksum']}")
