"""Summarize diagnostics across fixed banks, not across fresh-bank instances."""
import csv
import statistics
import sys
from collections import defaultdict

groups = defaultdict(list)
for path in sys.argv[1:]:
    with open(path) as source:
        for row in csv.DictReader(source):
            groups[int(row['bits']), int(row['bank_size']), row['family']].append(row)
metrics = ['max_xor_excess', 'max_add_excess', 'top_valuation_excess',
           'exact_bank_valuation_excess', 'valuation_prediction_error',
           'interval_peak_ge16', 'subspace_peak_ge16', 'bank_preimage_peak_ge16',
           'interval_overlap_tail', 'subspace_overlap_tail', 'bank_preimage_overlap_tail']
print('bits,bank_size,family,banks,' + ','.join(metrics))
for key, rows in sorted(groups.items()):
    means = [statistics.mean(float(row[m]) for row in rows) for m in metrics]
    print(','.join(map(str, key)) + f',{len(rows)},' + ','.join(f'{x:.7f}' for x in means))
