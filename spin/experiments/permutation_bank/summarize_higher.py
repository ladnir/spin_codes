"""Aggregate per-bank maxima; preserve the matched control's selection noise."""
import csv
import statistics as st
import sys
from collections import defaultdict

four, adaptive = defaultdict(list), defaultdict(list)
with open(sys.argv[1]) as source:
    for row in list(csv.reader(source))[1:]:
        if not row:
            continue
        key = int(row[1]), int(row[2]), row[4]
        if row[0] == 'four':
            four[key + (int(row[3]),)].append(list(map(float, row[7:10])))
        else:
            adaptive[key].append(list(map(float, row[7:11])))
groups = defaultdict(list)
for key, values in four.items():
    groups[key[:3]].append([max(v[j] for v in values) for j in range(3)] +
                           [st.mean(v[j] for v in values) for j in range(3)])
print('FOUR: bits,bank_size,family,banks,mean_max_xor_zero,mean_max_add_zero,mean_max_one_quarter,mean_xor_zero,mean_add_zero,mean_one_quarter')
for key, values in sorted(groups.items()):
    print(*key, len(values), *(f'{st.mean(v[j] for v in values):.8f}' for j in range(6)), sep=',')
print('ADAPTIVE: bits,bank_size,family,banks,train_pair_rate,holdout_pair_rate,holdout_peak_ge16,holdout_mean_peak')
for key, values in sorted(adaptive.items()):
    print(*key, len(values), *(f'{st.mean(v[j] for v in values):.8f}' for j in range(4)), sep=',')
