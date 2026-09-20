"""Equal-work wide-forward results: every row encodes 16 128-bit planes."""
from pathlib import Path
from statistics import median
from collections import defaultdict
import csv
import sys

root=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).parent/'measurements'
for name in ('wide_schedule','wide_vl','wide_range','natural_lengths'):
    groups=defaultdict(list)
    for path in sorted((root/name).glob('*.csv')):
        per_process=defaultdict(list)
        for r in csv.DictReader(path.open()):
            key=(path.name.split('-')[0],int(r['m']),int(r['lanes']),r['input_layout'],int(r['tile_rows']))
            per_process[key].append((float(r['encode_ms']),float(r['total_ms'])))
        for key,values in per_process.items():
            groups[key].append((median(v[0] for v in values),median(v[1] for v in values)))
    print(name)
    print('variant,m,lanes,layout,tile,processes,encode_ms,total_ms,encode_min,encode_max')
    for key,values in sorted(groups.items()):
        print(','.join(map(str,(*key,len(values),median(v[0] for v in values),median(v[1] for v in values),min(v[0] for v in values),max(v[0] for v in values)))))
