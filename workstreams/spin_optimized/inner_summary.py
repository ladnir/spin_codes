"""Use isolated-mode probes, never the superseded rotating-buffer diagnostic."""
import csv
import json
from collections import defaultdict
from pathlib import Path
from statistics import median
import sys

root=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).parent/'measurements'
for name,pattern in [('inner_probe_isolated','*.csv'),('forward_expansion','probe-*.csv')]:
    groups=defaultdict(list)
    for path in sorted((root/name).glob(pattern)):
        rows=list(csv.DictReader(path.open()))
        if not rows: continue
        assert len({r['mode'] for r in rows})==1
        groups[(int(rows[0]['m']),rows[0]['mode'])].append(median(float(r['ms']) for r in rows))
    print(name,'m,mode,processes,median_ms,min_process_ms,max_process_ms',sep=': ')
    for key,values in sorted(groups.items()):
        print(','.join(map(str,(*key,len(values),median(values),min(values),max(values)))))
groups=defaultdict(list)
for path in sorted((root/'forward_expansion').glob('m*-r*.json')):
    row=json.loads(path.read_text())
    groups[(row['m'],path.stem.split('-')[1])].append(row['median_ms'])
print('full_encoder: m,variant,processes,median_ms,min_process_ms,max_process_ms')
for key,values in sorted(groups.items()):
    print(','.join(map(str,(*key,len(values),median(values),min(values),max(values)))))
