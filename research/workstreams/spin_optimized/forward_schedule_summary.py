"""Single-stream schedule comparison; raw samples remain untracked."""
import json
from collections import defaultdict
from pathlib import Path
from statistics import median
import re
import sys

p=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).parent/'measurements/forward_schedule'
groups=defaultdict(list)
for path in sorted(p.glob('m*-r*.json')):
    row=json.loads(path.read_text())
    groups[(row['m'],path.stem.split('-')[1])].append(row['median_ms'])
print('m,schedule,processes,median_ms,min_process_ms,max_process_ms')
for key,values in sorted(groups.items()):
    print(','.join(map(str,(*key,len(values),median(values),min(values),max(values)))))
print('schedule,static_instruction_lines,static_stack_reference_lines')
for path in sorted(p.glob('*-assembly.txt')):
    lines=[line for line in path.read_text().splitlines() if re.match(r'^\s*[0-9a-f]+:\s',line) and len(line.split('\t'))>=3]
    # A rough static code-size/spill indicator, not retired instruction counts.
    print(path.stem.removesuffix('-assembly'),len(lines),sum('(%rsp)' in line or '(%rbp)' in line for line in lines),sep=',')
print('confirmation,m,direction,processes,median_ms,min_process_ms,max_process_ms')
confirmed=defaultdict(list)
for path in sorted(p.glob('confirm-m*-r*.json')):
    row=json.loads(path.read_text())
    confirmed[(row['m'],row['direction'])].append(row['median_ms'])
for key,values in sorted(confirmed.items()):
    print(','.join(map(str,('confirmation',*key,len(values),median(values),min(values),max(values)))))
