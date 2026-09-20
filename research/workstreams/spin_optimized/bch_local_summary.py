"""Compare exact-map BCH schedules using full-encoder timings and assembly."""
from pathlib import Path
import json
import re
import statistics
import sys
p=Path(sys.argv[1]);groups={};hashes={}
for file in sorted(p.glob('*.json')):
    r=json.loads(file.read_text())
    key=(r['m'],r['route_seed'],len(r['samples_ms']),r['configuration'].replace('imt_','').replace('_weight5_seed0_r1','').replace('_weight5',''))
    hashes.setdefault(key,set()).add(r['output_hash'])
    groups.setdefault(file.stem.rsplit('-r',1)[0],[]).append(r['median_ms'])
assert all(len(v)==1 for v in hashes.values()),hashes
for group,values in groups.items():print(group,'median_ms',statistics.median(values),'process_medians',values)
for path in sorted(p.glob('assembly-*.txt')):
    # Count instructions, not objdump's continuation bytes.
    rows=re.findall(r'^\s*[0-9a-f]+:\s+(?:[0-9a-f]{2}\s+)+\s*([a-z][a-z0-9]*\s+.*)$',path.read_text(),re.M)
    stack=[r for r in rows if re.search(r'\([^)]*%rsp[^)]*\)',r)]
    frames=[int(x,16) for x in re.findall(r'sub\s+\$0x([0-9a-f]+),%rsp',path.read_text())]
    print(path.stem,'instructions',len(rows),'stack_reference_instructions',len(stack),'largest_stack_sub',max(frames,default=0))
print('Matched output hashes: PASS')
