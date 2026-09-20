"""Inspect sampled assembly; sample locations are not causal stall attribution."""
from pathlib import Path
import re
import sys
groups={};current=None
for line in Path(sys.argv[1]).read_text().splitlines():
    header=re.search(r'[0-9a-f]+ <(.+)>:',line)
    if header:
        current=header[1];groups.setdefault(current,[]);continue
    match=re.match(r'\s*([\d.]+)\s+:\s+([0-9a-f]+):\s+(.*)',line)
    if match and current:groups[current].append((float(match[1]),match[2],match[3]))
for name,rows in groups.items():
    if not ('runFour<true>' in name or 'part0(' in name):continue
    print(name)
    stack=[r for r in rows if re.search(r'\([^)]*%rsp[^)]*\)',r[2])]
    print('Static instruction lines:',len(rows),'stack-address lines:',len(stack),
          'sample percent on stack-address instructions:',round(sum(r[0] for r in stack),2))
    print('Largest sampled instruction locations:')
    for row in sorted(rows,reverse=True)[:12]:print(row)
