"""Check final selected-default results and, optionally, regenerated source hashes."""
from pathlib import Path
import hashlib
import json
from statistics import median
import sys
result=Path(__file__).resolve().parent/'measurements'
for m in (16,18,20):
 for seed in ((1,17) if m==16 else (1,)):
  for kind in ('spin_benchmark','spin_bidirectional_benchmark'):
   groups={mode:[json.loads((result/f'production-{mode}-{kind}-m{m}-s{seed}-r{r}.json').read_text()) for r in (1,2,3)] for mode in ('control','default')}
   records=sum(groups.values(),[])
   assert len({d['output_hash'] for d in records})==1
   for d in records:
    assert d['m']==m and d['route_seed']==seed and d['backend']=='avx512'
    assert len(d['samples_ms'])==101 and abs(median(d['samples_ms'])-d['median_ms'])<1e-8
   if m!=16: assert len({d['retained_setup_bytes'] for d in records})==1
   old,new=[median(d['median_ms'] for d in groups[mode]) for mode in ('control','default')]
   print(f'{kind} m={m} seed={seed}: {old:.6f} -> {new:.6f} ms, reduction {100*(1-new/old):.2f}%')
for mode in ('control','default','off','masked','sanitize'):
 assert '100% tests passed' in (result/f'production-{mode}-tests.log').read_text()
assert (result/'production-isa.log').read_text().count('PASS')==9
if len(sys.argv)>1:
 local=Path(sys.argv[1])
 for row in (result/'production-sources.sha256').read_text().splitlines():
  expected,path=row.split(maxsplit=1)
  file=local/path.split('/',1)[1]
  assert hashlib.sha256(file.read_bytes()).hexdigest()==expected,file
 print('Regenerated local sources match the benchmarked sources.')
print('Outputs, medians, four release builds, sanitizers, and ISA isolation checked.')
