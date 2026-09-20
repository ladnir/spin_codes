"""Validate the K16 hand-tuning comparison without choosing by a single sample."""
from pathlib import Path
from statistics import median
import json
root=Path(__file__).resolve().parent/'measurements'
for kind in ('spin_benchmark','spin_bidirectional_benchmark'):
 for seed in (1,17):
  sets={mode:[json.loads((root/f'small-{mode}-{kind}-s{seed}-r{r}.json').read_text()) for r in (1,2,3)] for mode in range(5)}
  hashes=set()
  for mode,records in sets.items():
   for d in records:
    assert d['m']==16 and d['route_seed']==seed and d['backend']=='avx512'
    assert len(d['samples_ms'])==101 and abs(median(d['samples_ms'])-d['median_ms'])<1e-8
    hashes.add(d['output_hash'])
  assert len(hashes)==1
  control=median(d['median_ms'] for d in sets[0])
  for mode, records in sets.items():
   latency=median(d['median_ms'] for d in records)
   print(f'{kind} seed={seed} mode={mode}: {latency:.6f} ms ({100*(1-latency/control):.2f}% reduction), setup={records[0]["retained_setup_bytes"]}')
for m in (18,20):
 for kind in ('spin_benchmark','spin_bidirectional_benchmark'):
  allrows=[]
  for mode in (0,1):
   records=[json.loads((root/f'large-{mode}-{kind}-m{m}-r{r}.json').read_text()) for r in (1,2,3)]
   allrows+=records
   print(f'{kind} m={m} mode={mode}: {median(d["median_ms"] for d in records):.6f} ms')
  assert len({d['output_hash'] for d in allrows})==1
  assert len({d['retained_setup_bytes'] for d in allrows})==1
for mode in range(5):
 assert '100% tests passed' in (root/f'final-small-{mode}-tests.log').read_text()
print('All comparison hashes, medians, large-size storage checks, and correctness suites passed.')
