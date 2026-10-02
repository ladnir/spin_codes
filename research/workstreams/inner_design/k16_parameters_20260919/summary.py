"""Summarize serial process medians and check circuit-variant equality."""
from pathlib import Path
import json
import statistics
import sys
import hashlib

p=Path(sys.argv[1])
cells=['baseline']+(sys.argv[2:] or ['32_15_0','32_15_1','64_15_0','64_15_1','64_19_0','64_19_1'])
hashes={}
for seed in (1,17):
    baseline=None
    for cell in cells:
        records=[json.loads((p/f'perf-{cell}-s{seed}-r{r}.json').read_text()) for r in (1,2,3)]
        assert all(r['m']==16 and r['backend']=='avx512' and len(r['samples_ms'])==101 for r in records)
        assert len({r['output_hash'] for r in records})==1
        assert len({r['retained_setup_bytes'] for r in records})==1
        value=statistics.median(r['median_ms'] for r in records)
        if cell=='baseline': baseline=value
        else:
            key=(seed,cell.rsplit('_',1)[0])
            assert hashes.setdefault(key,records[0]['output_hash'])==records[0]['output_hash']
        print(seed,cell,f'{value:.6f} ms',f'{100*(1-value/baseline):.2f}% reduction',
              'setup',records[0]['retained_setup_bytes'])
for cell in cells[1:]:
    assert 'correctness=PASS' in (p/f'test-{cell.rsplit("_",1)[0]}-{cell[-1]}.log').read_text()
matched=0
for line in (p/'sources.sha256').read_text().splitlines():
    digest,name=line.split(maxsplit=1)
    if name.endswith('/AsymmetricMap.h'):
        _,_,t,s,shared,*batch=name.split('/')[0].split('-')
        header=(p.parent/'subspace' if batch else p.parent)/f'Map{t}S{s}.h'
        assert hashlib.sha256(header.read_bytes()).hexdigest()==digest
        matched+=1
print('Exact local/remote map header matches:',matched)
