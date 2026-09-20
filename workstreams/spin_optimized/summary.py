"""Check all integration comparisons, then print medians across three processes."""
from pathlib import Path
import hashlib
import json
from statistics import median

HERE = Path(__file__).resolve().parent
RESULT = HERE / 'measurements'
for prefix in ('','tuned-'):
 for m in (16,18,20):
    rows = {name:[json.loads((RESULT / f'{prefix}{name}-m{m}-r{r}.json').read_text())
                 for r in (1,2,3)] for name in ('avx2','auto')}
    for name, records in rows.items():
        for d in records:
            assert d['m']==m and d['route_seed']==1
            assert d['backend']==('avx512' if name=='auto' else 'avx2')
            assert len(d['samples_ms'])==101
            assert abs(median(d['samples_ms'])-d['median_ms'])<1e-8
    combined=rows['avx2']+rows['auto']
    for field in ('output_hash','retained_setup_bytes','workspace_bytes'):
        assert len({d[field] for d in combined})==1, (m,field)
    old,new=[median(d['median_ms'] for d in rows[name]) for name in ('avx2','auto')]
    print(f'{prefix or "generic-"}K=2^{m}: AVX2={old:.6f} ms, auto/AVX512={new:.6f} ms, reduction={100*(1-new/old):.2f}%')
for mode in ('on','off','masked','tuned'):
    assert '100% tests passed' in (RESULT / f'test-{mode}.log').read_text()
assert '100% tests passed' in (RESULT / 'sanitizer.log').read_text()
local = RESULT / 'generated-check'
if local.exists():
    for line in (RESULT / 'generated.sha256').read_text().splitlines():
        digest, remote = line.split(maxsplit=1)
        path = local / remote.split('/generated/',1)[1]
        assert hashlib.sha256(path.read_bytes()).hexdigest()==digest, path
    print('Local generated sources match measured sources.')
print('All output hashes, storage sizes, sample medians, and release test logs agree.')
