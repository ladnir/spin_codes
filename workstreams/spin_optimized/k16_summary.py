"""Check matched K16 integration timings and report medians of process medians."""
import json
from pathlib import Path
from statistics import median
import sys

root = Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).parent/'measurements'
rows = {}
for seed in (1, 17):
    for kind in ('spin_benchmark', 'spin_bidirectional_benchmark'):
        for cfg in ('12819', '6412'):
            records = [json.loads((root/f'k16-integrated-{kind}-{cfg}-s{seed}-r{r}.json').read_text()) for r in (1, 2, 3)]
            for record in records:
                assert record['m']==16 and record['route_seed']==seed and record['backend']=='avx512'
                assert len(record['samples_ms'])==101
                assert record['workspace_bytes']==3*1024*1024
            assert len({r['output_hash'] for r in records})==1
            assert len({r['retained_setup_bytes'] for r in records})==1
            rows[seed, kind, cfg] = records
        old = median(r['median_ms'] for r in rows[seed, kind, '12819'])
        new = median(r['median_ms'] for r in rows[seed, kind, '6412'])
        print(f'{kind} seed={seed}: old={old:.6f} ms new={new:.6f} ms reduction={100*(1-new/old):.2f}%')
    for cfg in ('12819', '6412'):
        a = rows[seed, 'spin_benchmark', cfg][0]
        b = rows[seed, 'spin_bidirectional_benchmark', cfg][0]
        assert a['output_hash']==b['output_hash'], 'transpose-only/bidirectional linear maps differ'
        print(f'  {cfg}: matching hash={a["output_hash"]}; setup={a["retained_setup_bytes"]}/{b["retained_setup_bytes"]} bytes')
