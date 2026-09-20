"""Summarize serial integration timings and check cross-library outputs."""
import argparse
import json
from pathlib import Path
import statistics

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('directory',type=Path)
a=p.parse_args()
for seed in (1,17):
    for cfg in ('12819','6412','6412r2'):
        hashes=set()
        for kind in ('spin_benchmark','spin_bidirectional_benchmark'):
            rows=[json.loads((a.directory/f'{kind}-{cfg}-s{seed}-r{r}.json').read_text()) for r in (1,2,3)]
            assert all(row['m']==16 and row['route_seed']==seed and row['backend']=='avx512' for row in rows)
            assert len({row['output_hash'] for row in rows})==1
            hashes.update(row['output_hash'] for row in rows)
            print(kind,cfg,'seed',seed,'median_ms',statistics.median(row['median_ms'] for row in rows),
                'setup',rows[0]['retained_setup_bytes'],'workspace',rows[0]['workspace_bytes'])
        assert len(hashes)==1,(seed,cfg,'cross-library output mismatch')
print('All configurations have matching output hashes across both libraries and repetitions.')
