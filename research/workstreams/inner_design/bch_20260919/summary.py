"""Validate local timing receipts and report exact-map BCH optimization results."""
from pathlib import Path
import hashlib
import json
import statistics as st

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def main():
    checked = 0
    for line in (HERE/'measurements/sources-final.sha256').read_text().splitlines():
        digest, name = line.split(maxsplit=1)
        if Path(name).suffix not in ('.h','.cpp'): continue
        path = ROOT/name
        assert path.is_file(), path
        assert hashlib.sha256(path.read_bytes()).hexdigest() == digest, path
        checked += 1
    print(f'{checked} measured C++ source/header hashes match.')
    records, hashes = {}, {}
    for path in (HERE/'measurements').glob('*.jsonl'):
        d = json.loads(path.read_text())
        samples = d['samples_ms']
        assert abs(st.median(samples)-d['median_ms']) < 1e-6
        if 'configuration' in d:
            assert d['inplace'] and d['outer_length']==256 and d['outer_dimension']==128
            assert len(samples)==d['trials']
            key = (d['m'],d['route_seed'],d['trials'])
            assert hashes.setdefault(key,d['output_hash'])==d['output_hash'], path
        records[path.stem]=d
    print('Full-output hashes agree across equal seed/size/iteration count.')
    for m, seed in [(16,1),(18,1),(20,1),(20,17)]:
        baseline = st.median(records[f'final-control-m{m}-s{seed}-r{r}']['median_ms'] for r in (1,2,3))
        print(f'\nK=2^{m}, route seed {seed}')
        for name in (['control','packedshare4'] if m!=20 else ['control','packedfour','packedshare3','packedshare4']):
            times = [records[f'final-{name}-m{m}-s{seed}-r{r}']['median_ms'] for r in (1,2,3)]
            median = st.median(times)
            print(f'{name:14s} {median:.6f} ms  range {min(times):.6f}..{max(times):.6f}  reduction {(1-median/baseline)*100:.2f}%')
    print('\nSynthetic BCH-only throughput, 8192 rows per call:')
    for mode in ('hot','stream'):
        for name in ('control','packedshare4'):
            times = [records[f'outer-{name}-{mode}-r{r}']['median_ms'] for r in (1,2,3)]
            print(mode,name,st.median(times),min(times),max(times))
    print('\nRelevant compiler stack-usage records:')
    for line in (HERE/'measurements/stack-usage-final.txt').read_text().splitlines():
        if any(key in line for key in ('generated/BchCircuit.cpp','generated/dfs.cpp','generated/packedfour.cpp','generated/packedshare4.cpp')):
            print(line)


if __name__=='__main__': main()
