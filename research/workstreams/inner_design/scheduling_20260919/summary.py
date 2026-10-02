"""Validate timing receipts and summarize this bounded scheduling exploration."""
import json
import statistics as st
import hashlib
from pathlib import Path

HERE = Path(__file__).resolve().parent


def main():
    root = HERE.parents[2]
    checked = 0
    for line in (HERE / 'measurements/sources-final.sha256').read_text().splitlines():
        digest, name = line.split(maxsplit=1)
        if Path(name).suffix not in ('.cpp', '.h'):
            continue  # Test targets and orchestration were extended after timing.
        target = root / name
        assert target.is_file(), target
        assert hashlib.sha256(target.read_bytes()).hexdigest() == digest, target
        checked += 1
    print(f'Final confirmation source check: {checked} C++ source/header hashes match.')
    records = {}
    hashes = {}
    for path in sorted((HERE / 'measurements').glob('*.jsonl')):
        d = json.loads(path.read_text())
        assert d['inplace'] and d['m'] == 20
        assert d['outer_length'] == 256 and d['outer_dimension'] == 128
        samples = d['samples_ms']
        assert len(samples) == d['trials']
        assert abs(st.median(samples)-d['median_ms']) < 1e-6
        key = (d['trials'], d['route_seed'])
        assert hashes.setdefault(key, d['output_hash']) == d['output_hash'], path
        records[path.stem] = d
    print('All complete-output hashes agree for equal input/seed/call count.')
    for seed in (1, 17):
        print(f'\nSeed {seed}: medians of three 101-call process medians')
        baseline = st.median(records[f'final-control-s{seed}-r{r}']['median_ms'] for r in (1, 2, 3))
        for name in ('control', 'pf16', 'tile1024'):
            times = [records[f'final-{name}-s{seed}-r{r}']['median_ms'] for r in (1, 2, 3)]
            median = st.median(times)
            print(f'{name:12s} {median:.6f} ms  range {min(times):.6f}..{max(times):.6f}  change {(median/baseline-1)*100:+.2f}%')
    print('\nPilot/refinement screens (single processes, not promotion evidence):')
    for name, d in records.items():
        if name.startswith(('pilot-', 'refine-', 'layout-')):
            print(f'{name:28s} {d["median_ms"]:.6f} ms')
    print('\nInstrumented means over warmups and measured calls:')
    for path in sorted((HERE / 'measurements').glob('profile-t*.log')):
        print(path.stem, path.read_text().strip())


if __name__ == '__main__':
    main()
