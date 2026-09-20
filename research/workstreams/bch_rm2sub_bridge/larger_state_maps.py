"""Isolated, write-once larger-state snapshots and independent exact map audits."""
import argparse
import math
from collections import Counter
from pathlib import Path
import bridge as base

SOURCE = Path('C:/Users/peter/.codex/worktrees/3061/permute_conv/workstreams/finite_asymptotic_theory/small_k_replay/rm2sub_calibration_constituents')
DIRECTORY = base.HERE / 'generated/larger_state_inputs_v1'
NAMES = ('t128_s19', 't64_s20')
SUFFIXES = ('selection', 'a_spectrum', 'b_kernel_spectrum')


def snapshot():
    assert not DIRECTORY.exists(), 'Snapshot already exists; do not overwrite'
    payloads = {f'{name}_{suffix}.json': base.read(SOURCE / f'{name}_{suffix}.json')
                for name in NAMES for suffix in SUFFIXES}
    hashes = {name: base.sha(SOURCE / name) for name in payloads}
    for name, payload in payloads.items():
        base.write_new(DIRECTORY / name, payload)
    base.write_new(DIRECTORY / 'manifest.json', dict(source_directory=str(SOURCE),
        source_sha256=hashes, snapshot_sha256={name: base.sha(DIRECTORY / name) for name in payloads}))


def load(name):
    assert name in NAMES
    manifest = base.read(DIRECTORY / 'manifest.json')
    for filename, digest in manifest['snapshot_sha256'].items():
        assert base.sha(DIRECTORY / filename) == digest
    selected = base.read(DIRECTORY / f'{name}_selection.json')
    t, s = (selected['parameters'][key] for key in ('step_bits', 'state_bits'))
    assert name == f't{t}_s{s}' and base.ROWS % t == 0
    generators = [int(v, 16) for v in selected['selected']['A_generator_words_hex']]
    columns = [int(v, 16) for v in selected['selected']['B_columns_hex']]
    assert len(generators) == s and len(columns) == t
    assert len(set(columns)) == t and all(0 < c < 1 << s for c in columns)
    assert generators == [sum(((c >> j) & 1) << i for i, c in enumerate(columns)) for j in range(s)]
    assert all((a & b).bit_count() % 2 == 0 for a in generators for b in generators)
    pivots = {}
    for word in generators:
        assert 0 < word < 1 << t
        while word:
            pivot = word.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = word
                break
            word ^= pivots[pivot]
    assert len(pivots) == s
    spectrum = Counter({0: 1})
    word = 0
    for index in range(1, 1 << s):
        word ^= generators[(index & -index).bit_length() - 1]
        spectrum[word.bit_count()] += 1
    claimed = {row['weight']: row['count'] for row in base.read(DIRECTORY / f'{name}_a_spectrum.json')['spectrum'] if row['count']}
    assert dict(spectrum) == claimed
    kernel = {row['total_weight']: row['kernel_words'] for row in base.read(DIRECTORY / f'{name}_b_kernel_spectrum.json')['by_total_weight']}
    for degree in range(t + 1):
        total = sum(count * sum((-1)**j * math.comb(w, j) * math.comb(t-w, degree-j)
                    for j in range(max(0, degree-t+w), min(degree, w)+1)) for w, count in spectrum.items())
        value, remainder = divmod(total, 1 << s)
        assert remainder == 0 and value == kernel.get(degree, 0) and value >= 0
    assert sum(kernel.values()) == 1 << (t-s)
    return t, s, {w: c for w, c in sorted(spectrum.items()) if w}, kernel


def sources():
    paths = [Path(__file__), Path(base.__file__), DIRECTORY / 'manifest.json']
    paths += [DIRECTORY / f'{name}_{suffix}.json' for name in NAMES for suffix in SUFFIXES]
    return {str(p.relative_to(base.HERE)): base.sha(p) for p in paths}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--snapshot', action='store_true')
    args = parser.parse_args()
    if args.snapshot:
        snapshot()
    rows = []
    for name in NAMES:
        t, s, spectrum, kernel = load(name)
        row = dict(configuration=name, minimum_A_distance=min(spectrum),
                   minimum_kernel_distance=min(w for w, c in kernel.items() if w and c),
                   A_words=1 << s, kernel_words=1 << (t-s), exact_audit_passed=True)
        rows.append(row)
        print(row, flush=True)
    if args.snapshot:
        base.write_new(DIRECTORY / 'independent_audit.json', dict(rows=rows, local_sha256=sources()))
