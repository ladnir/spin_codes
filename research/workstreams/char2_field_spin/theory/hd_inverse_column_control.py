"""FLOAT proof diagnostic for a new deterministic column-order schedule.

In column j, input RS row i enters mixer position inverse_GF16(i XOR j).
Use polynomial x^4+x+1 and define inverse(0)=0. This changes the outer
construction, but not its RS constituents, column mixer, or regional inner
operators. Reuse fixed-labelled-A column bounds from the frozen cache;
replace Q_A^16 by Q_{P_0 A} ... Q_{P_15 A}. Keep the original source row
set fixed throughout the word. This is not a runtime measurement or a
directed certificate. Only r18,F<=1,mixed arities7/8/9 are considered.
"""
from fractions import Fraction
from hashlib import sha256
import json
from math import log2
from pathlib import Path

import numpy as np

from hd_bit_syndrome_activity_control import constant_columns, marked_columns
from hd_post_refresh_quiet_rows import combine_row_bounds, log_matmul, log_matrix_power

DIRECTORY = Path(__file__).resolve().parent
INPUT = 'hd_mds_exact_empty_mixed789_9927_pinned_float.json'
INPUT_SHA = 'c89ef4973de29007fe293817d2e3b9f582246a45d19ca19c7fcce2554e7d33b8'


def gf16_multiply(a, b):
    answer = 0
    for _ in range(4):
        if b & 1:
            answer ^= a
        b >>= 1
        a <<= 1
        if a & 16:
            a ^= 0x13
    return answer


def schedule():
    inverse = [0] + [next(j for j in range(1, 16) if gf16_multiply(i, j) == 1)
                     for i in range(1, 16)]
    return tuple(tuple(inverse[i ^ column] for i in range(16)) for column in range(16))


def scheduled_words(column, rows, permutations):
    rows = tuple(tuple(int(i) for i in row) for row in rows)
    index = {row: i for i, row in enumerate(rows)}
    if len(index) != len(rows):
        raise ValueError('Distinct sorted row sets required')
    result = np.full(column.shape, -np.inf)
    result[:, 0, 0] = result[:, 1, 1] = 0.
    for permutation in permutations:
        if sorted(permutation) != list(range(16)):
            raise ValueError('Every column map must be a permutation')
        lookup = [index[tuple(sorted(permutation[i] for i in row))] for row in rows]
        result = log_matmul(result, column[lookup])
    return np.logaddexp2.reduce(result[:, 0, :], axis=1)


def run():
    raw = (DIRECTORY/INPUT).read_bytes()
    if sha256(raw).hexdigest() != INPUT_SHA:
        raise ValueError('Frozen input changed')
    saved = json.loads(raw)
    cache = Path(saved['cache'])
    if not cache.is_absolute():
        cache = DIRECTORY/cache
    if sha256(cache.read_bytes()).hexdigest() != saved['cache_sha256']:
        raise ValueError('Column cache changed')
    r = saved['groups']
    prefix = saved['common_prefix_log2']
    constant, _ = constant_columns(saved['state_measure'])
    values = []
    with np.load(cache, allow_pickle=False) as columns:
        for a in (7, 8, 9):
            marker = Fraction(5, 2) if a == 7 else Fraction(1)
            column = marked_columns(columns['column'+str(a)], constant, marker, r)
            rows = columns['rowsets'+str(a)]
            static = np.logaddexp2.reduce(log_matrix_power(column, 16)[:, 0, :], axis=1)
            moved = scheduled_words(column, rows, schedule())
            offset = prefix-8*a*r*log2(255)-9*r*log2(marker)
            values.append(dict(arity=a, marker=str(marker), rowset_count=len(rows),
                static_log2=offset+combine_row_bounds(static, r),
                scheduled_log2=offset+combine_row_bounds(moved, r)))
    mixed = prefix+r*np.logaddexp2.reduce([(v['scheduled_log2']-prefix)/r for v in values])
    return dict(construction='Hstar with inverse-GF16 input schedule by RS column',
        schedule=[list(p) for p in schedule()], input_receipt=INPUT, input_sha256=INPUT_SHA,
        input_cache_sha256=saved['cache_sha256'], groups=r, source_class=saved['source_class'],
        mark=saved['mark'], arities=values, scheduled_mixed_log2=mixed,
        all_labelled_original_rowsets_retained=True, original_inner=True,
        directed=False, whole_code_certificate=False, native_implementation_changed=False,
        runtime_measured=False, regional_provider_runs=0,
        source_sha256=sha256(Path(__file__).read_bytes()).hexdigest())


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Refuse to overwrite saved evidence')
    result = run()
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'schedule'}, indent=2), flush=True)
