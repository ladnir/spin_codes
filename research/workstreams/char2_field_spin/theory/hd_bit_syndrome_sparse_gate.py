"""Bounded FLOAT sparse-family test with prescribed nonquiet-column floors.

Every fixed original RS row set is retained through the inverse-column
schedule. Regional moments come from one authenticated r18/F<=1 cache.
Quiet, nonquiet and full branches are separately positive majorants.
No upper-bound subtraction, regional provider, or native change occurs.
The separately bounded bad-setup event is NOT included in these values.
"""
from hashlib import sha256
import json
from math import comb, log2
from pathlib import Path
from time import perf_counter

import numpy as np

from hd_bit_syndrome_counted_general import counted_column_matrices, quiet_group_values_general
from hd_bit_syndrome_gate import ordered_cover_products, weighted_cover_matrices
from hd_bit_syndrome_sparse_geometry import prepare_geometry
from hd_inverse_column_control import INPUT, INPUT_SHA, schedule
from hd_post_refresh_quiet_rows import combine_row_bounds, log_matmul, log_matrix_power, moment_atoms

DIRECTORY = Path(__file__).resolve().parent
FLOORS = {1: 9, 2: 9, 3: 9, 4: 5, 5: 5, 6: 4, 7: 3, 8: 2}


def terminals(caps, measure, groups):
    return np.stack([groups*quiet_group_values_general(
        [log2(v) if v else -np.inf for v in row], measure, 0.) for row in caps])


def column_branches(prepared, measure, groups, *, check_time=lambda: None, progress=lambda x: None):
    a = len(prepared['rowsets'][0])
    bank = weighted_cover_matrices(measure, groups, prepared['gap_values'])
    products = ordered_cover_products(bank['cover_matrices'], prepared['gap_parents'], prepared['gap_labels'])
    branches = {}
    for name, caps in (('full', prepared['outside_caps']), ('nonquiet', prepared['nonquiet_caps'])):
        branch, _ = counted_column_matrices(terminals(caps, measure, groups), prepared['classes'],
            measure, 16-a, products, prepared['gap_codes'], check_time=check_time)
        branches[name] = branch
        progress(name+'_column')
    _, empty = moment_atoms(measure)
    gaps = np.asarray([[row[0], *[row[i]-row[i-1]-1 for i in range(1, a)], 15-row[-1]]
                       for row in prepared['rowsets']], dtype=np.int32)
    powers = np.stack([log_matrix_power(empty, i) for i in range(17)])
    branches['quiet'], _ = counted_column_matrices(
        terminals(prepared['quiet_caps'], measure, groups), prepared['quiet_classes'],
        measure, a, powers, gaps, check_time=check_time)
    progress('quiet_column')
    return branches


def word_restriction(branches, rows, maps, minimum, *, one_copy=False, groups=18):
    """First-d flag DP, before the separately reported convexity factor."""
    if type(minimum) is not int or not 0 <= minimum <= len(maps):
        raise ValueError('An in-scope nonquiet floor is required')
    rows = tuple(tuple(int(i) for i in row) for row in rows)
    lookup = {row: i for i, row in enumerate(rows)}
    full = branches['full']
    quiet, nonquiet = branches['quiet'], branches['nonquiet']
    if one_copy and groups != 1:
        quiet = quiet/groups + full*((groups-1)/groups)
        nonquiet = nonquiet/groups + full*((groups-1)/groups)
    current = np.full((minimum+1, *full.shape), -np.inf)
    current[0, :, 0, 0] = current[0, :, 1, 1] = 0.
    for p in maps:
        if sorted(p) != list(range(16)):
            raise ValueError('A sixteen-coordinate permutation is required')
        index = [lookup[tuple(sorted(p[i] for i in row))] for row in rows]
        following = np.full_like(current, -np.inf)
        following[minimum] = log_matmul(current[minimum], full[index])
        for d in range(minimum):
            following[d] = np.logaddexp2(following[d], log_matmul(current[d], quiet[index]))
            following[d+1] = np.logaddexp2(following[d+1], log_matmul(current[d], nonquiet[index]))
        current = following
    return np.logaddexp2.reduce(current[minimum, :, 0, :], axis=1)


def normalized_controls(branches, rows, groups, a, floor, prefix):
    maps = schedule()
    full = word_restriction(branches, rows, maps, 0, groups=groups)
    restricted = word_restriction(branches, rows, maps, floor, groups=groups)
    cost = (groups-1)*log2(comb(len(maps), floor))
    restricted += cost
    one_copy = word_restriction(branches, rows, maps, floor, one_copy=True, groups=groups)
    best = np.minimum(np.minimum(full, restricted), one_copy)
    offset = prefix-8*a*groups*log2(255)
    values = dict(unrestricted=full, first_d_uniform=restricted, one_constrained_copy=one_copy,
                  per_rowset_minimum=best)
    return {name: offset+combine_row_bounds(word, groups) for name, word in values.items()}, cost


def run(arity, cache_path, *, max_seconds=120, progress=lambda x: None):
    if arity not in range(2, 8) or not 0 < max_seconds <= 240:
        raise ValueError('Arity2..7 and an explicit bounded budget required')
    if Path(cache_path).exists():
        raise ValueError('Do not overwrite a branch cache')
    raw = (DIRECTORY/INPUT).read_bytes()
    if sha256(raw).hexdigest() != INPUT_SHA:
        raise ValueError('Frozen common-measure receipt changed')
    saved = json.loads(raw)
    for name, expected in saved['reviewed_sources'].items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != expected:
            raise ValueError('Frozen dependency changed: '+name)
    started = perf_counter()
    stages = []
    def check_time():
        if perf_counter()-started > max_seconds:
            raise TimeoutError('Sparse contraction exceeded its bounded budget')
    def stage(label):
        row = dict(stage=label, elapsed_seconds=perf_counter()-started)
        stages.append(row)
        progress(row)
    prepared = prepare_geometry(arity, check_time=check_time)
    stage('all_labelled_geometry')
    groups = saved['groups']
    branches = column_branches(prepared, saved['state_measure'], groups,
                               check_time=check_time, progress=stage)
    check_time()
    controls, cost = normalized_controls(branches, prepared['rowsets'], groups, arity,
                                         FLOORS[arity], saved['common_prefix_log2'])
    stage('ordered_scheduled_controls')
    with Path(cache_path).open('xb') as handle:
        np.savez_compressed(handle, **branches, rows=np.asarray(prepared['rowsets']))
    local_sources = ['hd_bit_syndrome_sparse_gate.py', 'hd_bit_syndrome_sparse_geometry.py',
        'hd_bit_syndrome_flat_geometry.py', 'hd_bit_syndrome_counted_general.py',
        'hd_bit_syndrome_arity9_gate.py', 'hd_bit_syndrome_gate.py', 'hd_inverse_column_control.py',
        'hd_post_refresh_quiet_rows.py']
    return dict(construction='G16-Hstar-column-inverse-input-permutations-research', arity=arity,
        groups=groups, source_class=saved['source_class'], source_count=saved['source_count'],
        mark=saved['mark'], common_prefix_log2=saved['common_prefix_log2'],
        nonquiet_floor=FLOORS[arity], uniform_convexity_cost_bits=cost, controls=controls,
        geometry_metadata=prepared['metadata'], cache=str(Path(cache_path).resolve()),
        cache_sha256=sha256(Path(cache_path).read_bytes()).hexdigest(),
        input_receipt=INPUT, input_sha256=INPUT_SHA,
        source_pins={name: sha256((DIRECTORY/name).read_bytes()).hexdigest() for name in local_sources},
        stages=stages, elapsed_seconds=perf_counter()-started,
        separate_global_setup_charge_required=True, bad_setup_charge_included=False,
        input_marker=1, nonzero_column_marker=1, original_rowsets_fixed=True,
        directed=False, whole_code_certificate=False, native_changed=False, regional_provider_runs=0)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--arity', type=int, required=True)
    parser.add_argument('--max-seconds', type=float, default=120)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--cache', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.cache.exists():
        parser.error('Refuse to overwrite saved evidence')
    result = run(args.arity, args.cache, max_seconds=args.max_seconds,
                 progress=lambda x: print(json.dumps(x), flush=True))
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps({key: result[key] for key in ('arity', 'controls', 'elapsed_seconds')}), flush=True)
