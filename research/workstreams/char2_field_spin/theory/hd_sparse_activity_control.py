"""Cached FLOAT repair: preserve the RS nine-active-column requirement.

Compute quiet-NZ moments by deleting the quiet polynomial's exact constant
before its group power. Entrywise Minkowski then bounds the active column
kernel quietNZ+nonquiet, without subtracting any upper matrices. First-nine
active columns are expanded before the group power. The empty branch is
the constant functional S(0)^16 under the SAME common measure.
"""
from hashlib import sha256
import json
from math import comb, log2
from pathlib import Path
from time import perf_counter

import numpy as np

from hd_bit_syndrome_counted_general import counted_column_matrices
from hd_bit_syndrome_sparse_gate import terminals, word_restriction
from hd_bit_syndrome_sparse_geometry import quiet_data
from hd_inverse_column_control import INPUT, INPUT_SHA, schedule
from hd_post_refresh_quiet_rows import combine_row_bounds, log_matrix_power, moment_atoms

DIRECTORY = Path(__file__).resolve().parent


def minkowski_columns(left, right, groups):
    if type(groups) is not int or groups < 1 or left.shape != right.shape:
        raise ValueError('Matched moment matrices and positive group count required')
    return groups*np.logaddexp2(left/groups, right/groups)


def quiet_nonzero_columns(rows, measure, groups):
    a = len(rows[0])
    caps, classes, index = [], [], {}
    for source in rows:
        values = quiet_data(tuple(int(i) for i in source))['caps']
        values = (0,)+values[1:]
        cls = index.get(values)
        if cls is None:
            cls = len(caps)
            index[values] = cls
            caps.append(values)
        classes.append(cls)
    _, empty = moment_atoms(measure)
    gaps = np.asarray([[row[0], *[row[i]-row[i-1]-1 for i in range(1, a)], 15-row[-1]]
                       for row in rows], dtype=np.int32)
    powers = np.stack([log_matrix_power(empty, i) for i in range(17)])
    answer, _ = counted_column_matrices(terminals(caps, measure, groups),
        np.asarray(classes, dtype=np.int16), measure, a, powers, gaps)
    return answer


def controls(branches, quiet_nonzero, rows, measure, groups, arity, prefix):
    full = np.minimum(branches['full'], minkowski_columns(branches['quiet'], branches['nonquiet'], groups))
    active = minkowski_columns(quiet_nonzero, branches['nonquiet'], groups)
    _, empty = moment_atoms(measure)
    empty = np.broadcast_to(log_matrix_power(empty, 16), full.shape)
    active_branches = dict(quiet=empty, nonquiet=active, full=full)
    variants = {}
    for one_copy in (False, True):
        word = word_restriction(active_branches, rows, schedule(), 9,
                                one_copy=one_copy, groups=groups)
        if not one_copy:
            word += (groups-1)*log2(comb(16, 9))
        variants['one_constrained_copy' if one_copy else 'first9_uniform'] = word
    variants['per_rowset_minimum'] = np.minimum(*variants.values())
    offset = prefix-8*arity*groups*log2(255)
    return ({name: offset+combine_row_bounds(word, groups) for name, word in variants.items()},
            dict(full=full, active=active, quiet_nonzero=quiet_nonzero, empty=empty))


def run(input_path, cache_path):
    started = perf_counter()
    raw = Path(input_path).read_bytes()
    saved = json.loads(raw)
    if saved['arity'] not in range(2, 8) or saved['groups'] != 18 or saved['input_sha256'] != INPUT_SHA:
        raise ValueError('The frozen sparse gate scope is required')
    for name, digest in saved['source_pins'].items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Sparse source changed: '+name)
    common = (DIRECTORY/INPUT).read_bytes()
    if sha256(common).hexdigest() != INPUT_SHA:
        raise ValueError('Common regional receipt changed')
    measure = json.loads(common)['state_measure']
    branch_path = Path(saved['cache'])
    if sha256(branch_path.read_bytes()).hexdigest() != saved['cache_sha256']:
        raise ValueError('Sparse branch cache changed')
    with np.load(branch_path, allow_pickle=False) as data:
        rows = data['rows'].copy()
        branches = {name: data[name].copy() for name in ('full', 'quiet', 'nonquiet')}
    quiet_nz = quiet_nonzero_columns(rows, measure, saved['groups'])
    values, matrices = controls(branches, quiet_nz, rows, measure, saved['groups'],
                                 saved['arity'], saved['common_prefix_log2'])
    with Path(cache_path).open('xb') as handle:
        np.savez_compressed(handle, **matrices, rows=rows)
    return dict(arity=saved['arity'], groups=saved['groups'], source_class=saved['source_class'],
        construction=saved['construction'], mark=saved['mark'], controls=values,
        input_receipt=str(Path(input_path).resolve()), input_sha256=sha256(raw).hexdigest(),
        common_measure_sha256=INPUT_SHA, cache=str(Path(cache_path).resolve()),
        cache_sha256=sha256(Path(cache_path).read_bytes()).hexdigest(),
        source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
        active_column_floor=9, convexity_cost_bits=(saved['groups']-1)*log2(comb(16,9)),
        separate_nonquiet_floor_enforced=False, bad_setup_charge_required=False,
        entrywise_column_minkowski=True, upper_bound_subtractions=0,
        regional_provider_runs=0, directed=False, whole_code_certificate=False,
        elapsed_seconds=perf_counter()-started)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--cache', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.cache.exists():
        parser.error('Refuse to overwrite saved evidence')
    result = run(args.input, args.cache)
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps(result, indent=2), flush=True)
