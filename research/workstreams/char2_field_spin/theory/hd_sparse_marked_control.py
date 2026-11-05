"""Cached FLOAT sparse control preserving input and nonquiet activity.

On a quiet column input weight equals output weight. Mark it exactly by
X^weight; nonquiet input weight is at most a, so X^a is an upper bound.
A second marker Z counts nonquiet columns. First-nine-active expansion
is unchanged. Pay X^(-9ar) and Z^(-dr) once after the ordered word.
Z>1 uses a separate globally bounded bad-setup event, never conditional
independence. This selected-class diagnostic is not a certificate.
"""
from hashlib import sha256
import json
from math import comb, log2
from pathlib import Path
from time import perf_counter

import numpy as np

from hd_bit_syndrome_counted_general import counted_column_matrices, quiet_group_values_general
from hd_bit_syndrome_sparse_gate import FLOORS, word_restriction
from hd_bit_syndrome_sparse_geometry import quiet_data
from hd_inverse_column_control import INPUT, INPUT_SHA, schedule
from hd_post_refresh_quiet_rows import combine_row_bounds, log_matrix_power, moment_atoms
from hd_sparse_activity_control import minkowski_columns

DIRECTORY = Path(__file__).resolve().parent
SETUP = 'hd_bit_syndrome_nonquiet_expurgation_exact.json'
SETUP_SHA = 'e8c8b21914d57817989e9980327965410de34da330b7ff88fe1df16ff89568d2'


def quiet_geometry(rows, measure):
    a = len(rows[0])
    caps, classes, lookup = [], [], {}
    for source in rows:
        values = quiet_data(tuple(int(i) for i in source))['caps']
        if values not in lookup:
            lookup[values] = len(caps)
            caps.append(values)
        classes.append(lookup[values])
    gaps = np.asarray([[row[0], *[row[i]-row[i-1]-1 for i in range(1, a)], 15-row[-1]]
                       for row in rows], dtype=np.int32)
    _, empty = moment_atoms(measure)
    powers = np.stack([log_matrix_power(empty, i) for i in range(17)])
    return caps, np.asarray(classes, dtype=np.int16), gaps, powers


def marked_quiet(geometry, measure, groups, log2_x, *, nonzero):
    caps, classes, gaps, powers = geometry
    terminal = []
    for row in caps:
        logs = [log2(v) if v else -np.inf for v in row]
        if nonzero:
            logs[0] = -np.inf
        terminal.append(groups*quiet_group_values_general(logs, measure, log2_x))
    answer, _ = counted_column_matrices(np.stack(terminal), classes, measure,
        len(caps[0])-1, powers, gaps)
    return answer


def marked_words(branches, quiet, quiet_nz, rows, measure, groups, log2_x, log2_z, floor):
    if min(log2_x, log2_z) < 0 or not np.isfinite([log2_x, log2_z]).all():
        raise ValueError('Finite X,Z>=1 required')
    a = len(rows[0])
    nonquiet = branches['nonquiet']+groups*(a*log2_x+log2_z)
    full = np.minimum(branches['full']+groups*(a*log2_x+log2_z),
                      minkowski_columns(quiet, nonquiet, groups))
    active = minkowski_columns(quiet_nz, nonquiet, groups)
    _, empty = moment_atoms(measure)
    parts = dict(quiet=np.broadcast_to(log_matrix_power(empty, 16), full.shape),
                 nonquiet=active, full=full)
    word = word_restriction(parts, rows, schedule(), 9, groups=groups)
    return (word+(groups-1)*log2(comb(16, 9))
            -groups*(9*a*log2_x+floor*log2_z))


def run(input_path, cache_path, progress=lambda x: None):
    started = perf_counter()
    raw = Path(input_path).read_bytes()
    saved = json.loads(raw)
    common_raw = (DIRECTORY/INPUT).read_bytes()
    if sha256(common_raw).hexdigest() != INPUT_SHA or saved['input_sha256'] != INPUT_SHA:
        raise ValueError('Frozen common input changed')
    common = json.loads(common_raw)
    for name in ('groups', 'mark', 'source_class', 'source_count', 'common_prefix_log2'):
        if saved[name] != common[name]:
            raise ValueError('Scope mismatch: '+name)
    a, r, prefix = saved['arity'], saved['groups'], common['common_prefix_log2']
    if a not in (5, 6, 7) or r != 18:
        raise ValueError('This bounded diagnostic covers a5..7, r18 only')
    for name, digest in saved['source_pins'].items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Frozen sparse source changed: '+name)
    if sha256(Path(saved['cache']).read_bytes()).hexdigest() != saved['cache_sha256']:
        raise ValueError('Sparse branch cache changed')
    setup_raw = (DIRECTORY/SETUP).read_bytes()
    if sha256(setup_raw).hexdigest() != SETUP_SHA:
        raise ValueError('Setup receipt changed')
    setup = json.loads(setup_raw)
    setup_union = next(x for x in setup['unions'] if x['maximum_arity']==8
                       and x['variant']=='deterministic-four-row-event-removed')
    measure = common['state_measure']
    with np.load(saved['cache'], allow_pickle=False) as data:
        rows = data['rows'].copy()
        branches = {key:data[key].copy() for key in ('full','nonquiet')}
    geometry = quiet_geometry(rows, measure)
    best = np.full(len(rows), np.inf)
    unconditional = best.copy()
    controls = []
    for x in (0., 1., 2.):
        quiet = marked_quiet(geometry, measure, r, x, nonzero=False)
        quiet_nz = marked_quiet(geometry, measure, r, x, nonzero=True)
        for z in (0., 1., 2., 4., 8.):
            word = marked_words(branches, quiet, quiet_nz, rows, measure, r, x, z, FLOORS[a])
            value = prefix-8*a*r*log2(255)+combine_row_bounds(word, r)
            best = np.minimum(best, word)
            if z == 0:
                unconditional = np.minimum(unconditional, word)
            cell = dict(log2_x=x, log2_z=z, first9_log2=value)
            controls.append(cell)
            progress(cell)
    offset = prefix-8*a*r*log2(255)
    with Path(cache_path).open('xb') as handle:
        np.savez_compressed(handle, restricted_words=best, unconditional_words=unconditional, rows=rows)
    return dict(arity=a, groups=r, source_class=common['source_class'], mark=common['mark'],
        construction=saved['construction'], controls=controls,
        restricted_per_rowset_minimum_log2=offset+combine_row_bounds(best,r),
        unconditional_per_rowset_minimum_log2=offset+combine_row_bounds(unconditional,r),
        common_prefix_log2=prefix, source_count=common['source_count'],
        input_activity_floor=9*a, active_column_floor=9, nonquiet_floor=FLOORS[a],
        bad_setup_log2_upper=setup_union['log2_upper'], bad_setup_charge_included=False,
        input_receipt=str(Path(input_path).resolve()), input_sha256=sha256(raw).hexdigest(),
        common_measure_sha256=INPUT_SHA, setup_receipt=SETUP, setup_receipt_sha256=SETUP_SHA,
        source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
        cache=str(Path(cache_path).resolve()), cache_sha256=sha256(Path(cache_path).read_bytes()).hexdigest(),
        elapsed_seconds=perf_counter()-started, directed=False, whole_code_certificate=False,
        native_changed=False, regional_provider_runs=0)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--cache', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.cache.exists():
        parser.error('Refuse to overwrite saved evidence')
    result = run(args.input,args.cache,lambda row:print(json.dumps(row),flush=True))
    with args.output.open('x',encoding='utf-8') as handle:
        json.dump(result,handle,indent=2,allow_nan=False)
        handle.write('\n')
    print(json.dumps({k:result[k] for k in ('arity','restricted_per_rowset_minimum_log2',
        'unconditional_per_rowset_minimum_log2','elapsed_seconds')}),flush=True)
