"""Bounded exact SHEAR a3 FLOAT gate, with optional frozen a1/a2 replay.

Each adjacent-column pair has six original inputs and32 jointly determined
outputs. Its exact raw (input weight,joint mask) counts come from the80-form
catalog. For fixed original A, form each pair's scalar polynomial entrywise:
sum_(w,S) n(w,S) J_S^(1/18) X^w. Raise its entries to18, then multiply the
eight pair matrices chronologically. Retain total input weight>=486 using
positive reversed-deficit truncation. Divide the word root by255^24 once.
The source/output/projective prefix is applied once after fixed-A root sums.

The same original A is retained across every pair; no independent-column
support law is asserted. Final independent D2 follows the shear. The ideal
D1 comparison allows only one scalar gauge per pair, as in the frozen
a1/a2 controller; two separate column gauges are not silently reused.

Default scope is three prescribed original a3 sets, not all560 sets.
The explicit replay-a12 mode covers the frozen16+120 labelled sets and
compares complete fixed-A roots. All arithmetic after exact integer counts
is positive FLOAT, not a directed certificate. No provider or fitting runs.
"""
from hashlib import sha256
from itertools import combinations
import json
from math import isfinite, log2
from pathlib import Path
from time import perf_counter

import numpy as np

import hd_lambda_shear_categories as catalog_api
import hd_lambda_column_shear_control as frozen
import hd_lambda_double_middle_control as tail_backend


DIRECTORY = Path(__file__).resolve().parent
SOURCES = ((0, 1, 2), (0, 1, 8), (0, 4, 8))
REPLAY = 'hd_lambda_column_shear_control_9927_float.json'
REPLAY_SHA = '2ec123db4909fcf7e09c1e05658e1e76142b79b1829a6e4690416f22a349fda3'
PINS = {
    'hd_lambda_shear_categories.py': '7dfa8954fe9697a96232f49d1c61a5abafe675c3e526aafa7a598ac581b31a9b',
    'hd_lambda_column_shear_control.py': 'a1ede595d659ddd4c20d1b5ea33a98158e7c9484e34a958d2d682f985214de96',
    'hd_lambda_double_middle_control.py': '78655f83f08ae01118beb6d9f6258daaddd7a5035bc5837694b6243d762e9820',
}


def dependencies():
    for name, digest in PINS.items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Frozen dependency changed: '+name)
    return {**catalog_api.dependencies(), **frozen.dependencies(), **tail_backend.dependencies(),
            **PINS, Path(__file__).name: sha256(Path(__file__).read_bytes()).hexdigest()}


def pair_coefficients(weights, masks, counts, table, arity, *, groups=18, check_time=None):
    """Exact counts stay outside norm roots; gather joint moments in bounded chunks."""
    catalog_api._arity(arity)
    weights, masks, counts = np.asarray(weights), np.asarray(masks), np.asarray(counts)
    if (weights.ndim != 1 or not len(weights) or masks.shape != weights.shape or counts.shape != weights.shape
            or not np.issubdtype(weights.dtype, np.integer) or not np.issubdtype(masks.dtype, np.integer)
            or counts.dtype != np.uint64 or np.any(counts == 0) or np.any(weights < 0)
            or np.any(weights > 2*arity) or np.any(masks < 0) or np.any(masks > 0xffffffff)
            or type(groups) is not int or not 1 <= groups <= 18):
        raise ValueError('Bounded exact joint pair arrays and groups1..18 required')
    result = np.full((2*arity+1, 2, 2), -np.inf)
    for start in range(0, len(weights), 8192):
        if check_time: check_time()
        end = min(start+8192, len(weights))
        selected = masks[start:end]
        matrices = frozen.branch_log_matrices((selected & 65535).astype(np.uint16),
                                              (selected >> 16).astype(np.uint16), table)
        terms = matrices/groups+np.log2(counts[start:end].astype(float))[:, None, None]
        for weight in range(2*arity+1):
            take = weights[start:end] == weight
            if np.any(take):
                result[weight] = np.logaddexp2(result[weight], np.logaddexp2.reduce(terms[take], axis=0))
    return result


def prepare_source(source, catalog, table, *, check_time=None):
    source = tuple(source)
    columns, terms = [], []
    for j in range(0, 16, 2):
        if check_time: check_time()
        weights, masks, counts = catalog_api.pair_arrays(source, j, catalog)
        columns.append(pair_coefficients(weights, masks, counts, table, len(source), check_time=check_time))
        terms.append(len(counts))
    return np.stack(columns), dict(original_source=source, chronological_pairs=[[j, j+1] for j in range(0, 16, 2)],
        merged_joint_terms_per_pair=terms, exact_raw_pair_mass=256**(2*len(source)),
        distinct_universal_forms=len(catalog_api.universal_forms(len(source))),
        merged_universal_categories=len(catalog.counts), positive_mass_kernels=catalog.positive_mass_kernels)


def load_replay():
    raw = (DIRECTORY/REPLAY).read_bytes()
    if sha256(raw).hexdigest() != REPLAY_SHA: raise ValueError('Frozen replay receipt changed')
    saved = json.loads(raw)
    if (saved.get('complete') is not True or saved.get('requested_arities') != [1, 2]
            or saved.get('groups') != 18 or saved.get('source_class') != 'Fmax<=1'
            or saved.get('mark') != '9927/10000' or saved.get('normal') != 1
            or saved.get('schedule') != frozen.SCHEDULE or saved.get('source_pins') != frozen.dependencies()):
        raise ValueError('Frozen replay is not the expected complete a1/a2 control')
    rows = {tuple(row['source']): row for row in saved['fixed_A_rows']}
    expected = set(combinations(range(16), 1)) | set(combinations(range(16), 2))
    if set(rows) != expected or len(saved['fixed_A_rows']) != len(expected):
        raise ValueError('Frozen replay lacks unique full16+120 source coverage')
    return saved, rows


def run(*, mode='selected-a3', sources=SOURCES, max_seconds=120, catalog_path=None, catalog_sha256=None, progress=None):
    if mode not in ('selected-a3', 'replay-a12') or not isfinite(max_seconds) or not 0 < max_seconds <= 120:
        raise ValueError('Use a prescribed mode and wall budget at most120 seconds')
    if (catalog_path is None) != (catalog_sha256 is None) or mode == 'replay-a12' and catalog_path is not None:
        raise ValueError('A catalog path/hash pair is allowed only for selected a3')
    if mode == 'selected-a3':
        sources = tuple(map(tuple, sources))
        if (not 1 <= len(sources) <= 16 or len(set(sources)) != len(sources)
                or any(len(source) != 3 or len(set(source)) != 3 or tuple(sorted(source)) != source
                       or any(type(i) is not int or not 0 <= i < 16 for i in source) for source in sources)):
            raise ValueError('Select at most16 distinct increasing original a3 sets; no full560 sweep')
    else:
        sources = tuple(combinations(range(16), 1))+tuple(combinations(range(16), 2))
    pins, started = dependencies(), perf_counter()
    def check_time():
        if perf_counter()-started >= max_seconds: raise TimeoutError('Bounded exact SHEAR control exhausted its budget')
    bank, measure, prefix = frozen.load_measure()
    if bank['groups'] != 18 or bank['source_class'] != 'Fmax<=1' or bank['mark'] != '9927/10000':
        raise ValueError('Wrong frozen regional scope')
    old, expected = load_replay() if mode == 'replay-a12' else (None, {})
    if old is not None and (old['common_measure'] != measure or old['common_prefix_log2'] != prefix):
        raise ValueError('Replay needs the identical common matrix measure')
    rows, catalogs, diagnostics, error = [], {}, {}, None
    try:
        table = frozen.mask_log_matrices(measure, 18)
        for source in sources:
            check_time()
            a = len(source)
            if a not in catalogs:
                catalogs[a] = (catalog_api.load_catalog(catalog_path, catalog_sha256) if catalog_path is not None
                               else catalog_api.build_catalog(a, 256, check_time=check_time))
                if catalogs[a].arity != a or catalogs[a].q != 256:
                    raise ValueError('Wrong catalog field or arity')
                diagnostics[a] = dict(arity=a, universal_forms=len(catalog_api.universal_forms(a)),
                    positive_mass_kernels=catalogs[a].positive_mass_kernels, merged_categories=len(catalogs[a].counts))
                if progress: progress(dict(stage='catalog-ready', **diagnostics[a]))
            columns, data = prepare_source(source, catalogs[a], table, check_time=check_time)
            tail = float(tail_backend.positive_deficit_tail(columns, 18, 9*a*18, check_time=check_time))
            root = tail/18-8*a*log2(255)
            replay_error = abs(root-expected[source]['normalized_group_root_log2']) if old is not None else None
            if replay_error is not None and replay_error > 1e-7:
                raise ArithmeticError('Merged exact categories did not reproduce the frozen complete-A root')
            rows.append(dict(original_source=source, arity=a, enumeration=data, word_tail_log2=tail,
                normalized_group_root_log2=root, fixed_A_repeated18_bound_log2=prefix+18*root,
                aggregate_original_input_floor=9*a*18, RS_factor_per_group=f'255^(-{8*a})',
                frozen_root_absolute_error=replay_error))
            if progress: progress(dict(stage='completed-fixed-A', original_source=source,
                bound_log2=rows[-1]['fixed_A_repeated18_bound_log2'], elapsed_seconds=perf_counter()-started))
    except TimeoutError as exception:
        error = str(exception)
    if dependencies() != pins: raise ValueError('Sources or common bank changed during the control')
    if old is not None: load_replay()
    if catalog_path is not None and sha256(Path(catalog_path).read_bytes()).hexdigest() != catalog_sha256:
        raise ValueError('Catalog changed during the control')
    complete = len(rows) == len(sources)
    bound = prefix+18*float(np.logaddexp2.reduce([row['normalized_group_root_log2'] for row in rows])) if complete else None
    return tail_backend._jsonable(dict(construction='scheduled-lambda-bc-staggered-column-shear-final-D2',
        mode=mode, complete=complete, error=error, requested_sources=len(sources), completed_sources=len(rows),
        schedule=frozen.SCHEDULE, normal=1, groups=18, source_class=bank['source_class'],
        source_count=bank['source_count'], mark=bank['mark'], common_measure=measure, common_prefix_log2=prefix,
        source_and_projective_prefix_once=True, original_A_fixed_across8pairs=True,
        exact_joint_32output_supports=True, input_weights_original6variables=mode == 'selected-a3',
        independent_column_support_assumption=False, relation_kernel_shared_across_pairs=False,
        raw_counts_outside_norm_roots=True, full_independent_D1_comparison=True,
        permitted_gauges_per_pair=1, restored_D1_products_per_superblock=8, final_D2_after_shear=True,
        rows=rows, catalogs=list(diagnostics.values()), selected_rowset_sum_bound_log2=bound,
        complete_all560_original_a3=False, all_labelled_a1_a2_replayed=mode == 'replay-a12' and complete,
        partial_rows_are_not_complete_class_bound=not complete,
        replay_receipt=REPLAY if old is not None else None, replay_sha256=REPLAY_SHA if old is not None else None,
        catalog_path=str(catalog_path) if catalog_path is not None else None, catalog_sha256=catalog_sha256,
        common_measure_receipt=frozen.INPUT, common_measure_sha256=frozen.INPUT_SHA,
        reviewed_source_sha256=pins, source_pins_unchanged=True, regional_provider_runs=0,
        optimizations=0, directed=False, whole_code_certificate=False,
        elapsed_seconds=perf_counter()-started, max_seconds=max_seconds))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=('selected-a3', 'replay-a12'), default='selected-a3')
    parser.add_argument('--source', nargs=3, type=int, action='append')
    parser.add_argument('--catalog', type=Path)
    parser.add_argument('--catalog-sha256')
    parser.add_argument('--max-seconds', type=float, default=120)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists(): parser.error('Refuse to overwrite evidence')
    result = run(mode=args.mode, sources=args.source or SOURCES, catalog_path=args.catalog,
        catalog_sha256=args.catalog_sha256, max_seconds=args.max_seconds,
        progress=lambda row: print(json.dumps(row, allow_nan=False), flush=True))
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
