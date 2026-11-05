"""All560 labelled a3 sets through the frozen exact SHEAR evaluator.

One pinned six-input/80-form catalog and one fixed regional moment table
serve every source set. Each original A stays fixed through its eight ordered
pairs. The unchanged evaluator retains total input weight>=486 after the
entrywise18th powers and chronological matrix product. Each complete-A root
is tail/18-24log2(255). Sum those roots once, raise to18, and apply the common
source/output/projective prefix once.

The full-class result exists only after all560 unique original three-row sets
complete with unchanged pins. Timeout retains completed rows but no full-class
bound. There is no orbit quotient, extra factor560, fresh provider, fitting,
or transfer of an old outer certificate. This is FLOAT evidence for r18 and
Fmax<=1, not whole-code certification. Root owns the production run.
"""
from hashlib import sha256
from itertools import combinations
import json
from math import isfinite, log2
from pathlib import Path
from time import perf_counter

import numpy as np

import hd_lambda_shear_exact3_control as core


DIRECTORY = Path(__file__).resolve().parent
CATALOG = 'hd_lambda_shear_categories_a3_q256.npz'
CATALOG_SHA = '4a515aca09ca731da1c171190cc9f1c8c7bed000c220d2d1e55cba8398c64116'
CORE_SHA = '6f0261a5e5eb1c24cabd5fdbe3ead485acf109c21668043471eec64c14921eba'
SOURCES = tuple(combinations(range(16), 3))
GROUPS, ARITY, FLOOR = 18, 3, 486


def dependencies():
    if sha256((DIRECTORY/'hd_lambda_shear_exact3_control.py').read_bytes()).hexdigest() != CORE_SHA:
        raise ValueError('Frozen selected evaluator changed')
    return {**core.dependencies(), Path(__file__).name: sha256(Path(__file__).read_bytes()).hexdigest()}


def coverage(rows, prefix):
    """Validate completed-A records; no extrapolation from a proper subset."""
    if not isfinite(prefix): raise ValueError('A finite common prefix is required')
    expected, seen, roots = set(SOURCES), set(), []
    for row in rows:
        source = tuple(row['original_source'])
        if (any(type(i) is not int for i in source) or source not in expected or source in seen or row.get('arity') != ARITY
                or row.get('status') != 'complete' or row.get('aggregate_original_input_floor') != FLOOR):
            raise ValueError('Duplicate, foreign, or incomplete original-A record')
        root, tail = row['normalized_group_root_log2'], row['word_tail_log2']
        if (not (isfinite(root) or root == -np.inf) or not (isfinite(tail) or tail == -np.inf)
                or (root == -np.inf) != (tail == -np.inf)):
            raise ValueError('Invalid complete-A log bound')
        if isfinite(root) and abs(root-(tail/GROUPS-8*ARITY*log2(255))) > 1e-9:
            raise ValueError('Actual arity-three RS normalization changed')
        seen.add(source)
        roots.append(root)
    complete = seen == expected
    bound = prefix+GROUPS*float(np.logaddexp2.reduce(roots)) if complete else None
    worst = rows[int(np.argmax(roots))]['original_source'] if roots else None
    return dict(complete=complete, unique_original_rowsets=len(seen), expected_original_rowsets=560,
        missing_original_rowsets=560-len(seen), all_a3_bound_log2=bound,
        worst_completed_original_source=worst,
        largest_completed_fixed_A_bound_log2=(prefix+GROUPS*max(roots) if roots else None))


def run(*, max_seconds=180, catalog_path=None, progress=None):
    if type(max_seconds) not in (int, float) or not isfinite(max_seconds) or not 0 < max_seconds <= 300:
        raise ValueError('Use a positive wall budget of at most300 seconds')
    path = DIRECTORY/CATALOG if catalog_path is None else Path(catalog_path)
    pins, started = dependencies(), perf_counter()
    if sha256(path.read_bytes()).hexdigest() != CATALOG_SHA:
        raise ValueError('Complete a3 catalog differs from the explicit frozen pin')
    def check_time():
        if perf_counter()-started >= max_seconds: raise TimeoutError('All-a3 gate exhausted its wall budget')
    bank, measure, prefix = core.frozen.load_measure()
    if bank['groups'] != GROUPS or bank['source_class'] != 'Fmax<=1' or bank['mark'] != '9927/10000':
        raise ValueError('The common regional scope changed')
    rows, error, catalog_metadata, table_builds = [], None, None, 0
    try:
        check_time()
        catalog = core.catalog_api.load_catalog(path, CATALOG_SHA)
        if (catalog.arity != ARITY or catalog.q != 256 or len(catalog.counts) != 66795
                or catalog.positive_mass_kernels != 564357):
            raise ValueError('Unexpected complete six-input category inventory')
        catalog_metadata = dict(arity=3, field_order=256, merged_categories=66795,
            positive_mass_kernels=564357, universal_forms=80)
        table = core.frozen.mask_log_matrices(measure, GROUPS)
        table_builds = 1
        for source in SOURCES:
            check_time()
            columns, diagnostics = core.prepare_source(source, catalog, table, check_time=check_time)
            tail = float(core.tail_backend.positive_deficit_tail(columns, GROUPS, FLOOR, check_time=check_time))
            root = tail/GROUPS-8*ARITY*log2(255)
            rows.append(dict(original_source=source, arity=ARITY, status='complete',
                word_tail_log2=tail, normalized_group_root_log2=root,
                fixed_A_repeated18_bound_log2=prefix+GROUPS*root,
                aggregate_original_input_floor=FLOOR,
                merged_joint_terms_per_pair=diagnostics['merged_joint_terms_per_pair']))
            if progress and (len(rows) % 32 == 0 or len(rows) == len(SOURCES)):
                progress(core.tail_backend._jsonable(dict(stage='completed-original-rowsets',
                    completed=len(rows), requested=560, last_original_source=source,
                    largest_fixed_A_bound_log2=max(row['fixed_A_repeated18_bound_log2'] for row in rows),
                    elapsed_seconds=perf_counter()-started)))
    except TimeoutError as exception:
        error = str(exception)
    if dependencies() != pins:
        raise ValueError('Reviewed sources or regional bank changed during all-a3 gate')
    if sha256(path.read_bytes()).hexdigest() != CATALOG_SHA:
        raise ValueError('Pinned catalog changed during all-a3 gate')
    summary = coverage(rows, prefix)
    return core.tail_backend._jsonable(dict(construction='scheduled-lambda-bc-staggered-column-shear-final-D2',
        arithmetic='positive-FLOAT-log2-not-directed', **summary, timeout=error,
        groups=GROUPS, arity=ARITY, source_class=bank['source_class'], source_count=bank['source_count'],
        mark=bank['mark'], schedule=core.frozen.SCHEDULE, normal=1,
        common_measure=measure, common_prefix_log2=prefix, rows=rows,
        complete_all560_original_a3=summary['complete'], partial_rows_are_not_full_class_bound=not summary['complete'],
        original_A_fixed_across8pairs=True, labelled_rowset_sum_once_before18th_power=True,
        extra_rowset_factor_applied=False, source_and_projective_prefix_once=True,
        aggregate_original_input_floor=FLOOR, retained_deficit=378,
        RS_factor_per_group='255^(-24)', RS_factor_whole_word='255^(-432)',
        exact_joint_32output_supports=True, independent_column_support_assumption=False,
        full_independent_D1_comparison=True, permitted_gauges_per_pair=1,
        restored_D1_products_per_superblock=8, final_D2_after_shear=True,
        catalog_path=str(path), catalog_sha256=CATALOG_SHA, catalog_metadata=catalog_metadata,
        common_measure_receipt=core.frozen.INPUT, common_measure_sha256=core.frozen.INPUT_SHA,
        reviewed_source_sha256=pins, source_pins_unchanged=True,
        regional_provider_runs=0, catalog_builds=0, moment_table_builds=table_builds,
        optimizations=0, directed=False, whole_code_certificate=False,
        elapsed_seconds=perf_counter()-started, max_seconds=max_seconds))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--max-seconds', type=float, default=180)
    parser.add_argument('--catalog', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists(): parser.error('Refuse to overwrite proof evidence')
    result = run(max_seconds=args.max_seconds, catalog_path=args.catalog,
                 progress=lambda row: print(json.dumps(row, allow_nan=False), flush=True))
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps({key: result[key] for key in ('complete', 'unique_original_rowsets', 'all_a3_bound_log2', 'timeout', 'elapsed_seconds')}, allow_nan=False), flush=True)
