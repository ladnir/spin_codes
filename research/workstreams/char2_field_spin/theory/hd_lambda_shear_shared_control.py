"""Bounded selected a4 SHEAR gate using the shared-ratio majorant.

Fix original A throughout the eight ordered adjacent-column pairs. The
helper returns a branch-root polynomial: each product of raw column counts
multiplies (254/255)||M_OR||_18+(1/255)||M_XOR||_18. Those ratio weights
and raw counts are already outside the roots. Do NOT take another root
of a returned coefficient. Raise each entry to18, compose chronologically,
and retain original input degree>=9*|A|*18. Normalize each completed-A
word root by255^(-8|A|), then use the common prefix once.

One frozen positive regional measure is prepared once. Its scalar atoms
are in[0,1], as required by the pointwise cancellation-class comparison.
The final independent D2 remains after the shear; only one pair gauge is
removed from ideal independent D1. No independent-column support law is used.

Defaults are A012, A0123, and A0148. The A012 comparison reports the old
complete bound independently. The new norm majorant need not exceed the old
branchwise Minkowski bound; no such ordering or extra minimum is assumed.
At most eight selected sets and120 seconds are allowed. There is no all1820
enumeration, new provider, optimization, native test, or whole-code claim.
"""
from hashlib import sha256
import json
from math import isfinite, log2
from pathlib import Path
from time import perf_counter

import numpy as np

import hd_lambda_shear_shared_ratio as ratio
import hd_lambda_column_shear_control as frozen
import hd_lambda_double_middle_control as tail_backend


DIRECTORY = Path(__file__).resolve().parent
SOURCES = ((0, 1, 2), (0, 1, 2, 3), (0, 1, 4, 8))
REFERENCE = 'hd_lambda_shear_exact3_selected_float.json'
REFERENCE_SHA = '7df7cbfdc687a8e561f5d2c6ea1d3303bdae6e837920b5ead1e78b698c9abac6'
PINS = {
    'hd_lambda_shear_shared_ratio.py': '5b60c22fa6e120c9e4faee276de0f029dc52121c21dc53c2d7c12eafdac9b363',
    'hd_lambda_column_shear_control.py': 'a1ede595d659ddd4c20d1b5ea33a98158e7c9484e34a958d2d682f985214de96',
    'hd_lambda_double_middle_control.py': '78655f83f08ae01118beb6d9f6258daaddd7a5035bc5837694b6243d762e9820',
}


def dependencies():
    for name, digest in PINS.items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Reviewed helper changed: '+name)
    return {**ratio.dependencies(), **frozen.dependencies(), **tail_backend.dependencies(), **PINS,
            Path(__file__).name: sha256(Path(__file__).read_bytes()).hexdigest()}


def validate_sources(sources):
    sources = tuple(map(tuple, sources))
    if (not 1 <= len(sources) <= 8 or len(set(sources)) != len(sources)
            or any(len(source) not in (3, 4) or tuple(sorted(source)) != source
                   or len(set(source)) != len(source)
                   or any(type(i) is not int or not 0 <= i < 16 for i in source) for source in sources)):
        raise ValueError('One to eight distinct increasing original a3/a4 sets required')
    return sources


def prepare_source(source, prepared, *, check_time=None):
    source = validate_sources((source,))[0]
    if prepared.groups != 18 or prepared.rows != 16:
        raise ValueError('One r18/sixteen-row common measure is required')
    pairs = []
    for column in range(0, 16, 2):
        if check_time: check_time()
        polynomial = ratio.pair_polynomial(source, column, prepared, q=256,
            schedule=frozen.SCHEDULE, chunk_size=8192, check_time=check_time)
        if np.shape(polynomial) != (2*len(source)+1, 2, 2):
            raise ValueError('Shared-ratio helper returned the wrong original-degree polynomial')
        pairs.append(polynomial)
    return np.stack(pairs)


def load_reference(measure, prefix):
    raw = (DIRECTORY/REFERENCE).read_bytes()
    if sha256(raw).hexdigest() != REFERENCE_SHA: raise ValueError('Pinned a3 reference changed')
    saved = json.loads(raw)
    if (saved.get('construction') != 'scheduled-lambda-bc-staggered-column-shear-final-D2'
            or saved.get('complete') is not True or saved.get('groups') != 18
            or saved.get('source_class') != 'Fmax<=1' or saved.get('mark') != '9927/10000'
            or saved.get('schedule') != frozen.SCHEDULE or saved.get('normal') != 1
            or saved.get('common_measure') != measure or saved.get('common_prefix_log2') != prefix):
        raise ValueError('Reference needs the same construction, measure, and prefix')
    rows = {tuple(row['original_source']): row for row in saved['rows']}
    if len(rows) != len(saved['rows']): raise ValueError('Duplicate reference row set')
    return rows


def run(*, sources=SOURCES, max_seconds=120, progress=None):
    sources = validate_sources(sources)
    if type(max_seconds) not in (int, float) or not isfinite(max_seconds) or not 0 < max_seconds <= 120:
        raise ValueError('Use a positive wall budget of at most120 seconds')
    pins, started = dependencies(), perf_counter()
    def check_time():
        if perf_counter()-started >= max_seconds: raise TimeoutError('Shared-ratio selected gate exhausted its budget')
    bank, measure, prefix = frozen.load_measure()
    if bank['groups'] != 18 or bank['source_class'] != 'Fmax<=1' or bank['mark'] != '9927/10000':
        raise ValueError('Wrong common regional scope')
    references = load_reference(measure, prefix) if any(len(source) == 3 for source in sources) else {}
    rows, error, prepared_count = [], None, 0
    try:
        check_time()
        prepared = ratio.prepare_measure(measure, groups=18)
        prepared_count = 1
        for source in sources:
            check_time()
            a = len(source)
            polynomials = prepare_source(source, prepared, check_time=check_time)
            tail = float(tail_backend.positive_deficit_tail(polynomials, 18, 9*a*18, check_time=check_time))
            root = tail/18-8*a*log2(255)
            bound = prefix+18*root
            reference = references.get(source)
            rows.append(dict(original_source=source, arity=a, status='complete',
                word_tail_log2=tail, normalized_group_root_log2=root, fixed_A_repeated18_bound_log2=bound,
                aggregate_original_input_floor=9*a*18, retained_deficit=7*a*18,
                RS_factor_per_group=f'255^(-{8*a})', RS_factor_whole_word=f'255^(-{144*a})',
                prior_branchwise_bound_log2=reference['fixed_A_repeated18_bound_log2'] if reference else None,
                new_minus_prior_bound_log2=bound-reference['fixed_A_repeated18_bound_log2'] if reference else None,
                ordering_between_two_valid_majorants_asserted=False))
            if progress:
                progress(tail_backend._jsonable(dict(stage='completed-fixed-A', original_source=source,
                    bound_log2=bound, elapsed_seconds=perf_counter()-started)))
    except TimeoutError as exception:
        error = str(exception)
    if dependencies() != pins: raise ValueError('Sources or bank changed during the control')
    if references: load_reference(measure, prefix)
    complete = len(rows) == len(sources)
    selected = prefix+18*float(np.logaddexp2.reduce([row['normalized_group_root_log2'] for row in rows])) if complete else None
    return tail_backend._jsonable(dict(construction='scheduled-lambda-bc-staggered-column-shear-final-D2',
        arithmetic='positive-FLOAT-log2-not-directed', complete=complete, timeout=error,
        requested_original_sources=sources, completed_sources=len(rows), groups=18,
        source_class=bank['source_class'], source_count=bank['source_count'], mark=bank['mark'],
        common_measure=measure, common_prefix_log2=prefix, schedule=frozen.SCHEDULE, normal=1,
        original_A_fixed_across8pairs=True, rows=rows, selected_rowset_sum_bound_log2=selected,
        full_a4_class_covered=False, representative_extrapolation=False,
        partial_rows_are_not_complete_class_bound=not complete,
        exact_single_column_counts=True, exact_joint_support_distribution_claimed=False,
        shared_ratio_weights=['254/255', '1/255'], shared_ratio_weights_outside_norm_roots=True,
        raw_count_products_outside_norm_roots=True, helper_result_rooted_again=False,
        entrywise_18_power_before_chronological_matrix_products=True,
        source_and_projective_prefix_once=True, independent_column_support_assumption=False,
        full_independent_D1_comparison=True, permitted_gauges_per_pair=1,
        restored_D1_products_per_superblock=8, final_D2_after_shear=True,
        reference_receipt=REFERENCE if references else None, reference_sha256=REFERENCE_SHA if references else None,
        common_measure_receipt=frozen.INPUT, common_measure_sha256=frozen.INPUT_SHA,
        reviewed_source_sha256=pins, source_pins_unchanged=True, measure_preparations=prepared_count,
        helper_metadata=ratio.metadata(), regional_provider_runs=0, input_marker_searches=0,
        directed=False, whole_code_certificate=False, native_operations=0,
        elapsed_seconds=perf_counter()-started, max_seconds=max_seconds))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', nargs='+', type=int, action='append')
    parser.add_argument('--include-a3', action='store_true')
    parser.add_argument('--max-seconds', type=float, default=120)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists(): parser.error('Refuse to overwrite evidence')
    sources = tuple(map(tuple, args.source)) if args.source else SOURCES
    if args.include_a3 and (0, 1, 2) not in sources: sources += ((0, 1, 2),)
    result = run(sources=sources, max_seconds=args.max_seconds,
                 progress=lambda row: print(json.dumps(row, allow_nan=False), flush=True))
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
