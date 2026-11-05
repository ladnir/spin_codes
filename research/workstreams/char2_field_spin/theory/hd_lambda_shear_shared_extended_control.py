"""Bounded selected-a5/a6 FLOAT control for the original Hlambda shear.

One exact single-column catalog per arity supplies the shared-ratio
majorant.  Its coefficients are already L18 norm bounds.  Raise each
entry polynomial to18, multiply eight factors chronologically, and
retain aggregate ORIGINAL input activity>=9*a*18.  No marker is fitted.

The new tail wrapper permits pair degree12, unlike the frozen degree11
dispatcher.  It reverses degrees and discards only excess deficits; all
convolutions are positive.  No discarded term can return to the prefix.
Only arities5/6 and at most eight prescribed original sets are evaluated.

For each completed original A, the normalized root is tail_log/18 minus
8*a*log2(255).  Sum these roots before the group power and apply the
cached source/Chernoff/projective prefix once.  No binomial row-set factor
is added.  Selected results do not cover missing labelled sets or other
source profiles.  The bank is the unchanged r18/Fmax<=1, z=9927/10000
common matrix measure.  No new regional provider or native work occurs.
"""
from hashlib import sha256
import json
from math import isfinite, log2
from pathlib import Path
from time import perf_counter

import numpy as np

import hd_lambda_shear_shared_ratio_extended as extended
import hd_lambda_shear_shared_ratio as ratio
import hd_lambda_column_shear_control as frozen
import hd_lambda_double_middle_control as generic
from hd_pair_state_tail import _entrywise_convolution


DIRECTORY = Path(__file__).resolve().parent
SOURCES = (tuple(range(5)), tuple(range(6)))
PINS = {
    'hd_lambda_shear_shared_ratio_extended.py': 'b51f0869028327137f58c95d260a86653cc2b4acfb50275b6ab346cd8db0c01f',
    'hd_lambda_column_shear_control.py': 'a1ede595d659ddd4c20d1b5ea33a98158e7c9484e34a958d2d682f985214de96',
    'hd_lambda_double_middle_control.py': '78655f83f08ae01118beb6d9f6258daaddd7a5035bc5837694b6243d762e9820',
    'hd_pair_state_tail.py': '2af63d37d228a76a26853d7550815de5f6c0bc5af05908a5b009e14aa9175858',
}


def dependencies():
    result = {**extended.dependencies(), **frozen.dependencies(), **generic.dependencies()}
    for name, digest in PINS.items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Frozen source changed: '+name)
        result[name] = digest
    result[Path(__file__).name] = sha256(Path(__file__).read_bytes()).hexdigest()
    return result


def positive_pair_deficit_tail(pairs, groups, minimum_input, *, check_time=None):
    """Positive degree<=12 tail for at most8 pairs, total degree<=1728.

    This is a new explicit dispatcher.  It neither edits nor calls the
    old degree<=11 dispatcher.  Its two reused convolution primitives
    have dimension-independent semantics and compatible workspace guards.
    """
    pairs = np.asarray(pairs, dtype=float)
    if (pairs.ndim < 4 or not 1 <= len(pairs) <= 8 or pairs.shape[-2:] != (2, 2)
            or not 2 <= pairs.shape[-3] <= 13 or np.any(np.isnan(pairs))
            or np.any(np.isposinf(pairs)) or type(groups) is not int or not 1 <= groups <= 18):
        raise ValueError('Up to8 ordered pair polynomials of degree1..12 and groups1..18 required')
    degree = pairs.shape[-3]-1
    total = len(pairs)*degree*groups
    if type(minimum_input) is not int or not 0 <= minimum_input <= total:
        raise ValueError('Input threshold must lie in the bounded full degree range')
    limit = total-minimum_input
    word = np.full(pairs.shape[1:-3]+(1, 2, 2), -np.inf)
    word[..., 0, 0, 0] = word[..., 0, 1, 1] = 0.
    for pair in pairs:
        if check_time is not None:
            check_time()
        base = pair[..., ::-1, :, :]
        factor = np.zeros(pair.shape[:-3]+(1, 2, 2))
        for _ in range(groups):
            if check_time is not None:
                check_time()
            factor = _entrywise_convolution(factor, base)[..., :limit+1, :, :]
        word = generic.matrix_product_truncated(word, factor, limit, check_time=check_time)
    return np.logaddexp2.reduce(np.logaddexp2.reduce(word[..., 0, :], axis=-1), axis=-1)


def validate_sources(sources):
    sources = tuple(map(tuple, sources))
    if (not 1 <= len(sources) <= 8 or len(set(sources)) != len(sources)
            or any(len(source) not in (5, 6) or tuple(sorted(source)) != source
                   or len(set(source)) != len(source)
                   or any(type(i) is not int or not 0 <= i < 16 for i in source) for source in sources)):
        raise ValueError('One to eight distinct increasing original a5/a6 sets required')
    return sources


def prepare_source(source, catalog, prepared, *, max_products=extended.MAX_PAIR_PRODUCTS,
                   check_time=None, progress=None):
    source = validate_sources((source,))[0]
    if catalog.arity != len(source) or catalog.q != 256 or prepared.groups != 18 or prepared.rows != 16:
        raise ValueError('Matching exact catalog and r18/sixteen-row common kernel required')
    pairs, diagnostics = [], []
    for column in range(0, 16, 2):
        if check_time is not None:
            check_time()
        left = extended.column_histogram(source, column, catalog)
        right = extended.column_histogram(source, column+1, catalog)
        polynomial = extended.pair_coefficients(left, right, prepared, q=256,
            chunk_size=8192, max_products=max_products, check_time=check_time)
        if polynomial.shape != (2*len(source)+1, 2, 2):
            raise ValueError('Wrong original-input degree from shared-ratio contraction')
        pairs.append(polynomial)
        info = dict(column=column, mate=column+1, left_histogram=len(left), right_histogram=len(right),
                    support_products=len(left)*len(right), exact_raw_pair_mass=str(256**(2*len(source))))
        diagnostics.append(info)
        if progress is not None:
            progress(dict(stage='completed-pair', original_source=source, **info))
    return np.stack(pairs), diagnostics


def run(*, sources=SOURCES, max_seconds=120, max_products=extended.MAX_PAIR_PRODUCTS, progress=None):
    sources = validate_sources(sources)
    if (type(max_seconds) not in (int, float) or not isfinite(max_seconds) or not 0 < max_seconds <= 120
            or type(max_products) is not int or not 1 <= max_products <= extended.MAX_PAIR_PRODUCTS):
        raise ValueError('Time budget(0,120] and pair-product budget1..1M required')
    pins, started = dependencies(), perf_counter()
    def check_time():
        if perf_counter()-started >= max_seconds:
            raise TimeoutError('Selected a5/a6 shared-ratio budget exhausted')
    bank, measure, prefix = frozen.load_measure()
    if bank['groups'] != 18 or bank['source_class'] != 'Fmax<=1' or bank['mark'] != '9927/10000':
        raise ValueError('Wrong frozen regional reference scope')
    rows, catalog_rows, error, stopped_source = [], [], None, None
    preparations, catalogs = 0, {}
    try:
        check_time()
        prepared = ratio.prepare_measure(measure, groups=18)
        preparations = 1
        for source in sources:
            stopped_source = source
            check_time()
            a = len(source)
            if a not in catalogs:
                if progress is not None:
                    progress(dict(stage='building-column-catalog', arity=a, plan=extended.catalog_plan(a)))
                catalogs[a] = extended.build_catalog(a, q=256, check_time=check_time)
                info = dict(arity=a, field_order=256, positive_mass_kernels=catalogs[a].positive_mass_kernels,
                    merged_categories=len(catalogs[a].counts), dense_forms=1 << a, exact_column_mass=str(256**a))
                catalog_rows.append(info)
                if progress is not None:
                    progress(dict(stage='completed-column-catalog', **info))
            polynomials, diagnostics = prepare_source(source, catalogs[a], prepared,
                max_products=max_products, check_time=check_time, progress=progress)
            tail = float(positive_pair_deficit_tail(polynomials, 18, 9*a*18, check_time=check_time))
            root = tail/18-8*a*log2(255)
            bound = prefix+18*root
            rows.append(dict(original_source=source, arity=a, status='complete', pairs=diagnostics,
                word_tail_log2=tail, normalized_group_root_log2=root, fixed_A_repeated18_bound_log2=bound,
                aggregate_original_input_floor=9*a*18, retained_deficit=7*a*18,
                RS_factor_per_group=f'255^(-{8*a})'))
            if progress is not None:
                progress(generic._jsonable(dict(stage='completed-fixed-A', original_source=source,
                    bound_log2=bound, elapsed_seconds=perf_counter()-started)))
        stopped_source = None
    except (TimeoutError, extended.BudgetExceeded) as exception:
        error = str(exception)
    if dependencies() != pins:
        raise ValueError('Frozen inputs changed during selected control')
    complete = len(rows) == len(sources)
    selected = prefix+18*float(np.logaddexp2.reduce([row['normalized_group_root_log2'] for row in rows])) if complete else None
    return generic._jsonable(dict(construction='scheduled-lambda-bc-staggered-column-shear-final-D2',
        arithmetic='positive-FLOAT-log2-not-directed', complete=complete, stop_reason=error,
        stopped_original_source=stopped_source, requested_original_sources=sources, rows=rows,
        catalogs=catalog_rows, groups=18, source_class=bank['source_class'], source_count=bank['source_count'],
        mark=bank['mark'], common_measure=measure, common_prefix_log2=prefix,
        selected_rowset_sum_bound_log2=selected, all_labelled_class_covered=False,
        representative_extrapolation=False, partial_rows_are_not_complete_class_bound=not complete,
        schedule=frozen.SCHEDULE, normal=1, original_A_fixed_across8pairs=True,
        source_and_projective_prefix_once=True, helper_result_rooted_again=False,
        shared_ratio_weights=['254/255', '1/255'], shared_ratio_weights_outside_norm_roots=True,
        raw_count_products_outside_norm_roots=True, exact_single_column_counts=True,
        exact_joint_output_distribution_claimed=False, full_independent_D1_comparison=True,
        permitted_gauges_per_pair=1, restored_D1_products_per_superblock=8, final_D2_after_shear=True,
        common_measure_receipt=frozen.INPUT, common_measure_sha256=frozen.INPUT_SHA,
        source_pins=pins, source_pins_unchanged=True, measure_preparations=preparations,
        helper_metadata=extended.metadata(), regional_provider_runs=0, marker_searches=0,
        directed=False, whole_code_certificate=False, native_operations=0,
        max_seconds=max_seconds, max_pair_products=max_products, elapsed_seconds=perf_counter()-started))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', nargs='+', type=int, action='append')
    parser.add_argument('--max-seconds', type=float, default=120)
    parser.add_argument('--max-products', type=int, default=extended.MAX_PAIR_PRODUCTS)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Refuse to overwrite evidence')
    result = run(sources=args.source if args.source else SOURCES, max_seconds=args.max_seconds,
        max_products=args.max_products, progress=lambda row: print(json.dumps(row, allow_nan=False), flush=True))
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
