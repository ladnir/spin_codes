"""Fresh bounded middle-arity FLOAT gates for the64-XOR double transvection.

One pinned, outer-independent r18/Fmax<=1 regional measure is reused.
All new outer counts come from the actual double-transvection matrix.
For each fixed original A, retain its sixteen chronological P_j/Q_j maps.
Zero-outside output degrees use the proved protected-coordinate surrogate;
nonzero-outside output degrees are inflated to|A| before coset domination.
This is upper-tail domination, not coefficientwise input-weight equality.

Positive polynomial powers retain aggregate original-input activity>=9ar.
The fixed-A norm is Tail_A^(1/r)/255^(8a). The source/output/projective
prefix occurs once. Representatives are not extrapolated to missing A.
The separate pair mode covers all120 original two-row sets and sums their
complete word roots before the rth power. No old outer receipt transfers.

The all-nonquiet word component is reported BEFORE the full deficit tail.
It has surrogate degree16ar and survives every input-degree threshold used
here. It is a LOWER bound on this proof majorant, not on failure probability.
A large value diagnoses the nonquiet relaxation rather than a construction
failure. There is no input-marker search and no new regional provider.
"""
from hashlib import sha256
from itertools import combinations
import json
from math import isfinite, log2
from pathlib import Path
from time import perf_counter

import numpy as np

import hd_lambda_double_fiber as fiber
from hd_lambda_fiber_control import moment_table
from hd_lambda_marked_fiber_control import nonquiet_norm
from hd_pair_scheduled_control import load_measure
from hd_pair_state_tail import _entrywise_convolution
from hd_post_refresh_quiet_rows import log_matmul
from hd_inverse_column_control import INPUT, INPUT_SHA
from hd_lambda_fiber_outer import _permute_mask
from hd_pair_scheduled_outer import column_permutations


DIRECTORY = Path(__file__).resolve().parent
GROUPS = 18
SOURCES = (tuple(range(6)), tuple(range(7)), tuple(range(8)), (0, 1, 2, 3, 4, 5, 8, 12),
           tuple(range(9)), tuple(range(10)), tuple(range(11)),
           tuple(i for i in range(16) if i not in (0, 7, 11, 13, 14)))
PINS = {
    'hd_lambda_double_fiber.py': '30108ff6fead25311ba9f9b704d0ba5594e052433fb1a64d11217ed7d88f4905',
    'hd_lambda_fiber_control.py': '9e0c5157053ced38562263c6af76aa640e935221f7365b0a772f4db036c0517d',
    'hd_lambda_marked_fiber_control.py': '95dd74e7df2f4d66760c11d03c213cf2aebb1ab5c73c35b136f99cecdee2235f',
    'hd_pair_scheduled_control.py': 'cd0ee6b03e794785dee26502940b3f68395cd9a1fd6f27dc6f3470a490730dd4',
    'hd_pair_scalar_control.py': '088f33f37dc70b9da7d4649d7f56e8657fd84d4945c909bf49235cc20d5c3f61',
    'hd_pair_scheduled_outer.py': '64dc50bf18171d371ae16f4484087e9d42d89cc7c9b55a7576be55b1251e4741',
    'hd_pair_identity_kernel.py': '93f3a9053a72bed39db6b00a32d2f2b7b060a35ce19ce814dd856d096a796961',
    'hd_pair_state_kernel.py': 'b28142e19e51a906ae742128325d1e7557dfb425624c6da45778a2a75753e195',
    'hd_post_refresh_quiet_rows.py': '75e38819d225e9771ade1ed90509b7d4582a2637b24b6522f5b62b0200eb49e7',
    'hd_pair_state_tail.py': '2af63d37d228a76a26853d7550815de5f6c0bc5af05908a5b009e14aa9175858',
    'hd_inverse_column_control.py': '40668f4bf673ef0859be429a82119769df169408d74aaaa52dfbb68787b9716b',
}


def dependencies():
    for name, digest in PINS.items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Reviewed dependency changed: '+name)
    return {**fiber.checked_dependencies(), **PINS,
            Path(__file__).name: sha256(Path(__file__).read_bytes()).hexdigest()}


def _logs(values):
    result = np.asarray(values, dtype=float)
    if np.any(np.isnan(result)) or np.any(np.isposinf(result)):
        raise ValueError('Finite log2 values or exact zero(-inf) required')
    return result


def matrix_product_truncated(left, right, limit, *, check_time=None):
    """Ordered positive convolution, truncating before any discarded allocation."""
    left, right = _logs(left), _logs(right)
    if (type(limit) is not int or not 0 <= limit <= 3168 or left.ndim < 3 or right.ndim < 3
            or left.shape[-2:] != (2, 2) or right.shape[-2:] != (2, 2)
            or not 1 <= left.shape[-3] <= 3169 or not 1 <= right.shape[-3] <= 3169):
        raise ValueError('Bounded2x2 positive matrix polynomials required')
    a, b = left.shape[-3], right.shape[-3]
    length = min(limit+1, a+b-1)
    shape = np.broadcast_shapes(left.shape[:-3], right.shape[:-3])
    result = np.full(shape+(length, 2, 2), -np.inf)
    if a <= b:
        for degree in range(min(a, length)):
            if check_time: check_time()
            take = min(b, length-degree)
            fixed, varying = left[..., degree:degree+1, :, :], right[..., :take, :, :]
            term = np.logaddexp2(fixed[..., :, 0, None]+varying[..., None, 0, :],
                                fixed[..., :, 1, None]+varying[..., None, 1, :])
            dest = result[..., degree:degree+take, :, :]
            np.logaddexp2(dest, term, out=dest)
    else:
        for degree in range(min(b, length)):
            if check_time: check_time()
            take = min(a, length-degree)
            fixed, varying = right[..., degree:degree+1, :, :], left[..., :take, :, :]
            term = np.logaddexp2(varying[..., :, 0, None]+fixed[..., None, 0, :],
                                varying[..., :, 1, None]+fixed[..., None, 1, :])
            dest = result[..., degree:degree+take, :, :]
            np.logaddexp2(dest, term, out=dest)
    return result


def positive_deficit_tail(columns, groups, minimum_input, *, check_time=None):
    """One to16 columns, arity<=11; every discarded deficit stays discarded."""
    columns = _logs(columns)
    if (columns.ndim < 4 or not 1 <= len(columns) <= 16 or columns.shape[-2:] != (2, 2)
            or not 2 <= columns.shape[-3] <= 12 or type(groups) is not int or not 1 <= groups <= 18):
        raise ValueError('Ordered column,batch?,degree,2,2 array and bounded groups required')
    arity = columns.shape[-3]-1
    total = len(columns)*arity*groups
    if type(minimum_input) is not int or not 0 <= minimum_input <= total:
        raise ValueError('Original-input floor must fit the full degree range')
    limit = total-minimum_input
    word = np.full(columns.shape[1:-3]+(1, 2, 2), -np.inf)
    word[..., 0, 0, 0] = word[..., 0, 1, 1] = 0.
    for column in columns:
        if check_time: check_time()
        base = column[..., ::-1, :, :]
        factor = np.zeros(column.shape[:-3]+(1, 2, 2))
        for _ in range(groups):
            if check_time: check_time()
            factor = _entrywise_convolution(factor, base)[..., :limit+1, :, :]
        word = matrix_product_truncated(word, factor, limit, check_time=check_time)
    return np.logaddexp2.reduce(np.logaddexp2.reduce(word[..., 0, :], axis=-1), axis=-1)


def prepare_source(source, table_r, table_2r, *, groups=18, check_time=None):
    """New exact fibre counts; same reviewed unmarked nonquiet-norm inequality."""
    source = fiber._source(source)
    a = len(source)
    if not 5 <= a <= 11:
        raise ValueError('This fibre preparation allows arities5..11')
    columns, nonquiet, diagnostics = [], [], []
    for j in range(16):
        if check_time: check_time()
        data = fiber.column_fiber_counts(source, j)
        bounds = nonquiet_norm(data, table_r, table_2r, groups, check_time=check_time)
        coefficients = np.full((a+1, 2, 2), -np.inf)
        for mask, count in data['zero_fiber_counts'].items():
            degree = fiber.protected_degree(mask, data)
            coefficients[degree] = np.logaddexp2(coefficients[degree], log2(count)+table_r[mask]/groups)
        coefficients[a] = np.logaddexp2(coefficients[a], bounds['selected'])
        columns.append(coefficients)
        nonquiet.append(bounds['selected'])
        diagnostics.append(dict(column=j, mapped_source=data['mapped_source'],
            cross_rank=data['cross_rank'], syndrome_rank=data['syndrome_rank'],
            protected_coordinates=data['physical_protected_mask'].bit_count(),
            physical_protected_mask=data['physical_protected_mask'],
            zero_fiber_pointwise_quiet=data['all_zero_fiber_coordinates_protected'],
            outside_supports=len(data['outside_counts']), zero_fiber_supports=len(data['zero_fiber_counts']),
            direct_union_terms=bounds['direct_terms'], factor_norm_log2=bounds['factor_2r'].tolist(),
            direct_norm_log2=bounds['direct_union'].tolist()))
    return np.stack(columns), np.stack(nonquiet), diagnostics


def all_nonquiet_diagnostic(nonquiet, arity, prefix, *, groups=18):
    nonquiet = _logs(nonquiet)
    if nonquiet.shape != (16, 2, 2):
        raise ValueError('Sixteen chronological nonquiet2x2 norm matrices required')
    word = np.array(((0., -np.inf), (-np.inf, 0.)))
    for column in nonquiet:
        word = log_matmul(word, groups*column)
    value = prefix+float(np.logaddexp2.reduce(word[0]))-8*arity*groups*log2(255)
    return dict(surrogate_degree=16*arity*groups, survives_floor=9*arity*groups,
        majorant_component_log2=value, lower_bound_on_method_expression=True,
        failure_probability_lower_bound=False, RS_normalization_included=True,
        fixed_A_only=True, no_quiet_or_floor_tuning_can_remove_this_component=True)


def exact_pair_columns(source, table_r, *, groups=18):
    """Fresh seven-branch counts for each actual scheduled two-column source."""
    source = fiber._source(source)
    if len(source) != 2:
        raise ValueError('Exactly two original source labels required')
    h, ratios = fiber.outer.matrix(), fiber.outer.core.original.exceptional_ratios(188, 256)
    columns = []
    for j in range(16):
        before, after = column_permutations(fiber.SCHEDULE, j)
        a, b = (before[i] for i in source)
        masks = tuple(sum(bool(fiber.outer.core.original._multiply(row[a], ratio, 256) ^ row[b]) << i
                          for i, row in enumerate(h)) for ratio in ratios)
        branches = ((0, 0, 1), (1, 65535, 255), (1, 65535, 255),
                    *((2, mask, 255) for mask in masks), (2, 65535, 255*252))
        column = np.full((3, 2, 2), -np.inf)
        for weight, mask, count in branches:
            physical = _permute_mask(mask, after)
            column[weight] = np.logaddexp2(column[weight], log2(count)+table_r[physical]/groups)
        columns.append(column)
    return np.stack(columns)


def _jsonable(value):
    if isinstance(value, dict): return {k: _jsonable(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)): return [_jsonable(v) for v in value]
    if isinstance(value, np.generic): return _jsonable(value.item())
    if isinstance(value, float) and not isfinite(value):
        if value == -np.inf: return None
        raise ArithmeticError('Invalid nonfinite FLOAT evidence')
    return value


def run(*, mode='representatives', max_seconds=120, progress=None):
    if mode not in ('representatives', 'pairs'):
        raise ValueError('Use the eight representatives or complete two-row safety mode')
    if type(max_seconds) not in (int, float) or not isfinite(max_seconds) or not 0 < max_seconds <= 120:
        raise ValueError('Use a bounded wall budget at most120 seconds')
    pins, started = dependencies(), perf_counter()
    def check_time():
        if perf_counter()-started >= max_seconds:
            raise TimeoutError('Bounded new-matrix control exhausted its wall budget')
    saved, measure, prefix = load_measure()
    sources = SOURCES if mode == 'representatives' else tuple(combinations(range(16), 2))
    controls, preflight_rows, roots, error = [], [], [], None
    try:
        table_r = moment_table(measure, 18)
        table_2r = moment_table(measure, 36) if mode == 'representatives' else None
        for source in sources:
            check_time()
            a = len(source)
            if mode == 'representatives':
                columns, nonquiet, details = prepare_source(source, table_r, table_2r, check_time=check_time)
                diagnostic = all_nonquiet_diagnostic(nonquiet, a, prefix)
                preflight_rows.append(dict(original_source=source, arity=a, diagnostic=diagnostic,
                                           column_diagnostics=details))
                if progress:
                    progress(_jsonable(dict(stage='all-nonquiet-before-tail', original_source=source,
                        arity=a, diagnostic=diagnostic, elapsed_seconds=perf_counter()-started)))
            else:
                columns, details, diagnostic = exact_pair_columns(source, table_r), [], None
            tail = float(positive_deficit_tail(columns, 18, 9*a*18, check_time=check_time))
            root = tail/18-8*a*log2(255)
            row = dict(original_source=source, arity=a, status='complete',
                normalized_fixed_A_root_log2=root, fixed_A_repeated18_bound_log2=prefix+18*root,
                aggregate_input_floor=9*a*18, retained_deficit=7*a*18,
                all_nonquiet_diagnostic=diagnostic, column_diagnostics=details)
            roots.append(root)
            controls.append(row)
            if progress:
                progress(_jsonable(dict(stage='completed-fixed-A', original_source=source,
                    arity=a, bound_log2=row['fixed_A_repeated18_bound_log2'], elapsed_seconds=perf_counter()-started)))
    except TimeoutError as exception:
        error = str(exception)
    pins_after = dependencies()
    if pins_after != pins:
        raise ValueError('Sources changed during the new-construction control')
    if sha256((DIRECTORY/INPUT).read_bytes()).hexdigest() != INPUT_SHA:
        raise ValueError('Saved common regional measure changed during the control')
    complete = len(controls) == len(sources)
    pair_bound = prefix+18*float(np.logaddexp2.reduce(roots)) if mode == 'pairs' and complete else None
    return _jsonable(dict(construction=fiber.outer.CONSTRUCTION, supports=fiber.outer.SUPPORTS,
        lambda_aes_byte=188, schedule=fiber.SCHEDULE, mode=mode, groups=18,
        source_class='Fmax<=1', source_count=saved['source_count'], mark=saved['mark'],
        common_measure=measure, common_prefix_log2=prefix, input_artifact=INPUT, input_sha256=INPUT_SHA,
        exact_new_outer_fiber_counts=True, original_A_fixed_across16columns=True,
        protected_test_matrix='Hnew-I, including overlapping row sets',
        original_input_upper_tail_domination=True, coefficientwise_input_equality_claimed=False,
        all_nonquiet_diagnostic_is_method_lower_bound_not_failure=True,
        full_independent_D1=True, final_independent_D2_restores_NZ_amplitudes=True,
        source_and_projective_normalizations_once=True, RS_factor_per_group='255^(-8*a)',
        representative_extrapolation=False, complete_all120_original_pairs=mode == 'pairs' and complete,
        all_pair_rowset_sum_bound_log2=pair_bound, complete=complete, requested_sources=len(sources),
        completed_sources=len(controls), error=error, controls=controls,
        all_nonquiet_preflights=preflight_rows,
        reviewed_source_sha256=pins, source_pins_unchanged=True,
        regional_provider_runs=0, input_marker_searches=0, directed=False,
        old_outer_receipts_transferred=False, native_performance_claim=False, whole_code_certificate=False,
        elapsed_seconds=perf_counter()-started, max_seconds=max_seconds))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=('representatives', 'pairs'), required=True)
    parser.add_argument('--max-seconds', type=float, default=120)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists(): parser.error('Refuse to overwrite evidence')
    result = run(mode=args.mode, max_seconds=args.max_seconds,
                 progress=lambda row: print(json.dumps(row, allow_nan=False), flush=True))
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
