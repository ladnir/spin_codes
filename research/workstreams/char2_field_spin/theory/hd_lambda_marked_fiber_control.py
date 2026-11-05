"""Bounded FLOAT input-mark controls for fixed scheduled Hlambda row sets.

For each column, verify rank H[O,A]=rank V_A. Then zero outside output
forces the affine syndrome to vanish, so Hx=x on that fibre. Quiet output
masks therefore have their exact original input weight. For input mark
X>=1, assign input weight a to every nonquiet branch and use the unmarked
coset comparison only on that branch:

 K_A(X,t) <= Z_A(X*t_A) + X**a (C_A(t_O)-1) Z_A(t_A),  0<=t_i<=1.

This is an increasing-input-weight (upper-tail) bound, not coefficientwise
domination by input degree. This controller uses only scalar Chernoff
evaluation. It never applies the coset comparison to input-marked fibres.

Under one common positive entry/path measure, the nonquiet r-norm is
bounded both by the 2r factor norms and by a direct disjoint-mask sum.
Take their entrywise minimum. Quiet masks are summed with their exact
degree before the entrywise r-th power. Compose the sixteen resulting
column matrices in order, then charge X**(-9*a*r), RS255**(-8*a*r), and
the existing source/output-mark/projective prefix, each once.

Only five prescribed fixed-A diagnostics are implemented. There is no
row-set extrapolation, new output mark, new provider, or certificate.
"""
from hashlib import sha256
import json
from math import isfinite, log2
from pathlib import Path
from time import perf_counter

import numpy as np

from hd_bit_syndrome_outer import affine_rank
from hd_inverse_column_control import INPUT, INPUT_SHA
from hd_lambda_fiber_control import moment_table, support_norm
from hd_lambda_fiber_outer import column_fiber_counts
from hd_pair_scheduled_control import load_measure
from hd_post_refresh_quiet_rows import log_matmul


DIRECTORY = Path(__file__).resolve().parent
SOURCES = ((0, 1, 2, 3, 4), (0, 1, 2, 4, 8), tuple(range(7)),
           tuple(range(8)), (0, 1, 2, 3, 4, 5, 8, 12))
LOG_INPUT_GRID = (0., .125, .25, .5, 1., 2., 4., 8.)
MAX_REFINEMENT_ROUNDS = 3


class UnsupportedQuietFiber(ValueError):
    """The zero-outside fibre is not proved to be fixed pointwise by H."""


def _check_logs(values):
    values = np.asarray(values, dtype=float)
    if np.any(np.isnan(values)) or np.any(np.isposinf(values)):
        raise ValueError('Finite log values or exact zero (-inf) required')
    return values


def nonquiet_norm(data, table_r, table_2r, groups, *, check_time=None):
    """Return the minimum of two r-norm upper bounds for (C-1)Z.

    Direct union masks keep every labelled outside/inside choice. Their
    supports must be disjoint. Zero is removed from C by exact omission;
    Z retains its exact zero coefficient. No upper-bound subtraction occurs.
    """
    if type(groups) is not int or not 1 <= groups <= 18:
        raise ValueError('A group count from 1 through 18 is required')
    outside, zero = data['outside_counts'], data['zero_fiber_counts']
    if outside.get(0) != 1 or zero.get(0) != 1:
        raise ValueError('Both exact zero coefficients must equal one')
    for counts in (outside, zero):
        if any(type(mask) is not int or mask < 0 or type(n) is not int or n <= 0
               for mask, n in counts.items()):
            raise ValueError('Exact positive support counts required')
    omask = 0
    for mask in outside:
        omask |= mask
    imask = 0
    for mask in zero:
        imask |= mask
    if omask & imask:
        raise ValueError('Outside and inside support masks must be disjoint')
    table_r, table_2r = _check_logs(table_r), _check_logs(table_2r)
    if (table_r.ndim != 3 or table_r.shape[1:] != (2, 2)
            or table_2r.shape != table_r.shape
            or max(omask, imask, omask | imask) >= len(table_r)):
        raise ValueError('Compatible mask tables for r and 2r required')
    if check_time:
        check_time()
    factored = (support_norm(outside, table_2r, 2*groups, omit_zero=True)
                +support_norm(zero, table_2r, 2*groups))
    positive = tuple((mask, n) for mask, n in outside.items() if mask)
    if not positive:
        direct = np.full((2, 2), -np.inf)
        terms = 0
    else:
        if len(positive)*len(zero) > 65536:
            raise ValueError('Direct disjoint-mask expansion exceeds 65536 terms')
        outside_masks = np.array([mask for mask, _ in positive], dtype=np.uint32)
        inside_masks = np.array(tuple(zero), dtype=np.uint32)
        masks = (outside_masks[:, None] | inside_masks[None, :]).reshape(-1)
        logs_o = np.array([log2(n) for _, n in positive])
        logs_i = np.array([log2(n) for n in zero.values()])
        counts = (logs_o[:, None]+logs_i[None, :]).reshape(-1)
        direct = np.logaddexp2.reduce(table_r[masks]/groups+counts[:, None, None], axis=0)
        terms = len(masks)
    if check_time:
        check_time()
    return dict(selected=np.minimum(factored, direct), factor_2r=factored,
                direct_union=direct, direct_terms=terms)


def quiet_degree_norms(counts, arity, table_r, groups):
    """Bin quiet monomial r-norms by exact input weight.

    The caller must first prove that H fixes the zero-outside fibre.
    This function does not infer original-input weights from general fibres.
    """
    if type(arity) is not int or not 1 <= arity <= 16 or type(groups) is not int or groups < 1:
        raise ValueError('Positive bounded arity and group count required')
    table_r = _check_logs(table_r)
    result = np.full((arity+1, 2, 2), -np.inf)
    if counts.get(0) != 1:
        raise ValueError('The quiet zero word must have exact count one')
    for mask, count in counts.items():
        if (type(mask) is not int or not 0 <= mask < len(table_r)
                or mask.bit_count() > arity or type(count) is not int or count <= 0):
            raise ValueError('Exact quiet support counts within the arity are required')
        degree = mask.bit_count()
        result[degree] = np.logaddexp2(result[degree], log2(count)+table_r[mask]/groups)
    return result


def prepare_source(source, table_r, table_2r, *, groups=18, check_time=None):
    """Prepare marker-independent column coefficients; no numerical fit."""
    source = tuple(source)
    arity = len(source)
    columns, diagnostics = [], []
    for j in range(16):
        if check_time:
            check_time()
        data = column_fiber_counts(source, j)
        rank_v = affine_rank(data['mapped_source_mask'])
        if data['cross_rank'] != rank_v:
            raise UnsupportedQuietFiber(
                f'Column {j}: cross rank {data["cross_rank"]} != affine rank {rank_v}')
        bounds = nonquiet_norm(data, table_r, table_2r, groups, check_time=check_time)
        coefficients = quiet_degree_norms(data['zero_fiber_counts'], arity, table_r, groups)
        coefficients[arity] = np.logaddexp2(coefficients[arity], bounds['selected'])
        columns.append(coefficients)
        diagnostics.append(dict(column=j, mapped_source=list(data['mapped_source']),
            physical_inside_mask=data['physical_inside_mask'], cross_rank=data['cross_rank'],
            affine_rank=rank_v, quiet_fixed_point_identity_verified=True,
            outside_supports=len(data['outside_counts']), quiet_supports=len(data['zero_fiber_counts']),
            direct_terms=bounds['direct_terms'], nonquiet_norm_log2=bounds['selected'].tolist(),
            factored_norm_log2=bounds['factor_2r'].tolist(),
            direct_norm_log2=bounds['direct_union'].tolist()))
    return np.stack(columns), diagnostics


def word_log_moment(columns, log2_input_mark, groups, *, check_time=None):
    """Compose chronological column bounds after their entrywise r powers."""
    columns = _check_logs(columns)
    if (columns.ndim != 4 or columns.shape[2:] != (2, 2)
            or not 1 <= len(columns) <= 16 or columns.shape[1] < 2
            or type(groups) is not int or not 1 <= groups <= 18
            or not isfinite(log2_input_mark) or log2_input_mark < 0):
        raise ValueError('Ordered degree-binned column norms and a mark at least one required')
    degrees = np.arange(columns.shape[1], dtype=float)
    current = np.full((2, 2), -np.inf)
    current[0, 0] = current[1, 1] = 0.
    for column in columns:
        if check_time:
            check_time()
        norm = np.logaddexp2.reduce(column+degrees[:, None, None]*log2_input_mark, axis=0)
        current = log_matmul(current, groups*norm)
    return float(np.logaddexp2.reduce(current[0]))


def optimize_prepared(columns, arity, prefix, *, groups=18,
                      log_grid=LOG_INPUT_GRID, refinement_rounds=3, check_time=None):
    """Eight prescribed markers, then at most six adjacent midpoint trials.

    Each trial uses one common input mark across all columns of this fixed
    A. Per-A choices are deterministic proof witnesses, not setup-adaptive.
    The output mark and regional measure never change during this search.
    """
    grid = tuple(float(value) for value in log_grid)
    if (not grid or grid != tuple(sorted(set(grid))) or grid[0] != 0
            or any(not isfinite(value) or not 0 <= value <= 8 for value in grid)
            or len(grid) > 8 or type(refinement_rounds) is not int
            or not 0 <= refinement_rounds <= MAX_REFINEMENT_ROUNDS
            or type(arity) is not int or not 1 <= arity <= 16
            or np.asarray(columns).shape[1] != arity+1 or not isfinite(prefix)):
        raise ValueError('A bounded ordered mark grid and compatible arity are required')
    trials = {}
    floor = 9*arity*groups
    def evaluate(value):
        if value not in trials:
            if check_time:
                check_time()
            word = word_log_moment(columns, value, groups, check_time=check_time)
            after_floor = word-floor*value
            root = after_floor/groups-8*arity*log2(255)
            trials[value] = dict(log2_input_mark=value, word_log2_moment=word,
                chernoff_floor_charge_log2=-floor*value,
                normalized_fixed_A_root_log2=root,
                fixed_A_repeated_r_bound_log2=prefix+groups*root)
    for value in grid:
        evaluate(value)
    for _ in range(refinement_rounds):
        best = min(trials, key=lambda value: trials[value]['fixed_A_repeated_r_bound_log2'])
        ordered = sorted(trials)
        position = ordered.index(best)
        neighbors = ordered[max(0, position-1):position]+ordered[position+1:position+2]
        for other in neighbors:
            evaluate((best+other)/2)
    best = min(trials, key=lambda value: trials[value]['fixed_A_repeated_r_bound_log2'])
    return dict(best=trials[best], trials=[trials[value] for value in sorted(trials)],
                input_mark_search_bounds_log2=[0., grid[-1]],
                best_at_search_boundary=best in (0., grid[-1]),
                evaluations=len(trials), aggregate_input_floor=floor,
                rowset_sum_included=False, individual_input_row_floors_dropped=True)


def _jsonable(value):
    if isinstance(value, dict):
        return {key: _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    if isinstance(value, float) and not isfinite(value):
        if value == -np.inf:
            return None
        raise ArithmeticError('Invalid positive infinity or NaN in FLOAT evidence')
    return value


def run(*, sources=SOURCES, max_seconds=120, progress=None):
    """Run only the selected fixed-A controls, retaining partial timeout output."""
    if not isfinite(max_seconds) or not 0 < max_seconds <= 120:
        raise ValueError('A positive time limit at most 120 seconds is required')
    sources = tuple(tuple(source) for source in sources)
    if len(sources) > len(SOURCES) or any(source not in SOURCES for source in sources):
        raise ValueError('Only the five prescribed representative row sets are in scope')
    start = perf_counter()
    def check_time():
        if perf_counter()-start > max_seconds:
            raise TimeoutError('Bounded marked-fibre control exceeded its budget')
    saved, measure, prefix = load_measure()
    rows, complete, error = [], True, None
    try:
        check_time()
        table_r, table_2r = moment_table(measure, 18), moment_table(measure, 36)
        check_time()
        for source in sources:
            try:
                columns, diagnostics = prepare_source(source, table_r, table_2r, check_time=check_time)
                row = dict(original_source=list(source), arity=len(source), status='complete',
                           column_diagnostics=diagnostics,
                           **optimize_prepared(columns, len(source), prefix, check_time=check_time))
            except UnsupportedQuietFiber as exception:
                row = dict(original_source=list(source), arity=len(source),
                           status='unsupported_zero_fiber', error=str(exception))
            rows.append(row)
            if progress:
                progress(_jsonable(dict(**row, elapsed_seconds=perf_counter()-start)))
    except TimeoutError as exception:
        complete, error = False, str(exception)
    names = ('hd_lambda_marked_fiber_control.py', 'hd_lambda_fiber_control.py',
        'hd_lambda_fiber_outer.py', 'hd_lambda_outer.py', 'hd_outer_pair_layers.py',
        'hd_bit_syndrome_outer.py', 'hd_bit_syndrome_flat_geometry.py',
        'hd_pair_scheduled_outer.py', 'hd_pair_scheduled_control.py', 'hd_pair_scalar_control.py',
        'hd_inverse_column_control.py', 'hd_pair_identity_kernel.py', 'hd_pair_state_kernel.py',
        'hd_post_refresh_quiet_rows.py', 'hd_scheduled_column_control.py',
        'hd_pair_state_tail_general.py', 'hd_pair_state_tail.py')
    return _jsonable(dict(construction='scheduled-Hlambda-bc-no-cross', lambda_aes_byte=0xbc,
        schedule='input-inverse-output-xor', groups=18, source_class='Fmax<=1',
        mark=saved['mark'], source_count=saved['source_count'],
        input_artifact=INPUT, input_sha256=INPUT_SHA,
        common_measure=measure, common_prefix_log2=prefix, full_independent_D1=True,
        exact_quiet_input_degree_requires_rank_equality=True,
        nonquiet_surrogate_input_degree='arity', coefficientwise_input_domination_claimed=False,
        scalar_chernoff_input_mark=True, atom_masses_powered=False,
        nonquiet_norm_methods=['same-measure-2r-Holder', 'direct-disjoint-mask-Minkowski'],
        source_normalization_once=True, projective_normalization_once=True,
        representative_extrapolation=False, rowset_sum_included=False,
        complete=complete, error=error, directed=False, whole_code_certificate=False,
        regional_provider_runs=0, max_seconds=max_seconds,
        controls=rows, elapsed_seconds=perf_counter()-start,
        source_pins={name: sha256((DIRECTORY/name).read_bytes()).hexdigest() for name in names}))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--max-seconds', type=float, default=120)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Refuse to overwrite evidence')
    result = run(max_seconds=args.max_seconds,
                 progress=lambda row: print(json.dumps(row, allow_nan=False), flush=True))
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
