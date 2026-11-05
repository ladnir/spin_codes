"""Bounded FLOAT marked-fiber controls for scheduled Hlambda, arities 4..12.

Fix a mapped source set A, O its complement, and C=H[O,A]. On the zero-
outside fiber ker C, coordinate i is fixed by H exactly when row i of
H-I belongs to rowspace C. Test this by adjoining i to the cross-map rows.

Let I contain these fixed coordinates. For nonzero quiet output support T,
the original input support is contained in (T intersect I) union (A minus I).
The zero output has only the zero input, since H is invertible. Hence the
surrogate input degree is a-|I|+|T intersect I| for T!=0, and zero for T=0.
Physical masks below use Q_j(I), not I; original labels remain before P_j.

For X>=1 and output atoms in [0,1], the marked quiet support sum plus
X**a*(C_A-1)*Z_A is a valid pointwise upper bound. The nonquiet comparison
does NOT insert X into a translated-fiber inequality. Degree inflation gives
upper-tail domination, not coefficientwise domination at equal degree.

The controller reuses one saved common measure and scalar Chernoff only.
It never calls a regional provider or the bounded polynomial-tail backend.
Selected original A remain fixed across all sixteen columns. No row-set
union, representative extrapolation, directed arithmetic, or certificate.
"""
from hashlib import sha256
import json
from math import isfinite, log2
from pathlib import Path
from time import perf_counter

import numpy as np

from hd_inverse_column_control import INPUT, INPUT_SHA
from hd_lambda_fiber_outer import FULL, SCHEDULE, rank_cross_masks, column_fiber_counts
from hd_lambda_marked_fiber_control import (
    LOG_INPUT_GRID, _check_logs, _jsonable, nonquiet_norm, optimize_prepared)
from hd_lambda_fiber_control import moment_table
from hd_pair_scheduled_control import load_measure

DIRECTORY = Path(__file__).resolve().parent
MIN_ARITY, MAX_ARITY = 4, 12
SOURCES = (tuple(range(6)), (0, 1, 2, 4, 8, 15),
           tuple(range(7)), (0, 1, 2, 4, 8, 12, 15),
           tuple(range(8)), tuple(range(9)), tuple(range(10)),
           tuple(range(11)), tuple(range(12)))
SOURCE_NAMES = ('hd_lambda_partial_marked_fiber_control.py',
    'hd_lambda_marked_fiber_control.py', 'hd_lambda_fiber_control.py',
    'hd_lambda_fiber_outer.py', 'hd_lambda_outer.py', 'hd_outer_pair_layers.py',
    'hd_bit_syndrome_flat_geometry.py', 'hd_pair_scheduled_outer.py',
    'hd_pair_scheduled_control.py', 'hd_pair_scalar_control.py',
    'hd_inverse_column_control.py', 'hd_pair_identity_kernel.py',
    'hd_pair_state_kernel.py', 'hd_post_refresh_quiet_rows.py')


def _source(source):
    source = tuple(source)
    if (not MIN_ARITY <= len(source) <= MAX_ARITY or len(set(source)) != len(source)
            or any(type(i) is not int or not 0 <= i < 16 for i in source)):
        raise ValueError('This controller requires 4..12 distinct ordered original row labels')
    return source


def fixed_coordinates(mapped_source_mask):
    """Return exact fixed-coordinate mask BEFORE Q, using actual cross ranks."""
    if (type(mapped_source_mask) is not int or not 0 <= mapped_source_mask <= FULL
            or not MIN_ARITY <= mapped_source_mask.bit_count() <= MAX_ARITY):
        raise ValueError('An actual mapped source mask of arity4..12 is required')
    outside = FULL ^ mapped_source_mask
    rank = rank_cross_masks(outside, mapped_source_mask)
    fixed = sum(1 << i for i in range(16) if mapped_source_mask >> i & 1
                and rank_cross_masks(outside | (1 << i), mapped_source_mask) == rank)
    return fixed, rank


def quiet_surrogate_degree(mask, physical_inside, physical_fixed):
    """Cardinality upper bound on original input support for a quiet output."""
    if (any(type(m) is not int or not 0 <= m <= FULL
            for m in (mask, physical_inside, physical_fixed))
            or not MIN_ARITY <= physical_inside.bit_count() <= MAX_ARITY
            or physical_fixed & ~physical_inside or mask & ~physical_inside):
        raise ValueError('Physical quiet/fixed masks must lie inside an arity4..12 set')
    if not mask:
        return 0
    return (physical_inside ^ physical_fixed).bit_count() + (mask & physical_fixed).bit_count()


def column_geometry(data):
    """Check mapped/physical/original labels and transport the exact fixed set."""
    source = _source(data['original_source'])
    before, after = tuple(data['input_permutation']), tuple(data['output_permutation'])
    if any(len(p) != 16 or any(type(i) is not int for i in p)
           or sorted(p) != list(range(16)) for p in (before, after)):
        raise ValueError('Both schedules must be permutations of sixteen labels')
    mapped = tuple(before[i] for i in source)
    amask = sum(1 << i for i in mapped)
    inside = sum(1 << after[i] for i in mapped)
    if (tuple(data['mapped_source']) != mapped or data['mapped_source_mask'] != amask
            or data['physical_inside_mask'] != inside or data['arity'] != len(source)):
        raise ValueError('The original, mapped, and physical source labels disagree')
    fixed, rank = fixed_coordinates(amask)
    if data['cross_rank'] != rank:
        raise ValueError('Cached cross rank disagrees with the actual mapped source')
    physical_fixed = sum(1 << after[i] for i in range(16) if fixed >> i & 1)
    original_fixed = sum(1 << i for i in source if fixed >> before[i] & 1)
    return dict(arity=len(source), cross_rank=rank, mapped_source_mask=amask,
        mapped_fixed_mask=fixed, physical_inside_mask=inside,
        physical_fixed_mask=physical_fixed, original_fixed_mask=original_fixed,
        fixed_coordinates=fixed.bit_count(), quiet_pointwise_identity=(fixed == amask))


def quiet_degree_norms(data, table_r, groups):
    """Positive degree-binned quiet monomial r-norms, with safe inflation."""
    if type(groups) is not int or not 1 <= groups <= 18:
        raise ValueError('A group count from1 through18 is required')
    geometry = column_geometry(data)
    table_r = _check_logs(table_r)
    if table_r.shape != (65536, 2, 2):
        raise ValueError('The full sixteen-coordinate 2x2 mask-moment table is required')
    counts = data['zero_fiber_counts']
    if counts.get(0) != 1:
        raise ValueError('The quiet zero vector must have exact count one')
    result = np.full((geometry['arity']+1, 2, 2), -np.inf)
    for mask, count in counts.items():
        degree = quiet_surrogate_degree(mask, geometry['physical_inside_mask'],
                                        geometry['physical_fixed_mask'])
        if type(count) is not int or count <= 0:
            raise ValueError('Exact positive quiet support counts are required')
        result[degree] = np.logaddexp2(result[degree], log2(count)+table_r[mask]/groups)
    return result, geometry


def prepare_source(source, table_r, table_2r, *, groups=18, check_time=None):
    """Prepare all sixteen ordered columns without changing any frozen provider."""
    source = _source(source)
    columns, diagnostics = [], []
    for j in range(16):
        if check_time:
            check_time()
        data = column_fiber_counts(source, j, schedule=SCHEDULE, max_subsets=4096)
        coefficients, geometry = quiet_degree_norms(data, table_r, groups)
        bounds = nonquiet_norm(data, table_r, table_2r, groups, check_time=check_time)
        coefficients[len(source)] = np.logaddexp2(coefficients[len(source)], bounds['selected'])
        columns.append(coefficients)
        diagnostics.append(dict(column=j, mapped_source=list(data['mapped_source']),
            **geometry, zero_degree_preserved=True, nonquiet_surrogate_degree=len(source),
            outside_supports=len(data['outside_counts']), quiet_supports=len(data['zero_fiber_counts']),
            direct_terms=bounds['direct_terms'], nonquiet_norm_log2=bounds['selected'].tolist(),
            factored_norm_log2=bounds['factor_2r'].tolist(),
            direct_norm_log2=bounds['direct_union'].tolist()))
    return np.stack(columns), diagnostics


def run(sources=SOURCES, *, max_seconds=120, log_grid=LOG_INPUT_GRID, refinement_rounds=0, progress=None):
    """At most nine selected fixed-A scalar controls, never an all-A sweep."""
    sources = tuple(_source(source) for source in sources)
    if not 1 <= len(sources) <= 9 or len(set(sources)) != len(sources):
        raise ValueError('Select one through nine distinct ordered source tuples')
    if not isfinite(max_seconds) or not 0 < max_seconds <= 120:
        raise ValueError('Use a positive time limit at most120seconds')
    start = perf_counter()
    def check_time():
        if perf_counter()-start > max_seconds:
            raise TimeoutError('Bounded partial-coordinate control exceeded its budget')
    saved, measure, prefix = load_measure()
    rows, complete, error = [], True, None
    try:
        check_time()
        table_r, table_2r = moment_table(measure, 18), moment_table(measure, 36)
        check_time()
        for source in sources:
            columns, diagnostics = prepare_source(source, table_r, table_2r, check_time=check_time)
            row = dict(original_source=list(source), arity=len(source), status='complete',
                column_diagnostics=diagnostics,
                **optimize_prepared(columns, len(source), prefix, groups=18,
                    log_grid=log_grid, refinement_rounds=refinement_rounds, check_time=check_time))
            rows.append(row)
            if progress:
                progress(_jsonable(dict(**row, elapsed_seconds=perf_counter()-start)))
    except TimeoutError as exception:
        complete, error = False, str(exception)
    return _jsonable(dict(construction='scheduled-Hlambda-bc-no-cross', lambda_aes_byte=0xbc,
        schedule=SCHEDULE, groups=18, arity_scope=[4, 12], source_class='Fmax<=1',
        mark=saved['mark'], source_count=saved['source_count'], input_artifact=INPUT, input_sha256=INPUT_SHA,
        common_measure=measure, common_prefix_log2=prefix, full_independent_D1=True,
        original_A_fixed_across_columns=True, quiet_fixed_coordinates_verified_per_column=True,
        quiet_surrogate_degree='0 if T=0, otherwise a-|Ifix|+|T intersect physical_Ifix|',
        nonquiet_surrogate_input_degree='arity', coefficientwise_input_domination_claimed=False,
        upper_tail_domination=True, scalar_chernoff_input_mark=True,
        polynomial_tail_backend_used=False, atom_masses_powered=False,
        nonquiet_norm_methods=['same-measure-2r-Holder', 'direct-disjoint-mask-Minkowski'],
        source_normalization_once=True, projective_normalization_once=True,
        representative_extrapolation=False, rowset_sum_included=False,
        complete=complete, error=error, directed=False, whole_code_certificate=False,
        regional_provider_runs=0, max_seconds=max_seconds, controls=rows,
        elapsed_seconds=perf_counter()-start,
        source_pins={name: sha256((DIRECTORY/name).read_bytes()).hexdigest() for name in SOURCE_NAMES}))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', action='append', help='Comma-separated original row labels; repeat at most nine times; otherwise use the nine prescribed representatives')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--max-seconds', type=float, default=120)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Refuse to overwrite evidence')
    sources = tuple(tuple(int(i) for i in text.split(',')) for text in args.source) if args.source else SOURCES
    result = run(sources, max_seconds=args.max_seconds,
                 progress=lambda row: print(json.dumps(row, allow_nan=False), flush=True))
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
