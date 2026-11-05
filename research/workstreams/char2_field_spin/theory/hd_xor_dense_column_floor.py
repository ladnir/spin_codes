"""One dense XOR2 group: positive transfer polynomial for nonempty columns.

Fix one active group with between3 and8 nonzero RS[32,16] constituents.
Every actual group output has at least17 nonempty eight-coordinate columns:
one nonzero RS row already meets17 columns, and D2 H D1 is invertible
separately in every column. This is a COLUMN floor, not an entry-weight floor.

The existing product domination bounds this group's expected message measure
by D_H(p) times iid f_p coordinates. Restrict that comparison to at least17
nonempty columns before evaluating the inner. Conditional on support, the
comparison's nonzero amplitudes are independent uniform GF256*. The regional
occupancy probability is eta=255p/256, on the group's one potential position.

Q0 and Q1 bound one chronological eight-region block, weighted respectively
by an empty and a nonempty column event under this comparison. Their positive
recurrence never subtracts one upper bound from another. Coefficients of
(Q0+Y Q1)^32 retain state continuity and count nonempty columns. Summing only
degrees17..32, taking the initial-zero row and terminal all-ones vector,
gives the marked comparison mass. Multiply by8192 D_H(p) z^-L /255 once.

This bounded evaluator is FLOAT only. No whole-code certificate is claimed.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from math import log2
from pathlib import Path
from time import perf_counter

import numpy as np

from hd_xor_column_potential import SAVED_WITNESSES
from hd_xor_multigroup import SAVED_ONE_ROW_128, Q, REGION_LENGTH, LIMIT
from hd_xor_mixed_gate import DENSE_P, dense_density_log2
from hd_xor_pair_route import pair_regions_screen
from hd_xor_paired_gate import CONSTRUCTION
from hd_shifted_swap_sparse import family_witness


R1_MARKS = (
    ('xor-one-row-r1', SAVED_ONE_ROW_128[0][2]),
    ('xor-two-row-r1', SAVED_WITNESSES[1][0]),
    ('swapped-sparse-r1', str(family_witness(1)[0])),
)


def log_identity():
    return np.array([[0., -np.inf], [-np.inf, 0.]])


def matrix_product(left, right):
    """Noncommuting2x2 multiplication of nonnegative matrices in log2 form."""
    return np.logaddexp2(left[..., :, 0, None]+right[..., None, 0, :],
                        left[..., :, 1, None]+right[..., None, 1, :])


def _checked_log_array(values, shape):
    values = np.asarray(values, dtype=float)
    if values.shape != shape or np.any(np.isnan(values)) or np.any(np.isposinf(values)):
        raise ValueError('Expected finite/negative-infinity log matrix entries with the declared shape')
    return values


def split_column(regions, eta):
    """Positive empty/nonempty-column split for labelled regional tables.

    Input shape is (positions, occupancy0-or1, row, column). A common table
    may be explicitly repeated by the caller; region order is not quotiented.
    The nonempty flag is Boolean OR, even if several positions are occupied.
    """
    regions = np.asarray(regions, dtype=float)
    if regions.ndim != 4 or len(regions) < 1:
        raise ValueError('At least one labelled regional table is required')
    regions = _checked_log_array(regions, (len(regions), 2, 2, 2))
    if not 0 <= eta <= 1:
        raise ValueError('Occupancy probability must lie in[0,1]')
    log_zero = -np.inf if eta == 1 else log2(1-eta)
    log_one = -np.inf if eta == 0 else log2(eta)
    empty, nonempty = log_identity(), np.full((2, 2), -np.inf)
    for table in regions:
        absent, occupied = table[0]+log_zero, table[1]+log_one
        nonempty = np.logaddexp2(matrix_product(nonempty, np.logaddexp2(absent, occupied)),
                                matrix_product(empty, occupied))
        empty = matrix_product(empty, absent)
    return empty, nonempty


def column_polynomial(empty, nonempty, columns=32):
    """Positive matrix coefficients of the ordered column-block product."""
    empty, nonempty = (_checked_log_array(matrix, (2, 2)) for matrix in (empty, nonempty))
    if type(columns) is not int or columns < 0:
        raise ValueError('A nonnegative integer column count is required')
    values = log_identity()[None, :, :]
    for _ in range(columns):
        following = np.full((len(values)+1, 2, 2), -np.inf)
        following[:-1] = matrix_product(values, empty)
        following[1:] = np.logaddexp2(following[1:], matrix_product(values, nonempty))
        values = following
    return values


def tail_log_mass(coefficients, minimum_columns):
    coefficients = np.asarray(coefficients, dtype=float)
    if (coefficients.ndim != 3 or coefficients.shape[1:] != (2, 2)
            or type(minimum_columns) is not int or not 0 <= minimum_columns < len(coefficients)):
        raise ValueError('Invalid matrix polynomial or nonempty-column floor')
    return float(np.logaddexp2.reduce(coefficients[minimum_columns:, 0, :].reshape(-1)))


def evaluate(mark, probability=DENSE_P):
    started = perf_counter()
    z, p = F(mark), F(probability)
    if not 0 < z < 1 or not 0 < p <= 1:
        raise ValueError('Rational marks must lie in their probability ranges')
    raw = pair_regions_screen(float(z), 1)
    regions = np.array([[raw[0, 0], raw[1, 0]]]*8)
    empty, nonempty = split_column(regions, float(p)*255/256)
    polynomial = column_polynomial(empty, nonempty)
    density = dense_density_log2(float(p))
    prefix = log2(REGION_LENGTH)+density-log2(Q-1)-LIMIT*log2(float(z))
    total, retained = tail_log_mass(polynomial, 0), tail_log_mass(polynomial, 17)
    by_degree = np.logaddexp2(polynomial[:, 0, 0], polynomial[:, 0, 1])
    return dict(mark=str(z), dense_probability=str(p), dense_group_density_log2=density,
        unrestricted_log2_screen=prefix+total, floor17_log2_screen=prefix+retained,
        retained_log2_mass=retained, unrestricted_log2_mass=total,
        all_zero_column_sequence_log2_mass=float(by_degree[0]),
        largest_retained_degree=int(17+np.argmax(by_degree[17:])),
        nonempty_column_degree_log2_masses=by_degree.tolist(),
        elapsed_seconds=perf_counter()-started)


def run():
    rows = [dict(witness_source=source, **evaluate(mark)) for source, mark in R1_MARKS]
    return dict(construction=CONSTRUCTION, active_groups=1, active_constituent_counts=[3, 4, 5, 6, 7, 8],
        source_supports=8192, source_matching_counts=[0, 0], field_order=256,
        forbidden_weight=LIMIT, nonempty_column_floor=17, total_columns=32,
        floor_unit='eight-coordinate-outer-columns', alpha=1,
        endpoint_potential_used=False, matrix_state_preserved_across_columns=True,
        positive_empty_nonempty_split=True, projective_normalizations=1,
        row_set_choices_inside_dense_density=True, directed=False,
        proposal_requires_fresh_directed_replay=True, whole_code_certificate=False,
        fitted_new_witnesses=0, entries=rows,
        best_log2_screen=min(row['floor17_log2_screen'] for row in rows))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite a previous diagnostic artifact')
    result = run()
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    with args.output.open('x', encoding='utf-8') as target:
        json.dump(result, target, indent=2, allow_nan=False)
        target.write('\n')
    print(json.dumps(result, allow_nan=False), flush=True)
