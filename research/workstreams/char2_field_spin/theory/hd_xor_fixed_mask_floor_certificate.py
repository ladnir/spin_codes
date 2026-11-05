"""Directed replay of the bounded r1 fixed-A mask/column-floor gate.

The counting-measure normalization and comparison law are described in
hd_xor_fixed_mask_floor.py and HD_XOR_FIXED_MASK_FLOOR.md. All mask counts
and beta factors are exact integers/rationals. This replays one fixed mark
for all219 A of sizes3..8; it neither fits nor imports FLOAT coefficients.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter

from flint import arb, ctx

from hd_shifted_coverage import rational_ceiling
from hd_xor_fixed_mask_floor import comparison_counts
from hd_xor_pair_route import arb_pair_regions
from hd_xor_paired_gate import CONSTRUCTION
from hd_xor_multigroup import Q, REGION_LENGTH, LIMIT


SAVED_MARK = '499451537/500000000'


def matrix_product_arb(left, right):
    a, b, c, d = left
    e, f, g, h = right
    return a*e+b*g, a*f+b*h, c*e+d*g, c*f+d*h


def _zero():
    return (arb(0),)*4


def _identity():
    return arb(1), arb(0), arb(0), arb(1)


def mask_operators_arb(regions):
    if not 2 <= len(regions) <= 8 or len(regions) % 2:
        raise ValueError('An even labelled group of size2..8 is required')
    if any(len(table) != 2 or any(len(matrix) != 4 for matrix in table) for table in regions):
        raise ValueError('Two flat-four matrices per labelled region required')
    if any(not entry.is_finite() or not entry >= 0 for table in regions for matrix in table for entry in matrix):
        raise ValueError('Nonnegative finite directed regional matrices required')
    result = []
    for mask in range(1 << len(regions)):
        value = _identity()
        for position, table in enumerate(regions):
            value = matrix_product_arb(value, table[(mask >> position) & 1])
        result.append(value)
    return result


def fixed_column_split_arb(operators, active_mask, group=8, q=256):
    if len(operators) != 1 << group:
        raise ValueError('A complete labelled mask-operator table is required')
    counts = comparison_counts(active_mask, group, q)
    normalizer = q**active_mask.bit_count()
    empty = tuple(entry/normalizer for entry in operators[0])
    nonempty = [arb(0) for _ in range(4)]
    for mask, count in enumerate(counts[1:], 1):
        if count:
            for index in range(4):
                nonempty[index] += count*operators[mask][index]
    return empty, tuple(entry/normalizer for entry in nonempty)


def column_polynomial_arb(empty, nonempty, columns=32):
    if type(columns) is not int or columns < 0 or len(empty) != 4 or len(nonempty) != 4:
        raise ValueError('Invalid matrices or column count')
    values = [_identity()]
    for _ in range(columns):
        following = [_zero() for _ in range(len(values)+1)]
        for degree, matrix in enumerate(values):
            first = matrix_product_arb(matrix, empty)
            second = matrix_product_arb(matrix, nonempty)
            following[degree] = tuple(a+b for a, b in zip(following[degree], first))
            following[degree+1] = tuple(a+b for a, b in zip(following[degree+1], second))
        values = following
    return values


def tail_mass_arb(coefficients, minimum_columns):
    if type(minimum_columns) is not int or not 0 <= minimum_columns < len(coefficients):
        raise ValueError('Invalid nonempty-column floor')
    return sum((entry[0]+entry[1] for entry in coefficients[minimum_columns:]), arb(0))


def receipt(precision=256, mark=SAVED_MARK):
    if type(precision) is not int or precision < 128:
        raise ValueError('At least128 precision bits required')
    z = F(mark)
    if not 0 < z < 1:
        raise ValueError('A rational mark in(0,1) is required')
    started = perf_counter()
    with ctx.workprec(precision):
        raw = arb_pair_regions(z, 1)
        operators = mask_operators_arb([(raw[0, 0], raw[1, 0])]*8)
        zz = arb(z.numerator)/z.denominator
        prefix = REGION_LENGTH*zz**(-LIMIT)/(Q-1)
        total = arb(0)
        arity_values = {active: arb(0) for active in range(3, 9)}
        rows = []
        for active_mask in range(1, 256):
            active = active_mask.bit_count()
            if active < 3:
                continue
            empty, nonempty = fixed_column_split_arb(operators, active_mask)
            coefficients = column_polynomial_arb(empty, nonempty)
            mass = tail_mass_arb(coefficients, 17)
            beta = arb(Q**(32*active))/(Q-1)**(16*active)
            value = prefix*beta*mass
            if not value > 0 or not value.is_finite():
                raise ArithmeticError('Invalid fixed-mask contribution')
            total += value
            arity_values[active] += value
            peak = max(range(17, 33), key=lambda degree: (coefficients[degree][0]+coefficients[degree][1]).upper())
            rows.append(dict(active_row_mask=active_mask, active_constituents=active,
                column_probability_normalizer=Q**active,
                failure_log2_upper=str(rational_ceiling(value.log()/arb(2).log())),
                largest_upper_mass_degree=peak))
        if len(rows) != 219:
            raise ArithmeticError('Incomplete fixed active-row-set coverage')
        arities = [dict(active_constituents=active,
                       fixed_row_sets=sum(row['active_constituents'] == active for row in rows),
                       failure_log2_upper=str(rational_ceiling(value.log()/arb(2).log())))
                   for active, value in arity_values.items()]
        return dict(construction=CONSTRUCTION, active_groups=1,
            active_constituent_counts=[3, 4, 5, 6, 7, 8], fixed_row_sets=219,
            exact_uniform_input_mask_mixture=True, uniform_input_includes_zero=True,
            comparison_is_actual_rs_law=False, beta_formula='q^(32a)/(q-1)^(16a)',
            conditional_nonzero_amplitude_operator=True, active_row_set_fixed_across_columns=True,
            source_supports=8192, source_matching_counts=[0, 0], field_order=Q,
            forbidden_weight=LIMIT, mark=str(z), nonempty_column_floor=17,
            total_columns=32, floor_unit='eight-coordinate-outer-columns', alpha=1,
            endpoint_potential_used=False, product_density_comparison_used=False,
            matrix_state_preserved_across_columns=True, positive_empty_nonempty_split=True,
            projective_normalizations=1, extra_row_set_binomial_factor=False,
            precision_bits=precision, directed=True, whole_code_certificate=False,
            class_reaches_40_bits=bool(total < arb(2)**-40), fitted_new_witnesses=0,
            failure_log2_upper=str(rational_ceiling(total.log()/arb(2).log())),
            arity_summaries=arities, entries=rows, elapsed_seconds=perf_counter()-started)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--precision', type=int, default=256)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite a previous directed receipt')
    result = receipt(args.precision)
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    with args.output.open('x', encoding='utf-8') as target:
        json.dump(result, target, indent=2, allow_nan=False)
        target.write('\n')
    print(json.dumps({key: value for key, value in result.items() if key != 'entries'}, allow_nan=False), flush=True)
