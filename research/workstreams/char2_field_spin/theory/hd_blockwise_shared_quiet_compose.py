"""Base-field quiet composition: a high-type gate and a support-only obstruction.

All probabilities are over independent regional permutations. The high-type
class covers every active-row count, but the low-top-half-type complement is
not bounded here. No mixer or field-amplitude independence is assumed.
"""

from fractions import Fraction as F
from itertools import permutations, product
from math import ceil, comb, factorial

from hd_blockwise_shared_maps import prescribed_mds_count


def quiet_weight_polynomial(width=64, state=4, field=256):
    """Exact weight enumerator of the MDS feedback kernel."""
    if not 0 <= state < width < field+1:
        raise ValueError('Invalid MDS feedback-kernel parameters')
    return tuple(comb(width, w)*prescribed_mds_count(w, width, width-state, field)
                 for w in range(width+1))


def orbit_size(counts):
    """Number of ordered columns with the specified value multiplicities."""
    if any(count < 0 for count in counts):
        raise ValueError('Negative multiplicity')
    result = factorial(sum(counts))
    for count in counts:
        result //= factorial(count)
    return result


def worst_capped_counts(rows, maximum, alphabet=256):
    """Most concentrated feasible type when every multiplicity is <=maximum."""
    if rows < 1 or maximum < 1 or maximum*alphabet < rows:
        raise ValueError('The maximum multiplicity cannot accommodate the column')
    full, remainder = divmod(rows, maximum)
    return (maximum,)*full+((remainder,) if remainder else ())


def quiet_ball_log2(precision=160):
    """Directed projective quiet-output ball count at strict40% of N=2^21."""
    from flint import arb, ctx

    previous = ctx.prec
    ctx.prec = precision
    try:
        field, width, packets = 256, 64, 32768
        limit = (2*width*packets-1)//5
        z = F(2, 765)
        polynomial = quiet_weight_polynomial(width, 4, field)
        value = sum((coefficient*z**w for w, coefficient in enumerate(polynomial)), F(0))
        last = (1+(field-1)*z)**width
        a = arb(value.numerator)/arb(value.denominator)
        b = arb(last.numerator)/arb(last.denominator)
        az = arb(z.numerator)/arb(z.denominator)
        # Counting zero as well is a harmless overcount before division by q-1.
        upper = ((packets-1)*a.log()+b.log()-limit*az.log()-arb(field-1).log())/arb(2).log()
        return upper
    finally:
        ctx.prec = previous


def high_type_receipt(maximum=47, failure_bits=80, precision=160):
    """Directed all-density class bound and the actual top-half-type threshold."""
    from flint import arb, ctx

    if failure_bits < 1:
        raise ValueError('Need a positive requested failure exponent')
    previous = ctx.prec
    ctx.prec = precision
    try:
        rows, parity_regions = 8192, 128
        ball = quiet_ball_log2(precision)
        # A floating conversion only proposes an integer; the next assertion
        # checks the final inequality with directed arithmetic.
        threshold = ceil(float(ball)+failure_bits)
        if not arb(threshold)-ball > failure_bits:
            raise ArithmeticError('The rounded type threshold was not certified')
        counts = worst_capped_counts(rows, maximum)
        information = arb(orbit_size(counts)).log()/arb(2).log()
        upper = ball-parity_regions*information
        return dict(field=256, rows=rows, outer=[256, 128], packet_width=64, state=4,
                    output_weight_limit=838860, rational_weight_mark='2/765',
                    quiet_projective_ball_log2_upper=str(ball),
                    requested_failure_bits=failure_bits,
                    sufficient_top128_information_bits=threshold,
                    threshold_failure_log2_upper=str(ball-threshold),
                    average_information_per_region_position=str(F(threshold, parity_regions*rows)),
                    sufficient_balanced_region_count=parity_regions,
                    maximum_value_multiplicity=maximum,
                    worst_allowed_column_information_bits=str(information),
                    bounded_multiplicity_class_log2_upper=str(upper),
                    bounded_multiplicity_class_below_2_to_minus_40=bool(upper < -40),
                    covers_all_active_row_counts=True,
                    low_top_half_type_complement_not_covered=True,
                    nonquiet_words_not_covered=True, whole_code_certificate=False)
    finally:
        ctx.prec = previous


def dense_input_sparse_output_witness(translation=64):
    """Local support-relaxation counterexample; outer value checks are NOT imposed."""
    from hd_blockwise_maps import matvec
    from hd_blockwise_shared_maps import matrices, rank

    feedthrough, _, feedback, _, _ = matrices(translation)
    support = tuple(range(0, 64, 8))
    output = tuple(int(i in support) for i in range(64))
    source = matvec(feedthrough, output)
    restricted = tuple(tuple(row[j] for j in support) for row in feedback)
    if (rank(restricted) != 4 or any(matvec(feedback, source))
            or any(matvec(feedback, output)) or not all(source)
            or matvec(feedthrough, source) != output):
        raise ArithmeticError('The beta-specific dense/sparse quiet witness failed')
    count = prescribed_mds_count(8, 8, 4, 256)
    return dict(translation=translation, input=source, output=output,
                input_weight=64, output_weight=8, feedback_zero=True,
                restricted_feedback_rank=4, quiet_support_space_dimension=4,
                full_output_support_vectors_per_packet=count,
                active_source_rows_after_arbitrary_region_route=8192,
                global_output_weight_if_repeated=262144,
                outer_value_constraints_not_imposed=True,
                actual_code_counterexample=False)


def small_information_set_check():
    """Exhaustive GF4 [3,2] outer, two rows, all fixed inverse columns/routes."""
    from hd_blockwise_inner import multiply

    routes = tuple(permutations(range(2)))
    checked, constant_membership_one, inconsistent_types = 0, 0, 0
    for values in product(range(4), repeat=6):
        columns = tuple(tuple(values[2*c:2*c+2]) for c in range(3))
        probabilities = [F(1) if column[0] == column[1] else F(1, 2) for column in columns]
        bound = min(probabilities)
        hits = 0
        for choices in product(routes, repeat=3):
            ordered = tuple(tuple(column[i] for i in order) for column, order in zip(columns, choices))
            valid = True
            for row in range(2):
                constant = ordered[0][row]
                linear = constant ^ ordered[1][row]
                if ordered[2][row] != constant ^ multiply(2, linear, 4):
                    valid = False
                    break
            hits += int(valid)
        actual = F(hits, 8)
        if actual > bound:
            raise ArithmeticError('The information-set type bound failed')
        constant_membership_one += int(actual == 1)
        inconsistent_types += int(actual == 0)
        checked += 1
    return dict(field=4, outer=[3, 2], rows=2, fixed_inverse_words_checked=checked,
                routes_per_word=8, words_with_membership_probability_one=constant_membership_one,
                words_with_zero_membership_probability=inconsistent_types)


def main():
    import json
    print(json.dumps(high_type_receipt()))
    print(json.dumps(dense_input_sparse_output_witness()))
    print(json.dumps(small_information_set_check()))


if __name__ == '__main__':
    main()
