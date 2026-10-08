"""Independent state diagonals around a fixed superregular mixer.

This is a new mixer, not an optimization of the rank-one distribution.
All joint spaces must be fixed before both current diagonal samplings.
No native implementation or whole-code distance claim is made.
"""

from fractions import Fraction as F
from itertools import combinations, product

from hd_blockwise_shared_maps import prescribed_mds_count
from hd_blockwise_shared_transfer import dot, multiply, power, rank


def cauchy_state_map(state=4, field=256):
    """Disjoint field labels 0..s-1 and s..2s-1; characteristic two."""
    if state < 1 or 2*state > field:
        raise ValueError('Need disjoint input/output field labels')
    return tuple(tuple(power(i ^ (state+j), field-2, field)
                       for j in range(state)) for i in range(state))


def point_probability(input_weight, target_weight, state=4, field=65536):
    """Exact Pr[D2 H D1 h=c] for fixed nonzero h,c and superregular H."""
    if not 1 <= input_weight <= state or not 1 <= target_weight <= state:
        raise ValueError('The point law concerns nonzero input and target')
    weight = input_weight+target_weight
    return F(prescribed_mds_count(weight, 2*state, state, field),
             (field-1)**weight)


def joint_return_bound(dimension, state=4, field=65536):
    """Cap once per fixed joint image, using at most its projective lines."""
    if not 0 <= dimension <= 2*state:
        raise ValueError('Invalid dimension in F^s x F^s')
    return min(F(1), F(field**dimension-1, (field-1)**(state+1)))


def small_mixers(field=4):
    """Bounded exhaustive diagonal sampler for s=2."""
    if field != 4:
        raise ValueError('This oracle intentionally uses only GF4')
    h = cauchy_state_map(2, field)
    for left in product(range(1, field), repeat=2):
        for right in product(range(1, field), repeat=2):
            yield tuple(tuple(multiply(left[i], multiply(h[i][j], right[j], field), field)
                              for j in range(2)) for i in range(2))


def k20_receipt(precision=128):
    from flint import arb, ctx

    previous = ctx.prec
    ctx.prec = precision
    try:
        p = joint_return_bound(2)
        bits = -(arb(p.numerator)/arb(p.denominator)).log()/arb(2).log()
        return dict(field=65536, state=4, joint_dimension=2,
                    exact_rational_upper=str(p), failure_bits_lower=str(bits),
                    fixed_point_upper=str(F(1, 65535**4)),
                    variable_field_products_per_packet=8,
                    additional_fixed_superregular_map=True,
                    joint_space_fixed_before_both_diagonals=True,
                    quiet_subcode_unchanged=True, whole_code_certificate=False)
    finally:
        ctx.prec = previous


if __name__ == '__main__':
    import json
    print(json.dumps(k20_receipt(), indent=2))
