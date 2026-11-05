"""Exact local bounds for b(theta+a); no whole-code certificate or file output.

Run with python -B. Only rational arithmetic on dimension-16 expressions is
used. The model fixes the message and input-lane span before the fresh
refresh randomness. See CHEAP_LANE_MIXING_GOAL.md for proofs and limitations.
"""

from fractions import Fraction
from math import comb, log2
import json


def return_bound(h=16, *, q=256, state=16):
    """Return to U^state, from quotient support h>0, averaging D1 and a."""
    assert q >= 2 and 1 <= h <= state
    return sum((Fraction(comb(state, g), (q - 1) ** (state - g) * q**g)
                for g in range(state + 1 - h, state + 1)), Fraction(0))


def low_support_bound(j, *, q=256, state=16):
    """Pr[next quotient support <=j], conditional on current support >0."""
    assert 0 <= j <= state and q >= 2
    if j == 0:
        return min(Fraction(1), return_bound(q=q, state=state))
    eta = Fraction(1, q - 1) + Fraction(1, q)
    return min(Fraction(1), comb(state, j) * eta ** (state - j))


def exact(value):
    return {"numerator": str(value.numerator), "denominator": str(value.denominator),
            "approximation": float(value)}


def self_check():
    q, state = 256, 16
    eta = Fraction(1, q - 1) + Fraction(1, q)
    assert return_bound(1) == Fraction(1, q**state)
    assert return_bound() == eta**state - Fraction(1, (q - 1)**state)
    assert all(return_bound(h) <= return_bound(h + 1) for h in range(1, state))
    # Independent expansion of the target-set split into zero/outside W.
    for m in range(state + 1):
        expanded = sum((Fraction(comb(m, r), q**r * (q - 1)**(m-r))
                        for r in range(m + 1)), Fraction(0))
        assert expanded == eta**m
    tails = [low_support_bound(j) for j in range(state + 1)]
    assert all(Fraction(0) <= x <= 1 for x in tails)
    assert all(x <= y for x, y in zip(tails, tails[1:]))
    assert low_support_bound(12) < Fraction(7, 10**6)
    quiet_dim = 4096 * 128 - 32767 * 16
    assert quiet_dim == 16
    quiet_average = Fraction(1048576 * 255 * 256 ** (quiet_dim - 1), 256 ** quiet_dim - 1)
    assert quiet_average.numerator // quiet_average.denominator == 1044480
    # This is a floor on an UPPER-BOUND EXPRESSION, not bad-code probability.
    assert Fraction(257, 256)**256 > 2
    assert (16 * 32767) // 256 - 8 * 129 == 1015
    assert Fraction(256**160, 255**33 * 256**256) > Fraction(1, 256**129)


def main():
    self_check()
    extinction = return_bound()
    # 13 ceil(h/4)=52-13*(1[h<=12]+1[h<=8]+1[h<=4]+1[h=0]).
    expected_floor = 52 - 13 * sum((low_support_bound(j) for j in (12, 8, 4, 0)), Fraction(0))
    result = {
        "status": "exact_local_bounds_not_a_whole_code_certificate",
        "q": 256, "state": 16, "width": 64,
        "conditioning": "fixed message/span; nonzero incoming quotient state; fresh D1 and a",
        "return_probability_upper_bound": exact(extinction),
        "return_bound_bits_approx": -log2(float(extinction)),
        "support_at_most_12_upper_bound": exact(low_support_bound(12)),
        "next_output_expected_weight_lower_bound": exact(expected_floor),
        "next_output_expected_fraction_lower_bound": exact(expected_floor / 64),
        "fixed_stream_return_over_32767_transitions_upper_bound": exact(32767 * extinction),
        "scope": "The last bound stops at first return; no union over messages is included.",
        "deterministic_quiet_distance_ceiling": {"weight": 1044480, "length": 2097152,
                                                 "fraction": "255/512"},
        "raw_scalar_shadow_obstruction": {
            "scope": "floor on the proposed upper-bound expression, NOT bad-event probability",
            "single_group_log2_expression_floor_approx": (
                1280 - 33*log2(255) - 2048 + 16*32767*log2(257/256)),
            "rigorous_expression_lower_bound": "strictly greater than 2^1015",
            "valid_for": "all 0<alpha<=1 and 0<z<=1; raw uniform-reference mask operator",
        },
        "self_checks": "passed",
    }
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
