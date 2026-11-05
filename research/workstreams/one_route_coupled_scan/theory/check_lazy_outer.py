"""Two adjacent full-column groups using unused outer Toeplitz coefficients.

First h-32 occupied columns use fixed invertible mixers and source freshness.
Each of the remaining32 columns per group is sampled/exposed lazily, with
one pivot partition per group. Conditioning back to invertible mixers costs
(256/255)^64. Directed positive arithmetic; not a whole-code certificate.
"""

from decimal import Decimal as D, ROUND_CEILING, localcontext
from fractions import Fraction as F

from check_empty_span_potential import check
from check_endpoint_route import average_tails, regions
from check_sparse_full_gate import multiply, upward
from check_two_groups import evaluate_fresh_columns
from lazy_outer_step import pivot_envelope


def prepare(vector, z, tau):
    masks = ((1, 2, 4, 8), (3, 12), (5, 6, 9, 10))
    lazy = [average_tails(vector, group, z, tau, 4, outer_randomized=True) for group in masks]
    fresh = [average_tails(vector, group, z, tau, 4, fresh=True, ranked=True) for group in masks]
    return lazy, fresh


def evaluate(precision, vector, z, tau, contraction, prepared):
    minimum = min(vector[r][h] for h in range(1, 16) for r in range(1, h.bit_count() + 1))
    pivot = pivot_envelope(vector, tau)
    with localcontext() as context:
        context.prec = precision
        context.rounding = ROUND_CEILING
        converted = [[[[[upward(x) for x in row] for row in matrix] for matrix in tail]
                       for tail in kind] for kind in prepared]
        ordinary = regions(contraction, *converted[0])
        fresh = regions(contraction, *converted[1])
        prefix = [[upward(x) for x in row] for row in pivot]
        # Sequential first-nonzero partitions for two groups; same global W.
        prefix = multiply(prefix, prefix)
        bounds = evaluate_fresh_columns(*ordinary, *fresh, minimum, z, tau, contraction,
                                        late_prefix=prefix)
        conditioning = upward(F(256, 255)**64)
        return {h: value * conditioning for h, value in bounds.items()}


def main():
    vector, z, tau, contraction = check(F(499, 500), F(993, 1000), tau=F(1024))
    prepared = prepare(vector, z, tau)
    for precision in (80, 120):
        bounds = evaluate(precision, vector, z, tau, contraction, prepared)
        with localcontext() as context:
            context.prec = precision
            context.rounding = ROUND_CEILING
            total = sum(bounds.values(), D(0))
            print(f"UPWARD precision={precision}: lazy-outer adjacent-group union <= {total}")
            print(f"h33={bounds[33]:.12E}; h48={bounds[48]:.12E}; h64={bounds[64]:.12E}")
            print(f"Separate h-classes below2^-40: {[h for h in bounds if bounds[h] * D(1 << 40) < 1]}")
    print("Scope remains shared full-column supports of adjacent groups, not the whole code.")


if __name__ == "__main__":
    main()
