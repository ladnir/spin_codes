"""Directed one-group bound for ALL choices of 33..64 fully populated columns.

Separates independent rank price tau=2^16 from coefficient root count q=256.
Conditional on any invertible group outer. Scalar holes inside occupied
columns and multiple groups are not covered. Source and small notes only.
"""

from decimal import Decimal as D, localcontext, ROUND_CEILING
from fractions import Fraction as F

from check_empty_span_potential import check
from check_sparse_full_gate import active_envelope, evaluate_columns


def main():
    vector, z, tau, contraction = check(F(499, 500), F(993, 1000), tau=F(65536))
    activation, reset, continuation, minimum = active_envelope(vector, z, tau)
    average = (1 - contraction**2048) / ((1 - contraction) * 2048)
    regional = [[F(0), activation * average],
                [reset * average, continuation * contraction**2047]]
    for precision in (80, 120):
        bounds = evaluate_columns(precision, regional, 1 / minimum, z, tau, contraction)
        with localcontext() as context:
            context.prec = precision
            context.rounding = ROUND_CEILING
            union = sum(bounds.values(), D(0))
            worst = max(bounds, key=bounds.get)
            print(f"UPWARD precision={precision}: all full-column supports union <= {union}")
            print(f"Worst column count={worst}; U33={bounds[33]}; U64={bounds[64]}")
            if not union * D(1 << 291) < 1:
                raise AssertionError("full-column union does not close 291 bits")
    print("PASS: 35%/291-bit union across all8192 groups and all33..64 full-column choices.")
    print("NOT full-code: partial occupied columns and multiple active groups remain.")


if __name__ == "__main__":
    main()
