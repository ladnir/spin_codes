"""Paired active-rank refinement composed with endpoint-aware actual routes.

Only two cyclically adjacent groups with the same fully populated columns.
No field-uniform assumption, independent insertion approximation, terminal
reset, or whole-code claim. Positive final arithmetic rounds upward.
"""

from decimal import Decimal as D, ROUND_CEILING, localcontext
from fractions import Fraction as F

from check_empty_span_potential import check
from check_endpoint_route import prepare, evaluate


def main():
    vector, z, tau, contraction = check(F(499, 500), F(993, 1000), tau=F(1024))
    prepared = prepare(vector, z, tau, ranked=True)
    for precision in (80, 120):
        bounds = evaluate(precision, vector, z, tau, contraction, prepared)
        with localcontext() as context:
            context.prec = precision
            context.rounding = ROUND_CEILING
            total = sum(bounds.values(), D(0))
            print(f"UPWARD precision={precision}: paired-rank shared-column union <= {total}")
            print(f"h33={bounds[33]:.12E}; h48={bounds[48]:.12E}; h64={bounds[64]:.12E}")
            print(f"Separate h-classes below2^-40: {[h for h in bounds if bounds[h] * D(1 << 40) < 1]}")
    print("Not a whole-code certificate: distinct column sets, partial columns and other group sets remain.")


if __name__ == "__main__":
    main()
