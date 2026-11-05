"""Bounded diagnosis of source-price versus continuation loss for two groups.

Four fixed marker choices compare a smaller rank price and a stronger output
tilt. Local envelopes are exactly verified and products round upward; printed
logarithms are diagnostics, not certificate comparisons. No file output.
"""

from decimal import Decimal as D, localcontext, ROUND_CEILING
from fractions import Fraction as F

from check_empty_span_potential import check
from check_endpoint_route import prepare, regions
from check_sparse_full_gate import multiply, upward
from check_two_groups import evaluate_fresh_columns


def logarithm(value):
    return float(value.ln() / D(2).ln()) if value else float("-inf")


def inspect(z, tau, contraction):
    vector, z, tau, contraction = check(z, contraction, tau=tau)
    prepared = prepare(vector, z, tau, length=4)
    minimum = min(vector[r][h] for h in range(1, 16) for r in range(1, h.bit_count() + 1))
    with localcontext() as context:
        context.prec = 80
        context.rounding = ROUND_CEILING
        converted = [[[[[upward(x) for x in row] for row in matrix] for matrix in tail]
                       for tail in kind] for kind in prepared]
        ordinary = regions(contraction, *converted[0])
        fresh = regions(contraction, *converted[1])
        bounds = evaluate_fresh_columns(*ordinary, *fresh, minimum, z, tau, contraction)
        print(f"DIAGNOSTIC log2 U_h: h33={logarithm(bounds[33]):.3f}, "
              f"h48={logarithm(bounds[48]):.3f}, h64={logarithm(bounds[64]):.3f}")
        print(f"h64 initial-rank cost={255 * logarithm(upward(tau)):.3f} bits; "
              f"output-threshold cost={-((7 * (1 << 21)) // 20) * logarithm(upward(z)):.3f} bits")
        for name, (a, b) in (("ordinary", ordinary), ("fresh", fresh)):
            product = multiply(a, b)
            column = multiply(product, product)
            print(f"{name} even-column log2 entries: "
                  f"qq={logarithm(column[0][0]):.3f}, ll={logarithm(column[1][1]):.3f}, "
                  f"ql*lq={logarithm(column[0][1] * column[1][0]):.3f}")
        return bounds


def main():
    candidates = ((F(499, 500), F(256), F(993, 1000)),
                  (F(499, 500), F(1024), F(993, 1000)),
                  (F(499, 500), F(65536), F(993, 1000)),
                  (F(999, 1000), F(65536), F(4981, 5000)))
    for z, tau, contraction in candidates:
        print(f"MARKERS z={z}, tau={tau}, lambda={contraction}", flush=True)
        try:
            inspect(z, tau, contraction)
        except RuntimeError as error:
            print(f"NO CLAIM: bounded proposal failed: {error}")
    print("Finite balance diagnostic only; no whole-code certificate or benchmark.")


if __name__ == "__main__":
    main()
