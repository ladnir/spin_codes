"""Actual two-group route averages with bounded endpoint-aware empty tails.

Positive matrix-series convolution retains short/zero gaps and batch
collisions. It never averages a nonlinear Bellman operator over routes.
Only shared full-column supports in two adjacent groups are covered.
"""

from decimal import Decimal as D, ROUND_CEILING, localcontext
from fractions import Fraction as F
from math import comb

from check_empty_span_potential import check
from check_endpoint_tail import endpoint_tails
from check_sparse_full_gate import multiply, upward
from check_two_groups import evaluate_fresh_columns


def zeros():
    return [[D(0), D(0)], [D(0), D(0)]]


def average_tails(vector, masks, z, tau, length, *, fresh=False, ranked=False, outer_randomized=False):
    """Exact per-mask tails, then their independent route-choice average."""
    tails = [endpoint_tails(vector, mask, z, tau, length, fresh=fresh, ranked=ranked,
                           outer_randomized=outer_randomized) for mask in masks]
    return [[[sum((t[length][i][j] for t in tails), F(0)) / len(masks)
              for j in range(2)] for i in range(2)] for length in range(length + 1)]


def gap_convolve(tails, previous, contraction):
    """Coefficients of K(x)previous(x) with geometric row tails after ell.

    K_L=diag(1,lambda**(L-ell)) K_ell for L>=ell.
    Arithmetic is positive under the caller's upward Decimal context.
    """
    ell = len(tails) - 1
    tail = zeros()
    result = []
    for n in range(len(previous)):
        if n >= ell:
            update = multiply(tails[ell], previous[n - ell])
            for i in range(2):
                factor = D(1) if i == 0 else contraction
                for j in range(2):
                    tail[i][j] = update[i][j] + factor * tail[i][j]
        value = [row[:] for row in tail]
        for gap in range(min(ell, n + 1)):
            term = multiply(tails[gap], previous[n - gap])
            for i in range(2):
                for j in range(2):
                    value[i][j] += term[i][j]
        result.append(value)
    return result


def regions(contraction, singleton, paired, cross, steps=2048):
    """Tail tables are Decimal; ordered event gaps sum to m-number_of_events."""
    contraction = upward(contraction)
    suffix = []
    live = D(1)
    for _ in range(steps):
        suffix.append([[D(1), D(0)], [D(0), live]])
        live *= contraction
    one_single = gap_convolve(singleton, suffix, contraction)
    two_single = gap_convolve(singleton, one_single, contraction)[steps - 2]
    pair = gap_convolve(paired, suffix, contraction)[steps - 1]
    collision = gap_convolve(cross, suffix, contraction)[steps - 1]
    split = zeros()
    for i in range(2):
        for j in range(2):
            pair[i][j] /= D(steps)
            split[i][j] = collision[i][j] / D(steps * (2 * steps - 1))
            split[i][j] += (two_single[i][j] * D(2 * steps - 2)
                            / D((2 * steps - 1) * comb(steps, 2)))
    return pair, split


def prepare(vector, z, tau, length=4, *, ranked=False):
    return [[average_tails(vector, masks, z, tau, length, fresh=fresh, ranked=ranked)
             for masks in ((1, 2, 4, 8), (3, 12), (5, 6, 9, 10))]
            for fresh in (False, True)]


def evaluate(precision, vector, z, tau, contraction, prepared):
    minimum = min(vector[r][h] for h in range(1, 16) for r in range(1, h.bit_count() + 1))
    with localcontext() as context:
        context.prec = precision
        context.rounding = ROUND_CEILING
        converted = [[[[[upward(x) for x in row] for row in matrix] for matrix in tail]
                       for tail in kind] for kind in prepared]
        pair, split = regions(contraction, *converted[0])
        fresh_pair, fresh_split = regions(contraction, *converted[1])
        return evaluate_fresh_columns(pair, split, fresh_pair, fresh_split,
                                      minimum, z, tau, contraction)


def main():
    vector, z, tau, contraction = check(F(999, 1000), F(4981, 5000), tau=F(65536))
    prepared = prepare(vector, z, tau)
    for precision in (80, 120):
        bounds = evaluate(precision, vector, z, tau, contraction, prepared)
        with localcontext() as context:
            context.prec = precision
            context.rounding = ROUND_CEILING
            total = sum(bounds.values(), D(0))
            print(f"UPWARD precision={precision}: endpoint-aware adjacent-group union <= {total}")
            print(f"h33={bounds[33]:.12E}; h64={bounds[64]:.12E}; worst h={max(bounds,key=bounds.get)}")
            closed = [h for h in bounds if bounds[h] * D(1 << 40) < 1]
            print(f"Separate h-classes below2^-40: {closed}")
            selected = sum((bounds[h] for h in (33, 34, 35)), D(0))
            if not selected * D(1 << 117) < 1:
                raise AssertionError("selected adjacent-group union does not close117 bits")
            print(f"PASS: h33..35 shared-full-column adjacent-group union <2^-117; bound={selected}")
    print("No whole-code claim. Actual short gaps, collisions, both parities and free terminal retained.")


if __name__ == "__main__":
    main()
