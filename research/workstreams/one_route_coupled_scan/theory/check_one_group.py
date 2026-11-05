"""All one-group support classes, using exact pair supports over the small field.

Rational-rank-two polynomial pairs have at most2(b-1) rank-one columns.
All final products are positive and rounded up. No file output or benchmarks.
This bound still excludes two or more active outer groups.
"""

from decimal import Decimal as D, ROUND_CEILING, localcontext
from fractions import Fraction as F

from check_empty_span_potential import check
from check_sparse_full_gate import active_envelope, multiply, power, upward


def column_matrices(regional, contraction):
    q = 256
    xi = [F(0)] + [F((q - 1)**k, (q - 1) * q**3) for k in range(1, 5)]
    theta = [F(0), F(0), F(1, (q - 1) * q**2), F(1, (q - 1) * q), F(1)]
    b0 = [[D(1), D(0)], [D(0), upward(contraction**2048)]]
    b1 = [[upward(x) for x in row] for row in regional]
    columns = []
    for mask in range(16):
        column = [[D(1), D(0)], [D(0), D(1)]]
        for a in range(4):
            column = multiply(column, b1 if mask & (1 << a) else b0)
        columns.append(column)
    active = []
    for weights in (xi, theta):
        active.append([[sum((upward(weights[mask.bit_count()]) * columns[mask][i][j]
                             for mask in range(1, 16)), D(0)) for j in range(2)] for i in range(2)])
    return columns[0], active[0], active[1]


def profile_moments(empty, rank_one, rank_two, terminal, columns=64):
    """[u^h v^r](empty+u(rank_two+v rank_one))^columns, starting quiet."""
    size = columns + 1
    rows = [[[D(0), D(0)] for _ in range(size)] for _ in range(size)]
    rows[0][0][0] = D(1)
    for column in range(columns):
        following = [[[D(0), D(0)] for _ in range(size)] for _ in range(size)]
        for h in range(column + 1):
            for r in range(h + 1):
                row = rows[h][r]
                for nh, nr, matrix in ((h, r, empty), (h + 1, r, rank_two), (h + 1, r + 1, rank_one)):
                    for j in range(2):
                        following[nh][nr][j] += row[0] * matrix[0][j] + row[1] * matrix[1][j]
        rows = following
    return [[row[0] + row[1] * terminal for row in degrees] for degrees in rows]


def rank_one_count(b, q=256):
    """Safe polynomial-pair count, allowing zeros at evaluation columns.

    Normalize the coprime scalar pair (a,b) up to a common nonzero scalar.
    Degree d has at most(q+1)q^(2d) possibilities; the vector factor has
    four polynomials of degree less than b-d. Extra coprimality/zero tests
    would only decrease this count.
    """
    return sum((q + 1) * q**(2 * d) * (q**(4 * (b - d)) - 1) for d in range(b))


def rank_two_count(b, r, q=256):
    """Two independently valid profile counts; use the smaller one.

    For a fixed projective line at each rank-one evaluation, the r row-pair
    constraints have full rank whenever their polynomial kernel contains
    two rationally independent rows. A dependency would give one common
    polynomial relation on every row by RS-dual interpolation.
    """
    if not 0 <= r <= 2 * b - 2:
        return 0
    s = min(r, b)
    interpolation = q**(8 * (b - s)) * ((q**4 - 1) * (q + 1))**s
    line_constraints = (q + 1)**r * q**(4 * (2 * b - r))
    return min(interpolation, line_constraints)


def evaluate(precision, vector, z, tau, contraction):
    if tau != 65536:
        raise ValueError("this exact-pair bridge uses tau=q^2")
    a, b, d, minimum = active_envelope(vector, z, tau)
    average = (1 - contraction**2048) / ((1 - contraction) * 2048)
    regional = [[F(0), a * average], [b * average, d * contraction**2047]]
    q = 256
    with localcontext() as context:
        context.prec = precision
        context.rounding = ROUND_CEILING
        matrices = column_matrices(regional, contraction)
        moments = profile_moments(*matrices, upward(1 / minimum))
        common = D(8192) * power(upward(1 / z), (7 * (1 << 21)) // 20) / D(q**2 - 256)
        bounds = {}
        components = {}
        for h in range(33, 65):
            b = h - 32
            r1 = D(rank_one_count(b)) * moments[h][h]
            r2 = D(0)
            for r in range(min(h, 2 * b - 2) + 1):
                count = rank_two_count(b, r, q)
                r2 += D(count) * moments[h][r]
            bounds[h] = common * (r1 + r2)
            components[h] = (common * r1, common * r2)
        return bounds, components


def main():
    candidates = ((F(999, 1000), F(4981, 5000)),
                  (F(7991, 8000), F(2489, 2500)),
                  (F(799, 800), F(9951, 10000)),
                  (F(7989, 8000), F(4973, 5000)))
    prepared = [check(z, contraction, tau=F(65536)) for z, contraction in candidates]
    for precision in (80, 120):
        results = [evaluate(precision, *item) for item in prepared]
        with localcontext() as context:
            context.prec = precision
            context.rounding = ROUND_CEILING
            best = {h: min(result[0][h] for result in results) for h in range(33, 65)}
            total = sum(best.values(), D(0))
            print(f"UPWARD precision={precision}: all one-group supports <= {total}")
            worst = max(best, key=best.get)
            print(f"worst h={worst}; h33={best[33]:.12E}; h64={best[64]:.12E}")
            if not total * D(1 << 147) < 1:
                raise AssertionError("all-one-group union does not close147 bits")
    print("PASS: 35%/147-bit union over EVERY one-group support and all8192 group choices.")
    print("Scope excludes two or more active groups; not a whole-code certificate.")


if __name__ == "__main__":
    main()
