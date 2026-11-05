"""Directed h=33 one-group gate including ALL partial scalar-column masks.

At tau=q^2, count ordered pairs over E in the SAME sampled outer. Their
rank is constant across these minimally supported RS words. This is not
an independent squaring of a first-moment outer bound. No generated files.
"""

from decimal import Decimal as D, ROUND_CEILING, localcontext
from fractions import Fraction as F

from check_empty_span_potential import check
from check_sparse_full_gate import active_envelope, multiply, power, upward


def coefficient_moment(empty, occupied, terminal, degree=33, column_count=64):
    """Coefficient of u^degree, with chronological noncommuting products."""
    rows = [[D(0), D(0)] for _ in range(degree + 1)]
    rows[0][0] = D(1)
    for column in range(column_count):
        following = [[D(0), D(0)] for _ in rows]
        for h in range(min(column + 1, degree) + 1):
            for j in range(2):
                following[h][j] = sum((rows[h][k] * empty[k][j] for k in range(2)), D(0))
                if h:
                    following[h][j] += sum((rows[h - 1][k] * occupied[k][j]
                                            for k in range(2)), D(0))
        rows = following
    return rows[degree][0] + rows[degree][1] * terminal


def evaluate(precision, regional, terminal, z, contraction):
    q = 256
    eta = [F(0)] + [min(F(1), F(q**k - 1, (q - 1) * q**3)) for k in range(1, 5)]
    theta = [F(0), F(0), F(1, (q - 1) * q**2), F(1, (q - 1) * q), F(1)]
    counts = [(q**4 - 1) * (q + 1), (q**4 - 1) * (q**4 - q)]
    with localcontext() as context:
        context.prec = precision
        context.rounding = ROUND_CEILING
        b0 = [[D(1), D(0)], [D(0), upward(contraction**2048)]]
        b1 = [[upward(x) for x in row] for row in regional]
        columns = []
        for mask in range(16):
            column = [[D(1), D(0)], [D(0), D(1)]]
            for a in range(4):
                column = multiply(column, b1 if mask & (1 << a) else b0)
            columns.append(column)
        terms = []
        for count, weights in zip(counts, (eta, theta)):
            occupied = [[sum((upward(weights[mask.bit_count()]) * columns[mask][i][j]
                              for mask in range(1, 16)), D(0)) for j in range(2)] for i in range(2)]
            moment = coefficient_moment(columns[0], occupied, upward(terminal))
            terms.append(moment * D(count) / D(q**2 - 1))
        common = D(8192) * power(upward(1 / z), (7 * (1 << 21)) // 20)
        return [term * common for term in terms]


def main():
    vector, z, tau, contraction = check(F(999, 1000), F(4981, 5000), tau=F(65536))
    activation, reset, continuation, minimum = active_envelope(vector, z, tau)
    average = (1 - contraction**2048) / ((1 - contraction) * 2048)
    regional = [[F(0), activation * average],
                [reset * average, continuation * contraction**2047]]
    for precision in (80, 120):
        terms = evaluate(precision, regional, 1 / minimum, z, contraction)
        with localcontext() as context:
            context.prec = precision
            context.rounding = ROUND_CEILING
            bound = sum(terms, D(0))
            print(f"UPWARD precision={precision}: all h33 masks <= {bound}")
            print(f"Rank-one pair term={terms[0]}; rank-two pair term={terms[1]}")
            if not bound * D(1 << 219) < 1:
                raise AssertionError("h33 partial-column union does not close 219 bits")
    print("PASS: 35%/219-bit union over all h33 column masks and all8192 groups.")
    print("NOT full-code: h>33 partial-column masks and multiple groups remain.")


if __name__ == "__main__":
    main()
