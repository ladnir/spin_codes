"""Directed positive two-state envelope for ONE full-support outer group.

Uses the EXACT positive-vector check and ordinary singleton-input rank
transfer. All matrix/scalar products in the final calculation round up.
Also screens the 33-full-column class, which this envelope does NOT close.
It does not cover all support classes or certify the complete code.
No generated files; bounded small vectors/matrices, serial use only.
"""

from decimal import Decimal as D, localcontext, ROUND_CEILING
from fractions import Fraction as F
from math import log2

from check_empty_span_potential import check
from coupled_rank_transfer import local_matrix, pair_antipodal_cycle


def active_envelope(vector, z, tau):
    width, size = 4, 16
    cycle = pair_antipodal_cycle(width)
    valid = [(r, h) for r in range(1, width + 1)
             for h in range(1, size) if r <= h.bit_count()]
    activation = max(z * vector[1][1 << j] for j in range(width))
    reset = F(0)
    continuation = F(0)
    for j in range(width):
        # Global rank is priced at tau, but dependent roots still cost 1/q.
        matrix = local_matrix(width, 1 << j, z, min(tau, F(256)), cycle)
        for rank, state in valid:
            reset = max(reset, matrix[state][0] / vector[rank][state])
            total = sum((matrix[state][following]
                         * max(vector[r][following] for r in range(1, min(rank + 1, following.bit_count()) + 1))
                         for following in range(1, size)), F(0))
            continuation = max(continuation, total / vector[rank][state])
    minimum = min(vector[r][h] for r, h in valid)
    return activation, reset, continuation, minimum


def upward(value):
    return D(value.numerator) / D(value.denominator)


def power(value, exponent):
    """Only nonnegative multiplications under the caller's rounding mode."""
    result = D(1)
    while exponent:
        if exponent & 1:
            result *= value
        exponent >>= 1
        if exponent:
            value *= value
    return result


def multiply(a, b):
    return [[sum((a[i][k] * b[k][j] for k in range(2)), D(0))
             for j in range(2)] for i in range(2)]


def column_moments(regional, terminal, contraction):
    """All64 chronological column counts, under the caller's upward context."""
    one = [[upward(v) for v in row] for row in regional]
    one = multiply(one, one)
    one = multiply(one, one)
    empty = [[D(1), D(0)], [D(0), upward(contraction**8192)]]
    coefficients = [[D(0), D(0)] for _ in range(65)]
    coefficients[0][0] = D(1)
    for column in range(64):
        updated = [[D(0), D(0)] for _ in range(65)]
        for degree in range(column + 2):
            for j in range(2):
                updated[degree][j] = sum((coefficients[degree][k] * empty[k][j] for k in range(2)), D(0))
                if degree:
                    updated[degree][j] += sum((coefficients[degree - 1][k] * one[k][j] for k in range(2)), D(0))
        coefficients = updated
    return [row[0] + row[1] * upward(terminal) for row in coefficients]


def evaluate_columns(precision, regional, terminal, z, tau, contraction):
    """Failure upper bounds for h=33..64 full columns, D_h=4(h-32)."""
    with localcontext() as context:
        context.prec = precision
        context.rounding = ROUND_CEILING
        moments = column_moments(regional, terminal, contraction)
        threshold = (7 * (1 << 21)) // 20
        common = D(8192) * power(upward(1 / z), threshold)
        return {h: moments[h] * common * upward(tau**(4 * (h - 32) - 1)) for h in range(33, 65)}


def evaluate(precision, regional, terminal, z, tau, contraction):
    with localcontext() as context:
        context.prec = precision
        context.rounding = ROUND_CEILING
        matrix = [[upward(v) for v in row] for row in regional]
        # B^256: eight positive matrix squarings.
        for _ in range(8):
            matrix = multiply(matrix, matrix)
        value = matrix[0][0] + matrix[0][1] * upward(terminal)
        threshold = (7 * (1 << 21)) // 20
        value *= D(8192) * upward(tau**127)
        value *= power(upward(1 / z), threshold)
        full = +value
        # Exactly33 active columns, four active regions each. Chronological
        # coefficient DP counts column choices itself; no extra binomial.
        partial = column_moments(regional, terminal, contraction)[33]
        partial *= D(8192) * upward(tau**3)
        partial *= power(upward(1 / z), threshold)
        return full, +partial


def main():
    vector, z, tau, contraction = check(F(499, 500), F(993, 1000))
    activation, reset, continuation, minimum = active_envelope(vector, z, tau)
    steps = 2048
    geometric_average = (1 - contraction**steps) / ((1 - contraction) * steps)
    regional = [[F(0), activation * geometric_average],
                [reset * geometric_average, continuation * contraction**(steps - 1)]]
    print("Exact active envelope checked for all singleton coordinates and32 live(rank,support) states.")
    print("Diagnostic log2(a*b), log2(d):", log2(float(activation * reset)), log2(float(continuation)))
    for precision in (80, 120):
        full, partial = evaluate(precision, regional, 1 / minimum, z, tau, contraction)
        print(f"UPWARD precision={precision}: U_full <= {full}")
        print(f"UPWARD precision={precision}: U_33 <= {partial}")
        with localcontext() as context:
            context.prec = precision
            context.rounding = ROUND_CEILING
            if not full * D(1 << 849) < 1:
                raise AssertionError("This envelope does not close the 849-bit full-support gate")
            if not partial * D(1 << 40) < 1:
                print("OPEN: the 33-full-column envelope does not establish40 bits; not a codeword witness.")
    print("PASS: 35%/849-bit full-support gate, summing all8192 groups.")
    print("NOT a full-code certificate: partial columns and multiple groups remain.")


if __name__ == "__main__":
    main()
