import itertools
import math
import unittest
from fractions import Fraction

from cauchy_inner_transfer import count_table, floating_transfer, transfer


def mul(a, b, polynomial=0xb):
    result = 0
    top = 1 << (polynomial.bit_length()-1)
    while b:
        if b & 1:
            result ^= a
        b >>= 1
        a <<= 1
        if a & top:
            a ^= polynomial
    return result


def xor(values):
    out = 0
    for value in values:
        out ^= value
    return out


class CauchyTransferTests(unittest.TestCase):
    def test_full_small_field_enumeration(self):
        for q, polynomial, width, state in ((8, 0xb, 2, 1), (16, 0x13, 2, 2)):
            size = width+state
            inverses = [0]+[next(y for y in range(1, q) if mul(x, y, polynomial) == 1)
                            for x in range(1, q)]
            matrix = [[inverses[i ^ (size+j)] for j in range(size)] for i in range(size)]
            expected = count_table(width, state, q)
            for weight in range(width+1):
                for support in itertools.combinations(range(width), weight):
                    actual = [[[0]*(width+1) for _ in range(2)] for _ in range(2)]
                    for amplitudes in itertools.product(range(1, q), repeat=weight):
                        data = [0]*width
                        for position, value in zip(support, amplitudes):
                            data[position] = value
                        for registers in itertools.product(range(q), repeat=state):
                            source = data+list(registers)
                            result = [xor(mul(a, b, polynomial) for a, b in zip(row, source)) for row in matrix]
                            actual[int(any(registers))][int(any(result[width:]))][sum(x != 0 for x in result[:width])] += 1
                    self.assertEqual(tuple(tuple(map(tuple, row)) for row in actual), expected[weight])

    def test_mass_and_branch_distance(self):
        for width, state, q in ((8, 2, 256), (16, 4, 256), (16, 6, 256)):
            counts = count_table(width, state, q)
            for w, matrix in enumerate(counts):
                for incoming, row in enumerate(matrix):
                    self.assertEqual(sum(map(sum, row)), (q-1)**w*(q**state-1 if incoming else 1))
                    for outgoing, poly in enumerate(row):
                        for h, number in enumerate(poly):
                            if number and (w or incoming):
                                self.assertGreaterEqual(w+h+state*(incoming+outgoing), width+state+1)

    def test_exact_and_floating_evaluation(self):
        import numpy as np
        for z in (Fraction(0), Fraction(1, 9), Fraction(9, 10), Fraction(1)):
            exact = np.array(transfer(z), dtype=float)
            np.testing.assert_allclose(floating_transfer(float(z)), exact, rtol=2e-14, atol=0)
        self.assertEqual(transfer(1)[0][0], (1, 0))

    def test_cauchy_domain(self):
        for params in ((16, 4, 16), (0, 4, 256), (16, 0, 256), (8, 2, 255)):
            with self.assertRaises(ValueError):
                count_table(*params)

    def test_middle_screen_at_unit_output_tilt(self):
        from cauchy_field_screen import geometry, middle, support_corrections
        from field_rate_tradeoff import log_choose
        b, k, t, s, q, ell, active = 33, 24, 16, 4, 256, 16, 50
        x, rho = .001, .99
        rows, n, _ = geometry(ell, b, k, t)
        eta = max(value-w*math.log(rho) for w, value in support_corrections(b, k, q))
        expected = (-(b-1)*log_choose(rows, active)-(b-k)*active*math.log(q)
                    +active*eta-math.log(q-1)-b*active*math.log(x)
                    +n*math.log1p((1+q*rho)*x))
        self.assertAlmostEqual(middle(ell, active, x, 1., rho), -expected/math.log(2), places=8)


if __name__ == '__main__':
    unittest.main()
