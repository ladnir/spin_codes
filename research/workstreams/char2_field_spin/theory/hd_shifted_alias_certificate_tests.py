"""Exact-reference checks for the directed shifted-pair alias replay."""
from fractions import Fraction as F
from itertools import product
import unittest
from unittest.mock import patch

import numpy as np
from flint import arb, ctx

from hd_shifted_alias_certificate import (
    _add, _mul, _sub, binary64_enclosure, interval_fft,
    interval_matrix_power, odd_coefficients, positive_alias_power,
    packet_table, rational, selected_inverse,
)


class AliasCertificateTests(unittest.TestCase):
    def test_real_interval_arithmetic_against_exact_fractions(self):
        points = [0., -1., 1., .1, np.nextafter(0., 1.), 2.**-500, 2.**500]
        for a, b in product(points, repeat=2):
            for operation, exact in ((_add, F(float(a))+F(float(b))),
                                     (_sub, F(float(a))-F(float(b))),
                                     (_mul, F(float(a))*F(float(b)))):
                lo, hi = operation((a, a), (b, b))
                self.assertLessEqual(F(float(lo)), exact)
                self.assertLessEqual(exact, F(float(hi)))

    def test_arb_to_binary64_checks_input_uncertainty(self):
        with ctx.workprec(256):
            for value in (arb(0), rational(F(1, 3)), rational(F(1, 10))**rational(F(3, 7)),
                          arb('1 +/- 0.01'), arb(2)**-1100):
                lo, hi = binary64_enclosure(value)
                self.assertTrue(arb(lo) <= value.lower())
                self.assertTrue(value.upper() <= arb(hi))

    def test_interval_fft_inverse_contains_exact_inputs(self):
        values = (np.arange(4*8**3, dtype=float).reshape(4, 8, 8, 8)%19-9)/16
        zero = np.zeros_like(values)
        interval = ((values.copy(), values.copy()), (zero.copy(), zero.copy()))
        with ctx.workprec(192):
            transformed = interval_fft(interval, (1, 2, 3))
            restored = interval_fft(transformed, (1, 2, 3), inverse=True)
        self.assertTrue(np.all(restored[0][0] <= values*512))
        self.assertTrue(np.all(values*512 <= restored[0][1]))
        self.assertTrue(np.all(restored[1][0] <= 0))
        self.assertTrue(np.all(0 <= restored[1][1]))

    def test_alias_contains_exact_cyclic_matrix_polynomial_power(self):
        grid, width = 4, 2
        polynomial = np.zeros((4, width, width, width), dtype=np.int64)
        # These positive coefficient matrices do not commute. The exact
        # cyclic reference must therefore retain matrix multiplication order.
        polynomial[:, 0, 0, 0] = (1, 2, 3, 1)
        polynomial[:, 1, 0, 1] = (2, 0, 1, 3)
        polynomial[:, 0, 1, 0] = (0, 5, 2, 1)
        base = {index: polynomial[(slice(None),)+index].reshape(2, 2).astype(object)
                for index in product(range(width), repeat=3)
                if polynomial[(slice(None),)+index].any()}
        result = {(0, 0, 0): np.eye(2, dtype=object)}
        for _ in range(8):
            following = {}
            for left_degree, left in result.items():
                for right_degree, right in base.items():
                    degree = tuple((a+b)%grid for a, b in zip(left_degree, right_degree))
                    following[degree] = following.get(degree, np.zeros((2, 2), dtype=object))+left@right
            result = following
        coefficients = [[arb(int(value)) for value in row.ravel()] for row in polynomial]
        with ctx.workprec(192):
            for target in ((0, 0, 0), (1, 2, 1), (2, 3, 2), (3, 1, 3)):
                bound, diagnostics = positive_alias_power(coefficients, target, grid, steps=8)
                expected = result.get(target, np.zeros((2, 2), dtype=object)).ravel()
                for actual, exact in zip(bound, expected):
                    self.assertTrue(actual >= arb(int(exact)))
                    self.assertLess(float(actual-arb(int(exact))), max(1., int(exact))*1e-9)
                self.assertTrue(diagnostics['directed'])

    def test_selected_inverse_matches_full_inverse_enclosures(self):
        values = np.arange(4*4**3, dtype=float).reshape(4, 4, 4, 4)%7
        zero = np.zeros_like(values)
        interval = ((values.copy(), values.copy()), (zero.copy(), zero.copy()))
        with ctx.workprec(192):
            frequency = interval_matrix_power(interval_fft(interval, (1, 2, 3)), 2)
            full = interval_fft(frequency, (1, 2, 3), inverse=True)
            target = (1, 2, 3)
            selected = selected_inverse(frequency, target)
            for part in range(2):
                lo = full[part][0][(slice(None),)+target]/64
                hi = full[part][1][(slice(None),)+target]/64
                self.assertTrue(np.all(selected[part][0] <= hi))
                self.assertTrue(np.all(lo <= selected[part][1]))

    def test_exact_half_profile_multiplicities_and_marker_axes(self):
        class ConstantPackets(dict):
            def items(self):
                from itertools import combinations_with_replacement
                return ((occupancy, (F(1),)*4)
                        for occupancy in combinations_with_replacement(range(17), 4))
        # With constant packets and alpha=1, the polynomial is exactly
        # (1 + x/8 + y/8 + u/16)^32 in every matrix entry.
        with ctx.workprec(192), patch('hd_shifted_alias_certificate.packet_table', return_value=ConstantPackets()):
            coefficients = odd_coefficients(F(1, 2), F(1), F(1, 8), F(1, 16))
        from math import factorial
        for a, b, c in ((0, 0, 0), (1, 2, 3), (7, 9, 16), (32, 0, 0), (0, 32, 0), (0, 0, 32)):
            count = factorial(32)//(factorial(32-a-b-c)*factorial(a)*factorial(b)*factorial(c))
            expected = rational(F(count, 8**(a+b)*16**c))
            for row in coefficients:
                self.assertTrue(row[a*33**2+b*33+c].contains(expected))

    def test_invalid_polynomial_and_grid_rejected(self):
        coefficients = [[arb(1)] for _ in range(4)]
        with self.assertRaises(ValueError):
            positive_alias_power(coefficients, (0, 0, 0), 3)
        coefficients[0][0] = arb(-1)
        with self.assertRaises(ValueError):
            positive_alias_power(coefficients, (0, 0, 0), 2)

    def test_packet_rounding_dominates_exact_rational_caps(self):
        from hd_common_refresh_packets import local_moments, matrix
        from hd_grouped_quiet import quiet_cap
        z = F(1, 2)
        local = local_moments(z, q=256, h=4)
        table = packet_table(z)
        density = F(256**16-1, 255**16)
        for occupancy in ((0, 0, 0, 0), (0, 1, 4, 16), (4, 5, 9, 13), (16, 16, 16, 16)):
            k, birth, ret, live = matrix(occupancy, z, local, q=256, h=4)
            exact = (k, density*birth, ret,
                     min(density*live, quiet_cap(z, q=256, h=4)) if not sum(occupancy) else density*live)
            for rounded, entry in zip(table[occupancy], exact):
                self.assertGreaterEqual(rounded, entry)
                self.assertLessEqual(rounded, entry*(1+F(1, 2**58)))


if __name__ == '__main__':
    unittest.main()
