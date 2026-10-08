"""Exact input-mark inequalities and directed derivative checks."""

from fractions import Fraction as F
from itertools import product
from math import isfinite
import unittest

from hd_randomized_outer_joint_quiet import joint_moment_enclosures
from hd_randomized_outer_quiet_exact import regional_coefficients
from hd_randomized_outer_quiet_support import (
    base_nonzero_logs, directed_endpoint, directed_input_derivative,
    weighted_moment_enclosures,
)


class QuietSupportTests(unittest.TestCase):
    def setUp(self):
        from flint import ctx
        self.previous = ctx.prec
        ctx.prec = 384

    def tearDown(self):
        from flint import ctx
        ctx.prec = self.previous

    def test_packet_surrogate_and_support_mark_by_exact_enumeration(self):
        # GF4 addition suffices for these linear maps. Three potential slots,
        # one feedback coordinate, and an invertible binary feedthrough.
        q, j, z = 4, 3, F(1, 3)
        records = []
        for x in product(range(q), repeat=j):
            if x[0] ^ x[1] ^ x[2]:
                continue
            y = (x[0] ^ x[1], x[1] ^ x[2], x[0] ^ x[1] ^ x[2])
            records.append((sum(t != 0 for t in x), z**sum(t != 0 for t in y)))
        base = sum((weight for unused, weight in records), F(0))/q**j
        for v in (F(1), F(3, 2), F(2)):
            moment = sum((v**w*weight for w, weight in records), F(0))/q**j
            surrogate = F(1, q**j)+v**j*(base-F(1, q**j))
            self.assertLessEqual(moment, surrogate)
            for floor in range(4):
                restricted = sum((weight for w, weight in records if w >= floor), F(0))/q**j
                self.assertLessEqual(restricted, moment/v**floor)

    def test_unit_input_mark_and_exact_zero_input_packets(self):
        from flint import arb
        z = F(1, 256)
        base = joint_moment_enclosures(z)
        marked = weighted_moment_enclosures(z, F(0))
        for original, new in zip(base, marked):
            self.assertTrue(all((a-b).contains(0) for a, b in zip(original, new)))
        quiet, free = weighted_moment_enclosures(z, F(8))
        for j in range(5):
            self.assertTrue((quiet[j]-arb(1)/65536**j).contains(0))
        expected = (1+arb(65535)*256*(arb(1)/256)**8)/65536
        self.assertTrue((free[1]-expected).contains(0))

    def test_small_nonzero_term_survives_arb_subtraction(self):
        from flint import arb
        z = F(1, 65536)
        quiet, unused = joint_moment_enclosures(z)
        zero = arb(1)/65536**5
        nonzero = quiet[5]-zero
        self.assertTrue(nonzero > 0)
        self.assertEqual(float(quiet[5]), float(zero))
        logs, unused = base_nonzero_logs(z)
        self.assertTrue(isfinite(logs[5]))
        marked, unused = weighted_moment_enclosures(z, F(32))
        self.assertTrue((marked[5]-(zero+arb(2)**160*nonzero)).contains(0))
        self.assertTrue(marked[5] > zero)

    def test_exact_coefficient_log_convexity_of_the_surrogate(self):
        from flint import arb
        values = []
        for v in (1, 2, 4):
            quiet = (F(1), F(1, 4)+F(1, 16)*v, F(1, 16)+F(1, 32)*v**2)
            free = (F(1), F(1, 4)+F(1, 8)*v, F(1, 16)+F(1, 8)*v**2)
            region = regional_coefficients(quiet, free, F(1, 2), F(1, 3), 3, 2)
            values.append(region['quiet']**2*region['last']/arb(v)**3)
        self.assertTrue(values[1]**2 < values[0]*values[2])

    def test_directed_derivative_matches_a_small_forward_difference(self):
        from flint import arb
        args = (2048, F(1, 4))
        derivative = directed_input_derivative(*args, F(1, 32))
        self.assertTrue(derivative['no_input_mark_improves_this_witness'])
        expected = arb(derivative['derivative_at_log2_v_zero'])
        self.assertTrue(expected > 8125)
        step = F(1, 1 << 20)
        first = directed_endpoint(*args, F(0), F(1, 32))
        next_value = directed_endpoint(*args, step, F(1, 32))
        slope = ((arb(next_value['exact_coefficient_log2_upper'])
                  -arb(first['exact_coefficient_log2_upper']))*(1 << 20))
        self.assertTrue(slope > expected)
        self.assertTrue(slope-expected < arb(1)/10)


if __name__ == '__main__':
    unittest.main()
