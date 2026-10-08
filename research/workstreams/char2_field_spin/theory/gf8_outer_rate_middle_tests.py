"""Geometry, directed-arithmetic, and coverage regressions for rate sweep."""
from fractions import Fraction as F
import math
import unittest

from gf8_outer_rate_middle import (WITNESSES, cell_upper, certify, conditioning,
                                   metadata, profile)
from return_cover_certificate import interval, terminal_power
from state8_density_certificate import transfer
from symmetric_density_screen import bound


class OuterRateMiddleTests(unittest.TestCase):
    def test_natural_padding_and_effective_rate(self):
        for dimension in range(87, 94):
            g = profile(16, dimension)
            self.assertEqual(g['rows'], 768)
            self.assertEqual(g['requested_rate'], '2/3')
            self.assertGreaterEqual(g['padded_k'], 2**16)
            self.assertLess(g['padded_k']-2**16, 64*dimension)
        self.assertEqual(profile(16, 94)['rows'], 704)
        g = profile(20, 88)
        self.assertEqual((g['rows'], g['n'], g['padded_k'], g['output_threshold']),
                         (11968, 1531904, 1053184, 153190))
        self.assertEqual(g['requested_rate'], '128/187')

    def test_modal_factor_bounds_every_column_mean(self):
        for rows in range(2, 20, 2):
            reciprocal = conditioning(rows)
            self.assertLessEqual(reciprocal, conditioning(rows, False))
            for hits in range(rows+1):
                probability = math.comb(rows, hits)*F(hits, rows)**hits*F(rows-hits, rows)**(rows-hits)
                self.assertGreaterEqual(probability*reciprocal, 1)

    def test_point_matches_floating_control_without_conditioning_change(self):
        g = profile(16, 88)
        _, witness = WITNESSES[4]
        weight = round(g['n']/10)
        prefactor = (interval(g['rows']+1)**128*terminal_power(transfer(witness), g['n']//16)
                     /(interval(witness.tau)*interval(witness.y)**g['threshold']))
        directed = cell_upper(weight, weight, g, witness, prefactor)
        point = [math.log(getattr(witness, name)) for name in ('x', 'y', 'tau', 'tilt')]
        floating = bound(point, weight, 16, 128, 88, 8, joint=True)
        self.assertAlmostEqual(float(directed.hi.ln())/math.log(2), floating, places=7)

    def test_interval_dominates_sum_of_single_weights(self):
        g = profile(16, 88)
        _, witness = WITNESSES[4]
        lo, hi = g['n']//10, g['n']//10+8
        covering = cell_upper(lo, hi, g, witness, interval(1))
        individual = sum((cell_upper(w, w, g, witness, interval(1))
                          for w in range(lo, hi+1)), interval(0))
        self.assertGreaterEqual(covering.lo, individual.hi)

    def test_complete_k20_dimension88_replay(self):
        receipt = certify()
        self.assertEqual(receipt['weights'], [45958, 367656])
        self.assertEqual(receipt['covered_weights'], 321699)
        self.assertGreaterEqual(receipt['class_margin_bits_lower'], 40)
        self.assertTrue(receipt['complete'] and receipt['rigorous'])
        self.assertFalse(receipt['whole_code'])
        self.assertEqual(receipt['partition_sha256'],
                         '1a91d49d27785379e794ae9161b1432afc5c213ab60728af5c1ac64c0bb169b5')
        for key, value in metadata(20, 88).items():
            self.assertEqual(receipt[key], value)

    def test_complete_k18_dimension88_replay(self):
        receipt = certify(18, 88, first=F(1, 25))
        self.assertEqual(receipt['weights'], [15401, 92405])
        self.assertEqual(receipt['covered_weights'], 77005)
        self.assertGreaterEqual(receipt['class_margin_bits_lower'], 129)
        self.assertEqual(receipt['partition_sha256'],
                         '2dbf084acd15c1e04cabb0a387e382e4667936e758a65f4e26d939746584a138')
        self.assertEqual(receipt['requested_rate'], '32/47')

    def test_reject_bad_scope_or_float_witnesses(self):
        for args in ((15, 88), (20, 128), (20, 0)):
            with self.assertRaises(ValueError):
                profile(*args)
        with self.assertRaises(ValueError):
            certify(witnesses=((.03, WITNESSES[0][1]),))
        with self.assertRaises(ValueError):
            certify(first=F(1, 4), last=F(1, 5))
        with self.assertRaises(ValueError):
            certify(cell_bits=39)


if __name__ == '__main__':
    unittest.main()
