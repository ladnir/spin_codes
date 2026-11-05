from fractions import Fraction as F
from math import log2
import unittest

from flint import arb, ctx

from hd_grouped_toeplitz import group_cap
from hd_shifted_outer_density import (DEFAULT_OUTER_MEASURE, TOEPLITZ_OUTER_MEASURE,
    bound_description, checked_outer_measure, density_exact, log_density_arb,
    log_density_float)


class OuterDensityTests(unittest.TestCase):
    def test_matches_established_column_atom_and_information_set(self):
        sharp = density_exact(TOEPLITZ_OUTER_MEASURE)
        self.assertEqual(sharp, 256**256*group_cap())
        inherited = density_exact()
        self.assertEqual(inherited/sharp, F(256, 255)**96)
        self.assertGreater(sharp, 1)
        self.assertLess(sharp, inherited)

    def test_default_float_expression_is_unchanged(self):
        self.assertEqual(log_density_float(), 2048-128*log2(255))
        self.assertAlmostEqual(log_density_float()-log_density_float(TOEPLITZ_OUTER_MEASURE),
                               96*log2(256/255), places=11)

    def test_directed_encloses_exact_density_log_at_both_precisions(self):
        for precision in (256, 384):
            with ctx.workprec(precision):
                for mode in (DEFAULT_OUTER_MEASURE, TOEPLITZ_OUTER_MEASURE):
                    density = density_exact(mode)
                    direct = (arb(density.numerator)/arb(density.denominator)).log()/arb(2).log()
                    self.assertTrue(log_density_arb(mode).overlaps(direct))
                improvement = log_density_arb()-log_density_arb(TOEPLITZ_OUTER_MEASURE)
                self.assertTrue(improvement > arb(542)/1000)
                self.assertTrue(improvement < arb(543)/1000)

    def test_explicit_supported_modes_only(self):
        for value in ('extension-field4', 'dense-random', '', None, 1024):
            with self.assertRaises(ValueError):
                checked_outer_measure(value)
        self.assertFalse(bound_description(TOEPLITZ_OUTER_MEASURE)['setup_changed'])


if __name__ == '__main__':
    unittest.main()
