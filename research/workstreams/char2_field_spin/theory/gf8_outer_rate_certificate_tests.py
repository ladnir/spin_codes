"""Aggregator contract tests use synthetic receipts, not proof certificates."""
from copy import deepcopy
from fractions import Fraction as F
import unittest

from gf8_outer_rate_certificate import certify_k18_dimension88, combine, metadata, probability


class RateCertificateContractTests(unittest.TestCase):
    def fixtures(self):
        meta = metadata(20, 88)
        n = meta['n']
        common = dict(**meta, complete=True, rigorous=True)
        # Deliberately synthetic positive bounds test the join, not the code.
        sparse = dict(**common, active_range=[1, 1200], failure_upper_dyadic=['1', str(2**45)])
        middle = dict(**common, weights=[49000, n//4], failure_upper=str(F(1, 2**46)))
        dense = dict(**common, weights=[n//4+1, n], failure_upper_dyadic_exponent=47)
        return sparse, middle, dense

    def test_geometry_and_padding_plateau(self):
        self.assertEqual({metadata(16, d)['n'] for d in range(87, 94)}, {98304})
        self.assertEqual(metadata(16, 94)['n'], 90112)
        self.assertEqual(metadata(20, 88)['rows'], 11968)
        self.assertEqual(metadata(20, 88)['n'], 1531904)
        for k_log in (16, 18, 20):
            for dimension in range(80, 100):
                m = metadata(k_log, dimension)
                self.assertGreaterEqual(m['padded_k'], m['requested_k'])
                self.assertEqual(m['rows'] % 64, 0)
                self.assertLess(F(m['output_threshold'], m['n']), F(1, 10))
                self.assertGreaterEqual(F(m['output_threshold']+1, m['n']), F(1, 10))

    def test_exact_probability_sum(self):
        receipts = self.fixtures()
        result = combine(20, 88, *receipts)
        total = sum(map(probability, receipts), F(0))
        self.assertGreaterEqual(F(result['failure_upper']), total)
        self.assertLessEqual(F(result['failure_upper']), F(1, 2**44))
        self.assertEqual(result['whole_code_margin_bits_lower'], 44)
        self.assertEqual(result['component_margin_bits_lower'], [45, 46, 47])

    def test_metadata_and_rigor_mismatch(self):
        for key in (*metadata(20, 88), 'complete', 'rigorous'):
            receipts = deepcopy(self.fixtures())
            receipts[1][key] = None
            with self.assertRaises(ValueError):
                combine(20, 88, *receipts)

    def test_missing_coverage(self):
        for position, key, replacement in (
                (0, 'active_range', [2, 1200]),
                (0, 'active_range', [1, 1190]),
                (1, 'weights', [49000, 50000]),
                (2, 'weights', [100000, 1531903])):
            receipts = deepcopy(self.fixtures())
            receipts[position][key] = replacement
            with self.assertRaises(ValueError):
                combine(20, 88, *receipts)

    def test_invalid_or_insufficient_probabilities(self):
        for endpoint in ('0', '-1/2', '1', 0.1):
            with self.assertRaises(ValueError):
                probability(dict(failure_upper=endpoint))
        receipts = deepcopy(self.fixtures())
        receipts[0]['failure_upper_dyadic'] = ['1', str(2**40)]
        with self.assertRaises(ArithmeticError):
            combine(20, 88, *receipts)

    def test_real_k18_dimension88_certificate(self):
        receipt = certify_k18_dimension88()
        self.assertEqual(receipt['component_margin_bits_lower'], [72, 129, 15031])
        self.assertEqual(receipt['whole_code_margin_bits_lower'], 72)
        self.assertEqual(receipt['effective_requested_rate'], '32/47')
        self.assertEqual(receipt['minimum_distance'], 38503)
        self.assertEqual(receipt['sparse']['active_range'], [1, 375])
        self.assertEqual(receipt['middle']['weights'], [15401, 92405])
        self.assertEqual(receipt['dense']['weights'], [92406, 385024])
        self.assertTrue(receipt['rigorous'] and receipt['whole_code'])
        self.assertLess(F(receipt['failure_upper']), F(1, 2**72))


if __name__ == '__main__':
    unittest.main()
