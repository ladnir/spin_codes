import importlib.util
from fractions import Fraction as F
import json
from math import comb
from pathlib import Path
import sys
import unittest

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('k18_rs8_certificate', HERE / 'k18_rs8_certificate.py')
certificate = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = certificate
spec.loader.exec_module(certificate)
outer = certificate.pinned_module('k18_rs8_test_outer', certificate.OUTER, certificate.OUTER_SHA)


class K18Rs8Tests(unittest.TestCase):
    def test_alternate_mds_counts_and_normalization(self):
        q, n, k = 1 << 32, 8, 4
        d = n - k + 1
        expected = [1] + [0] * n
        for w in range(d, n + 1):
            expected[w] = comb(n, w) * (q - 1) * sum(
                (-1)**j * comb(w - 1, j) * q**(w - d - j) for j in range(w - d + 1))
        self.assertEqual(list(outer.mds_symbol_weight_counts(q, n, k)), expected)
        counts = outer.expected_group_support_counts(8, 4, 4, 8)
        self.assertEqual(sum(counts), (1 << 128) - 1)
        self.assertEqual([v for v, count in enumerate(counts) if count], list(range(5, 65)))
        self.assertEqual(counts[5], F(comb(8, 5) * 120**5, (q - 1)**4))

    def test_uniform_envelope_each_column_support(self):
        q = 1 << 32
        counts = outer.mds_symbol_weight_counts(q, 8, 4)
        for h in range(5, 9):
            self.assertLessEqual(F(counts[h], comb(8, h) * (q - 1)**h), F(1, (q - 1)**4))
        self.assertLess(F(1 << 256, (q - 1)**4), F(1 << 256, 65535**8))

    def test_triangular_fold_matches_full_ordered_sum(self):
        counts = [F(v * v + 1, v + 1) for v in range(65)]
        bounds = [[F((v + u) % 13 + 1, 14) for u in range(v + 1)] for v in range(65)]
        expected = comb(17, 2) * sum(
            counts[v] * counts[u] * bounds[max(v, u)][min(v, u)]
            for v in range(65) for u in range(65))
        self.assertEqual(certificate.fold_two(counts, bounds, 17), expected)

    def test_one_and_two_group_factors(self):
        counts = [F()] * 65
        counts[5], counts[9] = F(3), F(7)
        one = [F(1, 11)] * 65
        two = [[F(1, 13)] * (v + 1) for v in range(65)]
        self.assertEqual(certificate.fold_one(counts, one, 7), F(7 * 10, 11))
        self.assertEqual(certificate.fold_two(counts, two, 7), F(comb(7, 2) * 100, 13))
        self.assertEqual(certificate.fold_two(counts, two, 7, low_support_only=True),
                         F(comb(7, 2) * (3 * 3 + 2 * 3 * 7), 13))

    def test_full_support_interface_rejects_omissions_and_bad_probabilities(self):
        one = dict(support_probability_uppers=[[1, 0] for _ in range(65)])
        self.assertEqual(len(certificate.probabilities(one, 1)), 65)
        one['support_probability_uppers'][5] = [2, 0]
        with self.assertRaises(ValueError):
            certificate.probabilities(one, 1)
        two = dict(support_pair_probability_uppers=[[[1, 0]] * (v + 1) for v in range(65)])
        self.assertEqual(sum(map(len, certificate.probabilities(two, 2))), 2145)
        two['support_pair_probability_uppers'][8].pop()
        with self.assertRaises(ValueError):
            certificate.probabilities(two, 2)

    def test_missing_occupancy_rejected(self):
        counts = outer.expected_group_support_counts(8, 4, 4, 8)
        selected = {str(q):[1, -100] for q in range(1, 2049)}
        del selected['107']
        with self.assertRaisesRegex(ValueError, 'every q1..2048'):
            certificate.calculate(counts, [F(1)] * 65,
                                  [[F(1)] * (v + 1) for v in range(65)], selected)

    def test_outward_dyadics_and_exact_serialization(self):
        for value in (F(1, 3), F(1, 7 << 4096), F(1 << 4096, 31)):
            endpoint = certificate.helpers.ceil_dyadic(value)
            self.assertGreaterEqual(certificate.helpers.dyadic(endpoint), value)
        tiny_sum = F(1, 1 << 30000) + F(1, 1 << 58)
        encoded = certificate.encode_exact(tiny_sum)
        self.assertEqual(certificate.helpers.dyadic([int(encoded['numerator_hex'], 16),
                                                    encoded['exponent']]), tiny_sum)

    def test_complete_authenticated_certificate(self):
        result = certificate.run()
        self.assertTrue(result['whole_code_certificate'])
        self.assertTrue(result['target_58_met'])
        self.assertFalse(result['target_60_met'])
        self.assertEqual(set(result['occupancy_uppers']), {str(q) for q in range(1, 2049)})
        self.assertEqual(result['q2_triangular_pairs_evaluated'], 2145)
        self.assertAlmostEqual(result['component_margin_bits_display']['q1'], 58.245633868649065, places=8)
        self.assertAlmostEqual(result['component_margin_bits_display']['q2'], 126.09293539386795, places=8)
        self.assertAlmostEqual(result['component_margin_bits_display']['tail'], 75.61227961565965, places=8)
        reference = json.loads(certificate.REFERENCE.read_text(encoding='utf-8'))
        self.assertTrue(all(result['occupancy_uppers'][str(q)] ==
                            reference['selected_occupancy_bounds'][str(q)] for q in range(3, 2049)))
        low = result['low_support_accounting']
        self.assertFalse(low['global_exclusion_used'])
        for q in ('q1', 'q2'):
            self.assertGreater(certificate.helpers.dyadic(low[f'{q}_charge_upper']), 0)
            self.assertLess(certificate.helpers.dyadic(low[f'{q}_charge_upper']),
                            certificate.helpers.dyadic(result['component_uppers'][q]))
        exact = result['exact_union']
        total = certificate.helpers.dyadic([int(exact['numerator_hex'], 16), exact['exponent']])
        self.assertEqual(total, sum(map(certificate.helpers.dyadic, result['occupancy_uppers'].values()), F()))
        self.assertLessEqual(total, certificate.helpers.dyadic(result['union_upper']))


if __name__ == '__main__':
    unittest.main()
