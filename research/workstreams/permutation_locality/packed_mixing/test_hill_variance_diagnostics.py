from fractions import Fraction as Q
import math
import unittest

import hill_variance_diagnostics as diagnostics


class VarianceDiagnosticTests(unittest.TestCase):
    def test_shared_offset_sum_and_removal(self):
        parts = [((Q(0), Q(1, 8)), (0, 0, 2)),
                 ((Q(1, 8), Q(1, 4)), (0, 0, -4))]
        result = diagnostics.summarize(parts, [math.log(2), math.log(3)],
                                       [math.log(5), math.log(7)], math.log(4), 8)
        self.assertAlmostEqual(result['log2_proposal'], math.log2(4*(2*5+3*7)))
        self.assertAlmostEqual(sum(r['fraction_of_proposal'] for r in result['parts']), 1)
        self.assertEqual(result['dominant_part'], 1)
        self.assertEqual(result['descending_parts'], [1, 0])
        self.assertAlmostEqual(result['parts'][0]['without_part_log2'], math.log2(84))
        self.assertAlmostEqual(result['parts'][1]['removal_gain_bits'], math.log2(31/10))
        self.assertAlmostEqual(result['parts'][0]['outer_endpoint_span_bits'], 2/math.log(2))
        self.assertAlmostEqual(result['parts'][1]['outer_endpoint_span_bits'], 4/math.log(2))
        self.assertTrue(result['diagnostic_only'])
        self.assertFalse(result['certificate'])

    def test_single_part_has_no_remaining_counterfactual(self):
        result = diagnostics.summarize([((0, Q(1, 4)), (0, 0, 0))],
                                       [-10000], [2], 3, 2048)
        self.assertAlmostEqual(result['log2_proposal'], -9995/math.log(2))
        self.assertIsNone(result['parts'][0]['without_part_log2'])
        self.assertIsNone(result['parts'][0]['removal_gain_bits'])
        self.assertEqual(result['parts'][0]['fraction_of_proposal'], 1)

    def test_stable_sum_matches_regional_formula(self):
        parts = [((Q(i, 12), Q(i+1, 12)), (0, 0, i)) for i in range(3)]
        outer, inner, common = [10000, 10001, 9999], [-2, -4, -6], 8
        result = diagnostics.summarize(parts, outer, inner, common, 2048)
        terms = [a+b for a, b in zip(outer, inner)]
        peak = max(terms)
        expected = (common+peak+math.log(sum(math.exp(x-peak) for x in terms)))/math.log(2)
        self.assertAlmostEqual(result['log2_proposal'], expected)

    def test_invalid_inputs_fail(self):
        parts = [((0, Q(1, 4)), (0, 0, 0))]
        for outer, inner, common, groups in [([], [0], 0, 1), ([math.inf], [0], 0, 1),
                                            ([0], [0], 0, 0)]:
            with self.assertRaises(ValueError):
                diagnostics.summarize(parts, outer, inner, common, groups)

    def test_zero_contribution(self):
        parts = [((0, Q(1, 8)), (0, 0, 0)),
                 ((Q(1, 8), Q(1, 4)), (0, 0, 0))]
        result = diagnostics.summarize(parts, [-math.inf, 2], [0, 3], 4, 1)
        self.assertEqual(result['parts'][0]['fraction_of_proposal'], 0)
        self.assertEqual(result['parts'][0]['removal_gain_bits'], 0)
        self.assertIsNone(result['parts'][1]['without_part_log2'])
        self.assertEqual(diagnostics._logsumexp([-math.inf]), -math.inf)


if __name__ == '__main__':
    unittest.main()
