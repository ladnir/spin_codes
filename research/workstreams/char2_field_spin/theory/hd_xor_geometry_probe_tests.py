"""Independent multiplicities and repeated-pair controls for FLOAT probes."""
import unittest

import numpy as np

from hd_xor_column_potential import PAIRS, case_plans, column_matrices
from hd_xor_multigroup import region_operators
from hd_xor_geometry_probe import (assignment_multiplicities,
    holder_geometry_screen, repeated_pair_matrices, repeated_pair_diagnostic)


class GeometryProbeTests(unittest.TestCase):
    def test_multiset_multiplicities_count_every_ordered_assignment(self):
        for single in range(4):
            for double in range(4-single):
                counts, denominator = assignment_multiplicities(single, double)
                self.assertEqual(sum(counts), 8**single*28**double)
                self.assertEqual(denominator, 8**single*28**double)
                self.assertTrue(all(count >= 1 for count in counts))

    def test_repeated_pair_matrices_match_ordered_two_group_checker(self):
        z = .99875
        exact = column_matrices(z, 2)[0, 2]
        controls = repeated_pair_matrices(region_operators(z, 2), 2)
        plans = case_plans(2)[0, 2]
        for index, pair in enumerate(PAIRS):
            position = next(i for i, (first, second, law) in enumerate(plans) if first == second == pair)
            np.testing.assert_allclose(controls[index], exact[position], rtol=1e-13, atol=1e-12)

    def test_holder_moment_does_not_exceed_universal_maximum(self):
        result = holder_geometry_screen(3, .998589552, -62.40625)
        for category, value in result['category_holder_caps'].items():
            self.assertLessEqual(value, result['category_maxima'][category]+1e-12)
        self.assertLess(result['holder_category_log2_bound'], result['maximum_category_log2_bound'])
        self.assertFalse(result['directed'])

    def test_repeated_pair_control_exposes_dominant_r17_category(self):
        result = repeated_pair_diagnostic()
        self.assertAlmostEqual(max(result['repeated_pair_row_caps']),
                               result['resource_envelope_pure_double_cap'], places=10)
        self.assertEqual(result['dominant_cancellation_counts'], [17, 0])
        self.assertFalse(result['directed'])


if __name__ == '__main__':
    unittest.main()
