from fractions import Fraction as Q
from itertools import product
import unittest

import joint_count_region as joint
from test_cell_search import source


class RegionTests(unittest.TestCase):
    def test_count_moments_cover_every_small_product_bernoulli_law(self):
        for probabilities in product((Q(0), Q(1, 4), Q(1, 2), Q(1)), repeat=3):
            mean = sum(probabilities)/3
            variance = sum(p*(1-p) for p in probabilities)/3
            if not 0 < mean < 1:
                continue
            mass = [Q(1)]
            for p in probabilities:
                new = [Q(0)]*(len(mass)+1)
                for j, value in enumerate(mass):
                    new[j] += value*(1-p)
                    new[j+1] += value*p
                mass = new
            constraints = joint.count_moments(3, (mean, mean), (variance, variance))
            for features, bound in constraints:
                self.assertLessEqual(sum(p*f for p, f in zip(mass, features)), bound)

    def test_interval_second_moment_accounts_for_mean_offset(self):
        constraints = joint.count_moments(4, (Q(1, 4), Q(3, 4)), (Q(0), Q(1, 4)))
        self.assertEqual(constraints[2][0], (Q(4), Q(1), Q(0), Q(1), Q(4)))
        self.assertEqual(constraints[2][1], 2)
        self.assertEqual(constraints[3][1], 0)

    def test_selection_does_not_silently_accept_duplicates_or_outside_parts(self):
        self.assertEqual(joint.selected_parts(None, 3), {0, 1, 2})
        self.assertEqual(joint.selected_parts([], 3), set())
        for values in ([0, 0], [-1], [3], [True]):
            with self.assertRaises(ValueError):
                joint.selected_parts(values, 3)

    def test_scope_rejects_other_code_or_claim(self):
        scope = source()
        scope.update(schema='packed-canonical-gl32-hill-dense-context-1',
            comparison='direct-expected-shell-majorant')
        joint.validate_scope(scope)
        for field, bad in [('updates', 3), ('minimum_groups', 34), ('K', 1 << 19),
                ('threshold', 209715), ('variance_bins', 65), ('schema', 'arbitrary')]:
            with self.assertRaises(ValueError):
                joint.validate_scope(dict(scope, **{field: bad}))


if __name__ == '__main__':
    unittest.main()
