from copy import deepcopy
from fractions import Fraction as F
from math import comb
import unittest

from hd_bit_syndrome_scheduled_expurgation import (
    CONSTRUCTION, EXPECTED, checked_histogram, deterministic_low_exclusions,
    event_union, log_ceiling,
)


def complete_summary():
    return dict(construction=CONSTRUCTION, field_order=256, columns=16,
        all_requested_arities_complete=True, source_rows_fixed_across_columns=True,
        records=[dict(arity=a, complete=True, labelled_sets_checked=comb(16, a),
                      expected_labelled_sets=comb(16, a),
                      quiet_dimension_sum_histogram=counts)
                 for a, counts in EXPECTED.items()])


class ScheduledExpurgationTests(unittest.TestCase):
    def test_fixed_array_product_dimension_and_projective_factor(self):
        # Three column spaces of dimensions 1,0,2 have 4^(1+0+2)-1
        # nonzero arrays, not 3*(4^3-1) or 4^(3*3)-1.
        per, total = event_union({1: {0: 1}, 2: {3: 2}}, 2,
                                groups=3, q=4, row_dimension=1)
        self.assertEqual(per[1], 0)
        self.assertEqual(total, F(3*2*(4**3-1), 3**3))
        _, excluded = event_union({1: {0: 1}, 2: {3: 2}}, 2,
                                  excluded_arities=(2,), groups=3, q=4, row_dimension=1)
        self.assertEqual(excluded, 0)

    def test_exact_histogram_scope_and_reject_incomplete(self):
        original = complete_summary()
        histogram = checked_histogram(original)
        self.assertEqual(histogram[8], {48: 12390, 49: 480})
        for a in range(1, 10):
            self.assertEqual(sum(histogram[a].values()), comb(16, a))
        bad = deepcopy(original)
        bad['records'][0]['quiet_dimension_sum_histogram'][0] += 1
        with self.assertRaises(ValueError):
            checked_histogram(bad)
        with self.assertRaises(ValueError):
            checked_histogram({**original, 'source_rows_fixed_across_columns': False})

    def test_exact_union_and_deterministic_variant(self):
        histogram = checked_histogram(complete_summary())
        per, total = event_union(histogram, 9)
        self.assertEqual(per[9], F(8192*11440*(256**64-1), 255**73))
        self.assertEqual(per[8], F(8192*(12390*(256**48-1)+480*(256**49-1)), 255**65))
        _, reduced = event_union(histogram, 9, excluded_arities=(1, 2, 3, 4, 5))
        self.assertEqual(total-reduced, per[4]+per[5])
        self.assertLess(reduced, total)
        self.assertLess(total, F(1, 2**40))

    def test_directed_ceiling_and_guards(self):
        for precision in (256, 384):
            ceiling = log_ceiling(F(1, 2**40), precision)
            self.assertGreaterEqual(ceiling, -40)
            self.assertLessEqual(ceiling, -40 + F(1, 10**12))
            upper = log_ceiling(F(1, 3), precision)
            self.assertGreater(upper, F(-1584962500722, 10**12))
            self.assertLess(upper, F(-1584962500719, 10**12))
        with self.assertRaises(ValueError):
            event_union({1: {-1: 1}}, 1)
        with self.assertRaises(ValueError):
            event_union({1: {0: 1}}, 2)

    def test_deterministic_exclusion_requires_separate_incidence_evidence(self):
        incidence = dict(
            four_row_summary=dict(construction=CONSTRUCTION,
                                  fixed_source_sets_checked=1820, maximum_quiet_columns=4),
            activity_incidence=dict(construction=CONSTRUCTION,
                all_requested_arities_complete=True,
                results=[dict(arity=5, complete=True, labelled_sets_checked=4368,
                              excluded_all_quiet_sets=4368, not_excluded_sets=0,
                              universal_nonquiet_column_floor=2)]))
        self.assertEqual(deterministic_low_exclusions(incidence), (1, 2, 3, 4, 5))
        incidence['activity_incidence']['results'][0]['not_excluded_sets'] = 1
        with self.assertRaises(ValueError):
            deterministic_low_exclusions(incidence)


if __name__ == '__main__':
    unittest.main()
