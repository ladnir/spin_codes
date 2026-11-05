from fractions import Fraction as F
from itertools import product
from math import log2
from unittest import TestCase, main
from unittest.mock import patch

import numpy as np

import hd_two_sum_late_alpha_word_order as gate


class WholeWordOrderTests(TestCase):
    def compare(self, alpha=.5):
        return gate.order_comparison(np.log2([1., 2., 3.]), columns=3, floor=4,
            alpha=alpha, log_source_count=log2(7), log_rowset_count=log2(5),
            log_global_marker=log2(2), log_word_normalizer=log2(64))

    def test_whole_word_tail_matches_exact_words(self):
        raw = [F(1), F(2), F(3)]
        total = sum((raw[a]*raw[b]*raw[c] for a, b, c in product(range(3), repeat=3)
                     if a+b+c >= 4), F(0))
        result = self.compare()
        self.assertAlmostEqual(result['raw_word_tail_log2'], log2(total))
        mean = 5*2*total/64
        self.assertAlmostEqual(result['fixed_source_all_rowsets_mean_log2'], log2(mean))
        self.assertAlmostEqual(result['direct_first_moment_alpha1_log2'], log2(7*mean))

    def test_both_fractional_relaxations_are_in_upper_direction(self):
        result = self.compare()
        self.assertGreater(result['degree_order_gain_bits'], 0)
        self.assertAlmostEqual(result['rowset_order_gain_bits'], .5*log2(5))
        self.assertGreater(result['coefficientwise_alpha_rowset_union_log2'],
                           result['whole_word_alpha_rowset_union_log2'])
        self.assertGreater(result['whole_word_alpha_rowset_union_log2'],
                           result['whole_word_and_rowset_alpha_log2'])

    def test_alpha1_all_orders_agree(self):
        result = self.compare(1.)
        expected = result['direct_first_moment_alpha1_log2']
        for key in ('coefficientwise_alpha_rowset_union_log2',
                    'whole_word_alpha_rowset_union_log2', 'whole_word_and_rowset_alpha_log2'):
            self.assertAlmostEqual(result[key], expected)

    def test_source_count_is_never_fractionalized(self):
        args = dict(columns=2, floor=1, alpha=.2, log_rowset_count=3.,
                    log_global_marker=5., log_word_normalizer=8.)
        first = gate.order_comparison([0., 1.], log_source_count=0., **args)
        second = gate.order_comparison([0., 1.], log_source_count=17., **args)
        for key in ('coefficientwise_alpha_rowset_union_log2', 'whole_word_alpha_rowset_union_log2',
                    'whole_word_and_rowset_alpha_log2', 'direct_first_moment_alpha1_log2'):
            self.assertAlmostEqual(second[key]-first[key], 17.)

    def test_class_sum_inside_fractional_bound_and_component_split(self):
        # Positive components may be split after their complete class sums.
        # min(1,sum x_i)<=sum min(1,x_i)<=sum x_i^alpha_i.
        components = [F(1, 16), F(1, 9), F(3, 2)]
        self.assertLessEqual(min(F(1), sum(components)), sum(min(F(1), x) for x in components))
        self.assertLessEqual(sum(min(F(1), x) for x in components),
                             float(components[0])+float(components[1])**.5+float(components[2])**.2)
        # A source has several deterministic row-set contributions, no new sampling.
        self.assertLess((float(F(1, 16)+F(1, 9)))**.5, float(F(1, 16))**.5+float(F(1, 9))**.5)

    def test_guards_and_scope(self):
        with self.assertRaises(ValueError):
            self.compare(1.1)
        with self.assertRaises(ValueError):
            gate.run(max_seconds=31)
        with self.assertRaises(ValueError):
            gate.checked_artifact({}, 1)
        curve = np.arange(19, dtype=float)*-2
        source = dict(two_atom_envelope=dict(log2_nonzero_atom=-2., log2_probability=0.),
                      moment_control=dict(majorant_term_log2=99.))
        with patch.object(gate, 'read_source', return_value=(source, curve, 7)):
            with patch.object(gate, 'quiet_mixture_coefficients', return_value=np.array([0., 1.])):
                with patch.object(gate, 'order_comparison', return_value=dict(coefficientwise_alpha_rowset_union_log2=99.)):
                    result = gate.run()
        self.assertTrue(result['complete_requested_control'])
        self.assertTrue(result['source_count_outside_alpha'])
        self.assertTrue(result['rowset_sum_inside_alpha'])
        self.assertEqual(result['raw_regional_banks_built'], 0)
        self.assertFalse(result['whole_code_certificate'])


if __name__ == '__main__':
    main()
