from fractions import Fraction as F
from itertools import product
from math import comb
from unittest import TestCase, main

import hd_shifted_common_source_bad_count as lemma


class BadCountTests(TestCase):
    def test_labelled_independent_flags(self):
        probabilities = (F(1, 3), F(2, 7), F(1, 3), F(2, 7))
        exact = [F(0)]*5
        for flags in product((0, 1), repeat=4):
            mass = F(1)
            for flag, probability in zip(flags, probabilities):
                mass *= probability if flag else 1-probability
            exact[sum(flags)] += mass
        self.assertEqual(lemma.count_distribution(probabilities), tuple(exact))
        self.assertEqual(sum(exact), 1)

    def test_factorial_moment_bounds_each_tail(self):
        p, q = F(1, 5), F(2, 9)
        distribution = lemma.count_distribution((p, q)*3)
        for threshold in range(1, 7):
            moment = sum(comb(k, threshold)*mass for k, mass in enumerate(distribution) if k >= threshold)
            calculated = lemma.two_parity_factorial_bound(p, q, 3, threshold)
            self.assertEqual(calculated, moment)
            self.assertGreaterEqual(calculated, sum(distribution[threshold:]))
            self.assertLessEqual(calculated, comb(6, threshold)*max(p, q)**threshold)

    def test_total_cut_is_not_a_per_column_cut(self):
        # Keeping <=1 in EACH two-region column incorrectly retains [1,1].
        one_column = lemma.count_distribution((F(1, 3), F(1, 4)))
        global_distribution = lemma.count_distribution((F(1, 3), F(1, 4))*2)
        self.assertGreater(sum(one_column[:2])**2, sum(global_distribution[:2]))

    def test_exact_production_lemma_and_scope(self):
        receipt = lemma.exact_lemma()
        bound = F(receipt['exact_union_bound'])
        self.assertLess(bound, F(1, 2**63))
        self.assertGreater(bound, F(1, 2**64))
        self.assertEqual(receipt['retained_bad_regions_at_most'], 12)
        self.assertEqual(receipt['source_full_pair_max_exactly'], 1)
        self.assertFalse(receipt['distance_bound'])
        self.assertFalse(receipt['whole_code_certificate'])
        self.assertFalse(receipt['rowset_or_message_count_charged'])

    def test_threshold_twelve_is_not_forty_five_bits(self):
        probability = max(lemma.load_inputs())
        bound = lemma.union_bound(lemma.SOURCE_COUNT, probability, 256, 12)
        self.assertLess(bound, F(1, 2**43))
        self.assertGreater(bound, F(1, 2**44))

    def test_invalid_probability_or_threshold(self):
        for value in (F(-1, 2), F(3, 2)):
            with self.assertRaises(ValueError):
                lemma.count_distribution((value,))
        with self.assertRaises(ValueError):
            lemma.union_bound(1, F(1, 2), 256, 0)
        with self.assertRaises(ValueError):
            lemma.two_parity_factorial_bound(F(1, 2), F(1, 3), 2, 5)


if __name__ == '__main__':
    main()
