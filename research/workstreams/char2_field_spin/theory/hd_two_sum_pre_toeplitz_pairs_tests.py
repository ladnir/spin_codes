"""Exact tiny law checks and the pre-H majorant interface."""
from fractions import Fraction as F
from math import comb
import unittest

from hd_grouped_toeplitz import binary_field, toy_counts
from hd_two_sum_toeplitz_pairs import pair_kernel, column_contraction, xor_matching
import hd_two_sum_pre_toeplitz_pairs as target


class PreToeplitzTests(unittest.TestCase):
    def test_uniform_input_pair_law_is_flat_on_each_support(self):
        q = 4
        counts, accepted = toy_counts(q, 2)
        for active in range(4):
            for marker in (F(1), F(3, 2)):
                masses = target.pre_pair_support_masses(active, q=q, input_marker=marker, mode='exact-support')
                for output in range(q*q):
                    out_support = int(output % q != 0)+2*int(output // q != 0)
                    actual = F(0)
                    for incoming in range(q*q):
                        support = int(incoming % q != 0)+2*int(incoming // q != 0)
                        if support & ~active == 0:
                            actual += marker**support.bit_count()*F(counts[incoming][output], accepted)
                    self.assertEqual(actual, masses[out_support]/(q-1)**out_support.bit_count())

    def test_majorant_can_follow_arbitrary_positive_H_observable(self):
        q, mul = 4, binary_field(4)
        counts, accepted = toy_counts(q, 2)
        majorant = pair_kernel(q, 'iid-majorant')
        # Fixed invertible H, nonuniform positive coordinate markers.
        values = []
        for output in range(q*q):
            x, y = output % q, output // q
            hx, hy = x ^ y, x ^ mul[2][y]
            values.append(F(2, 3)**bool(hx)*F(3, 5)**bool(hy))
        for incoming in range(q*q):
            support = int(incoming % q != 0)+2*int(incoming // q != 0)
            actual = sum(F(counts[incoming][output], accepted)*values[output] for output in range(q*q))
            upper = sum(majorant[support][int(output % q != 0)+2*int(output // q != 0)]
                        / (q-1)**(int(output % q != 0)+int(output // q != 0))*values[output]
                        for output in range(q*q))
            self.assertLessEqual(actual, upper)

    def test_exact_pair_mass_normalization_and_domination(self):
        for active in range(4):
            for marker in (F(1), F(3, 2), F(300)):
                exact = target.pre_pair_support_masses(active, input_marker=marker, mode='exact-support')
                safe = target.pre_pair_support_masses(active, input_marker=marker)
                self.assertEqual(sum(exact), (1+255*marker)**active.bit_count())
                self.assertTrue(all(x <= y for x, y in zip(exact, safe)))

    def test_uniform_space_mixture_matches_NZ_expansion(self):
        group, z = 8, F(1, 2)
        for active in (0, 1, 3, 15, 255):
            matching = xor_matching(4, group)
            mixture = target.uniform_space_mixture(active, matching, group=group)
            self.assertEqual(sum(mixture), 1)
            direct = target.pre_column_scalar(active, matching, z, group=group, mode='exact-support')
            alternative = F(0)
            for arity, probability in enumerate(mixture):
                uniform_moment = sum(F(comb(arity, k)*255**k, 256**arity)
                                     *target.h_nz_scalar_moment(k, z, group) for k in range(arity+1))
                alternative += probability*uniform_moment
            self.assertEqual(direct, 256**active.bit_count()*alternative)

    def test_full_arity_and_diagonal_controls(self):
        for group in (8, 16):
            matching, z = xor_matching(1, group), F(1, 2)
            full = (1 << group)-1
            self.assertEqual(target.pre_column_scalar(full, matching, z, group=group, mode='exact-support'),
                             (1+255*z)**group)
            for active in (1, 3, full):
                self.assertEqual(target.pre_column_scalar(active, matching, z, group=group, mode='diagonal'),
                                 column_contraction(active, matching, (z,)*group, mode='diagonal'))

    def test_marked_expansion_and_guards(self):
        coefficients = target.expanded_support_polynomial(15, xor_matching(4, 8), group=8,
                                                        input_marker=F(3, 2), mode='exact-support')
        self.assertEqual(sum(coefficients), (1+255*F(3, 2))**4)
        for active in (-1, True, 4):
            with self.assertRaises(ValueError):
                target.pre_pair_support_masses(active)
        with self.assertRaises(ValueError):
            target.pre_pair_support_masses(1, input_marker=F(1, 2))


if __name__ == '__main__':
    unittest.main()
