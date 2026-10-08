"""Exact regressions and the directed K20 permitted-singleton gate."""

from fractions import Fraction as F
from itertools import combinations
import unittest

import hd_component_rank as proof
from hd_additive_two_state import support_law as four_law
from hd_three_product import support_law as three_law
from hd_affine_proof import projective_vectors


class ComponentRankTests(unittest.TestCase):
    def test_shared_control_dominates_both_support_laws(self):
        for q in (4, 16, 65536):
            for mark in (F(1, 2), F(9, 10), F(1)):
                controls = proof.collapsed_majorant(mark, q)
                for law in (four_law, three_law):
                    for active in (0, 1):
                        for old in range(4):
                            for new_live in (0, 1):
                                actual = sum(p*mark**out for (new, out), p in law(old, active, q).items()
                                             if int(new != 0) == new_live)
                                self.assertLessEqual(actual, controls[active][int(old != 0)][new_live])

    def test_run_sum_equals_every_small_polynomial_coefficient(self):
        for rows in range(1, 9):
            for q, mark in ((4, F(1, 2)), (16, F(9, 10))):
                kernels = proof.exact_kernels(rows, q, mark, collapsed=True)
                for h, kernel in enumerate(kernels):
                    self.assertEqual(proof.run_kernel_exact(rows, h, q, mark), kernel)

    def test_directed_run_upper_against_exact_small_kernels(self):
        from flint import arb, ctx
        previous = ctx.prec
        ctx.prec = 128
        try:
            for rows in (2, 5, 8):
                for h in range(rows+1):
                    exact = proof.run_kernel_exact(rows, h, 16, F(9, 10))
                    upper = proof.run_kernel_upper(rows, h, 16, F(9, 10))
                    for i in range(2):
                        for j in range(2):
                            self.assertTrue(upper[i][j]+arb(2)**-100 > proof.arb_fraction(exact[i][j]))
        finally:
            ctx.prec = previous

    def test_component_transfer_equals_all_component_flag_assignments(self):
        for cuts in (((1,), (2,), (0,)), ((0, 2), (1, 3)), ((0, 1), (0, 1))):
            kernels = proof.exact_kernels(4, 4, F(1, 2))
            exact = proof.component_moment(cuts, 4, 4, kernels.__getitem__)
            direct = proof.direct_component_moment(cuts, 4, 4, kernels)
            self.assertEqual(exact, direct)

    def test_mark_one_counts_the_entire_component_space(self):
        kernels = proof.exact_kernels(4, 4, F(1))
        for cuts in (((1,), (2,), (0,)), ((0, 2), (1, 3))):
            moment = proof.component_moment(cuts, 4, 4, kernels.__getitem__)
            self.assertEqual(moment, 4**sum(map(len, cuts))-1)

    def test_positive_permitted_row_decomposition(self):
        rows, regions, q, mark = 2, 2, 4, F(1, 2)
        z = F(1, 3)
        kernels = proof.exact_kernels(rows, q, mark)
        direct = F(0)
        for difference in projective_vectors(rows*regions):
            active = sum(any(difference[r*regions:(r+1)*regions]) for r in range(rows))
            carry = 0
            vector = [F(1), F(0), F(0), F(0)]
            for c in range(regions):
                occupancy = 0
                for r in range(rows):
                    carry ^= difference[r*regions+c]
                    occupancy += carry != 0
                matrix = kernels[occupancy]
                vector = [sum(vector[i]*matrix[i][j] for i in range(4)) for j in range(4)]
            direct += z**active*sum(vector)
        expansion = F(0)
        for size in range(1, rows+1):
            for allowed in combinations(range(rows), size):
                moment = proof.component_moment([allowed]*regions, rows, q, kernels.__getitem__)
                expansion += z**size*(1-z)**(rows-size)*moment/(q-1)
        self.assertEqual(expansion, direct)

    def test_uniform_occupancy_bounds_at_all_small_h(self):
        for rows in range(1, 9):
            q, mark = 16, F(9, 10)
            p = F(1, q-1)
            base = mark+(1-mark)*p
            eta = (1+p*p)**rows-1
            birth = mark*(1-base**rows)/(rows*(1-base))
            bound = [[p*p*(1+p*p)**(rows-1), birth+eta],
                     [p*p*(birth+eta), base**rows+eta]]
            for h in range(1, rows+1):
                kernel = proof.run_kernel_exact(rows, h, q, mark)
                for i in range(2):
                    for j in range(2):
                        self.assertLessEqual(kernel[i][j], bound[i][j])
                self.assertEqual(kernel[1][0], p*p*kernel[0][1])

    def test_common_mark_obstruction_is_not_a_class_split_claim(self):
        receipt = proof.common_mark_obstruction()
        self.assertEqual(receipt['small_mark_log2_lower_strict'], 805580)
        self.assertEqual(receipt['large_mark_log2_lower'], 14470333)
        self.assertTrue(receipt['applies_to_common_mark_only'])
        self.assertFalse(receipt['code_counterexample'])

    def test_directed_k20_all_singleton_term(self):
        receipt = proof.singleton_uniform_receipt()
        self.assertTrue(receipt['all_singletons_below_2_to_minus_40'])
        self.assertTrue(receipt['any_fixed_regionwise_row_order'])
        self.assertFalse(receipt['whole_distance_certificate'])
        self.assertIn('-84.9378624997', receipt['all_singletons_log2_upper'])

    def test_selected_pair_profiles_remain_only_selected_profiles(self):
        receipt = proof.pair_receipts()
        self.assertTrue(all(p['below_2_to_minus_40'] for p in receipt['selected_profiles']))
        self.assertFalse(receipt['all_pairs_bound'])
        self.assertFalse(receipt['whole_distance_certificate'])


if __name__ == '__main__':
    unittest.main()
