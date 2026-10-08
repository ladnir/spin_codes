"""Bounded exact tests of the affine source/rank proof objects."""

from fractions import Fraction as F
from itertools import product
import unittest

import hd_affine_proof as proof


class AffineProofTests(unittest.TestCase):
    def test_all_gf4_row_words_and_nonzero_diagonals(self):
        code = {tuple(a ^ proof.mul(b, x) for x in range(4))
                for a, b in product(range(4), repeat=2)}
        self.assertEqual(len(code), 16)
        for word in product(range(4), repeat=4):
            count = 0
            for diagonal in product(range(1, 4), repeat=4):
                scaled = tuple(proof.mul(x, d) for x, d in zip(word, diagonal))
                count += scaled in code
            self.assertEqual(F(count, 81), proof.source_probability(word))

    def test_prescribed_support_enumerator_and_check_bound(self):
        for b, k, q in ((4, 2, 4), (64, 32, 65536)):
            self.assertEqual(sum(proof.prescribed_mds_count(b, k, w, q)*proof.comb(b, w)
                                 for w in range(b+1)), q**k)
            for w in range(1, b+1):
                probability = proof.row_probability(b, k, w, q)
                self.assertGreaterEqual(probability, 0)
                self.assertLessEqual(probability, F(1, (q-1)**(b-k)))

    def test_global_prefix_boundary_is_not_reset(self):
        middle = (0, 1, 1, 1, 0, 0, 0, 0)
        differences = proof.difference_rows(middle, 2, (1, 1, 1, 1))
        # Region1 begins at the preceding region's live carry, hence no jump.
        self.assertEqual(differences, [0, 0, 1, 0, 1, 0, 0, 0])

    def test_full_map_first_moment_and_generic_rank(self):
        matrix = proof.tiny_output_matrix()
        result = proof.exact_preweight_ranks(matrix)
        self.assertEqual(result['histogram'], {2: 459, 3: 6102})
        self.assertEqual(result['failure'], F(17, 243))
        self.assertEqual(result['generic_rank'], 3)
        self.assertEqual(proof.exact_source_first_moment(matrix), result['projective_first_moment'])
        self.assertEqual(proof.block_rank_moment(matrix), F(29, 81))
        self.assertLessEqual(result['projective_first_moment'], proof.block_rank_moment(matrix))

    def test_block_rank_identity_directly(self):
        matrix = proof.tiny_output_matrix()
        z = F(1, 9)
        direct = F(0)
        for word in proof.projective_vectors(3):
            values = proof.matvec(matrix, word)
            active = sum(any(values[start:start+4]) for start in (0, 4))
            direct += z**active
        self.assertEqual(direct, proof.block_rank_moment(matrix))

    def test_positive_two_sided_rank_contraction(self):
        matrix = proof.tiny_output_matrix()
        z, mark = F(1, 9), F(1, 2)
        direct = F(0)
        for word in proof.projective_vectors(3):
            values = proof.matvec(matrix, word)
            active = sum(any(values[start:start+4]) for start in (0, 4))
            direct += mark**sum(x != 0 for x in word)*z**active
        self.assertEqual(direct, proof.block_rank_weight_moment(matrix, mark))
        self.assertEqual(proof.block_rank_moment(matrix),
                         proof.block_rank_weight_moment(matrix, F(1)))

    def test_component_contraction_nullity_including_global_carries(self):
        matrix = proof.tiny_output_matrix(support=tuple(range(8)))
        backward = proof.inverse_matrix(matrix)
        for row_mask in range(4):
            constrained = [r*4+c for r in range(2) if row_mask >> r & 1 for c in range(4)]
            free = [i for i in range(8) if i not in constrained]
            for output_mask in range(256):
                live = [i for i in range(8) if output_mask >> i & 1]
                zero = [i for i in range(8) if not (output_mask >> i & 1)]
                a = [[matrix[i][j] for j in live] for i in constrained]
                h = [[backward[i][j] for j in free] for i in zero]
                self.assertEqual(len(live)-proof.rank(a), len(free)-proof.rank(h))

    def test_same_support_profile_different_value_profile(self):
        repeated = proof.regional_source_average(4, (1, 1, 1, 1), (0, 0, 1, 1))
        distinct = proof.regional_source_average(4, (1, 1, 1, 1), (0, 0, 1, 2))
        self.assertEqual(repeated['probability'], F(5645, 86093442))
        self.assertEqual(distinct['probability'], F(21797, 1377495072))
        self.assertNotEqual(repeated['probability'], distinct['probability'])

    def test_affine_layout_cannot_be_replaced_by_independent_rows(self):
        values = (0, 0, 0, 0, 0, 0, 1, 1)
        common = proof.regional_source_average(8, (1, 1, 1, 1), values)
        affine = proof.regional_source_average(8, (1, 3, 5, 7), values)
        self.assertEqual(common['probability'], F(5627, 108884466432))
        self.assertEqual(affine['probability'], F(113, 1296243648))
        self.assertLess(affine['feasible_arrangements'], common['feasible_arrangements'])
        self.assertGreater(affine['probability'], common['probability'])

    def test_k20_dense_requirement_is_not_a_possible_uniform_lemma(self):
        receipt = proof.dense_profile()
        self.assertEqual(receipt['uniform_active_rows_needed'], 30192)
        self.assertEqual(receipt['impossible_uniform_premise_by_rows'], 19661)
        self.assertGreater(receipt['uniform_active_rows_needed'],
                           receipt['impossible_uniform_premise_by_rows'])


if __name__ == '__main__':
    unittest.main()
