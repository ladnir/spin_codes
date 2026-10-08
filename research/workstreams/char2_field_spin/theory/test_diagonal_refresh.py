"""Exact algebra and directed-class regressions for the separate new inner."""

from fractions import Fraction as F
from math import comb
import unittest

from diagonal_refresh_rank import exact_checks, rank_potential_bound
from diagonal_refresh_certificate import Witness, certify, transfer


class DiagonalRefreshTests(unittest.TestCase):
    def test_exhaustive_rank_potential(self):
        self.assertEqual(exact_checks(), {'exact_potential_cases': 882,
                         'exhausted_diagonal_subspace_pairs': 50421})

    def test_integer_tilts_stay_exact(self):
        for r in range(1, 5):
            for z in range(5):
                value = rank_potential_bound(r, z, 255)
                self.assertIsInstance(value, F)
                self.assertTrue(0 <= value <= 1)

    def test_transfer_independent_intersection_count(self):
        witness = Witness(F(2), F(1, 4), F(5, 2), F(1))
        x, y, tau = witness.x, witness.y, witness.tau
        reference = [[F(1), F(0)], [F(0), F(0)]]
        for w in range(1, 17):
            term = comb(16, w)*(x*y)**w
            reference[0][1] += term
            if w >= 5:
                reference[0][0] += term/tau**4
        # Independently index by input/output weights and their intersection.
        for w in range(17):
            for v in range(17):
                for overlap in range(max(0, w+v-16), min(w, v)+1):
                    roots = 16-w-v+overlap
                    if roots >= 4:
                        continue
                    gamma = max(rank_potential_bound(r, roots, tau) for r in range(1, 5))
                    count = comb(16, w)*comb(w, overlap)*comb(16-w, v-overlap)
                    term = count*x**w*y**v/tau**(w-overlap)*gamma
                    reference[1][1] += term
                    if w and v and overlap and (v != 1 or w == 16):
                        reference[1][0] += term/tau**min(overlap, 4)
        actual = transfer(witness)
        for i in range(2):
            for j in range(2):
                self.assertLessEqual(F(actual[i][j].lo), reference[i][j])
                self.assertGreaterEqual(F(actual[i][j].hi), reference[i][j])

    def test_directed_k16_complete_class(self):
        result = certify(16)
        self.assertTrue(result['rigorous'])
        self.assertFalse(result['whole_code'])
        self.assertEqual(result['minimum_input_weight'], 10589)
        self.assertEqual(result['maximum_input_weight'], 23101)
        self.assertEqual(result['checked_cells'], 1300)
        self.assertEqual(result['class_margin_bits_lower'], 202)


if __name__ == '__main__':
    unittest.main()
