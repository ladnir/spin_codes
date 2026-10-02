"""Independent exact tiny-input checks for capped newborn measures."""
from fractions import Fraction as Q
from math import comb
import unittest

import numpy as np
from flint import arb, ctx
import kernel_maps as kernel
import kernel_birth_density as density
import sparse_kernel
import scalar_cover as sc


def endpoint(value):
    m, e = value.upper().man_exp()
    return Q(int(m))*Q(2)**int(e)


class BirthDensityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.images = tuple(sum(((s >> j)&1)*v for j, v in enumerate((1, 6, 120))) for s in range(8))
        cls.columns = [1, 2, 4, 3, 5, 7, 6, 1]
        cls.data = kernel.prepare_maps(cls.images, cls.columns, 3)

    def setUp(self):
        ctx.prec = 192

    def exact(self, z):
        result = [[Q(0)]*8 for _ in range(3)]
        for word in range(256):
            j = bool(word & 15)+bool(word >> 4)
            syndrome = 0
            for bit, column in enumerate(self.columns):
                if word >> bit & 1:
                    syndrome ^= column
            result[j][syndrome] += z**word.bit_count()/(comb(2, j)*15**j)
        return result

    def test_exact_walsh_laws_and_all_candidate_decompositions(self):
        for z in (Q(1), Q(3, 4), Q(1, 64)):
            za = kernel.aq(z)
            numerators, denominators = density.conditional_laws(self.data, za)
            exact = self.exact(z)
            local = sparse_kernel.outward_at_z(self.data, za)
            for j in range(3):
                for state in range(8):
                    upper = Q(int(numerators[state, j]), denominators[j])
                    self.assertGreaterEqual(upper, exact[j][state])
                    self.assertLess(upper-exact[j][state], Q(1, 2)**38)
                baseline = [local[j][0, k] for k in range(local[j].ncols())]
                for row in density.outward_rows(self.data, baseline, numerators[:, j], denominators[j]):
                    self.assertGreaterEqual(endpoint(row[0]), exact[j][0])
                    for target, level in enumerate(self.data['birth_class_levels'], 3):
                        residual = sum(max(Q(0), exact[j][s]-endpoint(row[2])/7)
                            for s in range(1, 8) if self.images[s].bit_count() == level)
                        self.assertGreaterEqual(endpoint(row[target]), residual)

    def test_iid_refinement_and_float_proposal(self):
        for p in (Q(1, 5), Q(3, 4)):
            probabilities = sc.probabilities(p)
            matrix = kernel.outward(self.data, probabilities, Q(1, 4))
            floating = kernel.floating(self.data, np.array(list(map(float, probabilities))), .25)
            n = matrix.nrows()
            np.testing.assert_allclose(floating,
                [[float(matrix[i, j]) for j in range(n)] for i in range(n)], rtol=1e-10, atol=2e-12)
            z = (-kernel.aq(Q(1, 4))).exp()
            local = sparse_kernel.outward_at_z(self.data, z)
            refined = density.refine_local(self.data, local, z, p)
            self.assertEqual(len(refined), len(local))
            for before, after in zip(local, refined):
                for i in range(1, n):
                    for j in range(n):
                        self.assertEqual(before[i, j], after[i, j])

    def test_proof_option_does_not_change_map_identity(self):
        baseline = kernel.prepare_maps(self.images, self.columns, 3, birth_density='classes')
        self.assertEqual(baseline['map_sha256'], self.data['map_sha256'])
        self.assertNotIn('birth_character_indices', baseline)


if __name__ == '__main__':
    unittest.main()
