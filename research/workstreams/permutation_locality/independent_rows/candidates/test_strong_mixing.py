import unittest

from flint import arb, arb_mat, ctx

import strong_mixing as candidate
from occupancy_memory import Z, M, C, U


class StrongMixingTests(unittest.TestCase):
    def test_retarget_precedes_lazy_mass_replacement(self):
        ctx.prec = 192
        spectrum = {48:3, 56:7, 64:11, 72:13, 80:17}
        base = [arb_mat(11, 11) for _ in range(33)]
        for i,v in enumerate(spectrum,U):
            base[0][i,i] = (-arb('.072')*v).exp()/4
        for j in range(1, 33):
            base[j][M, Z] = arb(3)/16  # Pure refresh before blending.
            base[j][C, Z] = arb(1)/8
            base[j][C, C] = arb(1)/16
            for k, size in enumerate(spectrum.values()):
                base[j][M, U+k] = arb(3)*size/16
        coefficients = {j: (arb(1)/2, arb(1)/4, arb(1)/8) for j in range(1,17)}
        result = dict(candidate.variants(base, coefficients, spectrum, '.072', 4))
        # Refresh ratio is (15/16)/(3/4)=5/4; lazy ratio is 1/4.
        self.assertEqual(result['both-6-16'][6][M, Z], arb(15)/64 + arb(1)/8)
        self.assertEqual(result['both-6-16'][6][M, C], arb(1)/8)
        self.assertEqual(result['both-6-16'][6][C, Z], 0)
        self.assertEqual(result['both-6-16'][6][C, C], 0)
        self.assertEqual(result['both-6-16'][5][M, C], 0)
        self.assertEqual(base[6][M, Z], arb(3)/16)

    def test_retarget_rejects_visible_coupled_density_column(self):
        base = arb_mat(11,11)
        base[M,C] = 1
        with self.assertRaisesRegex(ValueError, 'unsplit'):
            candidate.retarget([base], {}, '.072', 2, 4)


if __name__ == '__main__':
    unittest.main()
