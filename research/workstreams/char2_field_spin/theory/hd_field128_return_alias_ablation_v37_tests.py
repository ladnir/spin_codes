"""Tiny coefficient oracles; no bank/files."""
import unittest
from itertools import product
from math import log2
import numpy as np
import hd_field128_return_alias_ablation_v37 as pilot
import hd_field128_gamma_all_pair_screen as pairs


class Tests(unittest.TestCase):
    def test_copy(self):
        raw = np.arange(12.).reshape(3, 2, 2)
        raw.setflags(write=False)
        ablated = pilot.ablate(raw)
        self.assertTrue(np.all(np.isneginf(ablated[:, 1, 0])))
        np.testing.assert_array_equal(raw.reshape(3, 4)[:, [0, 1, 3]], ablated.reshape(3, 4)[:, [0, 1, 3]])
        self.assertTrue(np.isfinite(raw).all())

    def test_coefficients(self):
        matrices = np.array([[[1., 2.], [3., 4.]], [[2., 1.], [4., 3.]], [[3., 4.], [1., 2.]]])
        logs = pilot.ablate(np.log2(matrices))
        matrices[:, 1, 0] = 0
        result = np.exp2(pairs.all_pair_coefficients(logs, .37, slots=1))
        oracle = np.zeros((4, 3, 3))
        types = [(0, 0), (1, 0), (0, 1)]
        for a, b in product(range(3), repeat=2):
            s, f = np.add(types[a], types[b])
            oracle[:, s, f] += ((matrices[a]@matrices[b])**.37).ravel()
        np.testing.assert_allclose(result, oracle, rtol=2e-14)
        self.assertTrue(np.all(result[2] == 0))

    def test_fft(self):
        # Noncommuting coefficients, degree x only.
        a = np.array([[.2, .4], [0., .3]])
        b = np.array([[.5, .1], [0., .2]])
        coefficients = np.stack((a.ravel(), b.ravel()), axis=1)[:, :, None]
        with np.errstate(divide='ignore'):
            logs = np.log2(coefficients)
        x, blocks = 1.25, 4
        baseline = a+x*b+np.array([[0., .2], [.7, 0.]])
        exact = [np.eye(2)]
        for _ in range(blocks):
            out = [np.zeros((2, 2)) for _ in range(len(exact)+1)]
            for degree, term in enumerate(exact):
                out[degree] += term@a
                out[degree+1] += term@b*x
            exact = out
        for grid in (4, 8):
            table = pilot.triangular_atlas(logs, (log2(x), 0.), np.log2(baseline), grid=grid, blocks=blocks)
            wanted = np.zeros((4, grid, grid))
            for degree, term in enumerate(exact):
                wanted[:, degree % grid, 0] += term.ravel()
            # Allow the frozen diagnostic cushion.
            np.testing.assert_allclose(np.exp2(table.logs)[:, :, 0], wanted[:, :, 0], rtol=1e-9, atol=1e-10)
            np.testing.assert_allclose(np.exp2(table.positive_numerator).reshape(2, 2),
                                       np.linalg.matrix_power(a+x*b, blocks), rtol=1e-13)
            self.assertTrue(np.all(np.isneginf(table.logs[2])))
        self.assertFalse(np.allclose(a@b, b@a))

    def test_guard(self):
        with self.assertRaises(ValueError):
            pilot.triangular_atlas(np.zeros((4, 2, 2)), (0., 0.), np.zeros((2, 2)))


if __name__ == '__main__':
    unittest.main()
