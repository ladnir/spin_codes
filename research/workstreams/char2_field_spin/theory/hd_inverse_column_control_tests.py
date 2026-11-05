import unittest
import numpy as np

from hd_inverse_column_control import gf16_multiply, schedule, scheduled_words
from hd_post_refresh_quiet_rows import log_matrix_power


class InverseColumnTests(unittest.TestCase):
    def test_field_and_schedule(self):
        permutations = schedule()
        self.assertEqual(len(set(permutations)), 16)
        for j, p in enumerate(permutations):
            self.assertEqual(sorted(p), list(range(16)))
            for i in range(16):
                self.assertEqual(gf16_multiply(i ^ j, p[i]), int(i != j))

    def test_identity_recovers_old_power(self):
        rows = tuple((i,) for i in range(16))
        matrices = np.log2(np.arange(64, dtype=float).reshape(16, 2, 2)+1.)
        got = scheduled_words(matrices, rows, (tuple(range(16)),)*16)
        old = np.logaddexp2.reduce(log_matrix_power(matrices, 16)[:, 0, :], axis=1)
        np.testing.assert_allclose(got, old, atol=1e-12, rtol=1e-13)

    def test_independent_ordered_products(self):
        rows = tuple((i,) for i in range(16))
        values = (np.arange(64, dtype=float).reshape(16, 2, 2)+1.)/64
        got = scheduled_words(np.log2(values), rows, schedule())
        expected = []
        for i in range(16):
            row = np.array([1., 0.])
            for p in schedule():
                row = row @ values[p[i]]
            expected.append(np.log2(row.sum()))
        np.testing.assert_allclose(got, expected, atol=1e-12, rtol=1e-13)


if __name__ == '__main__':
    unittest.main()
