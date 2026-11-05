import unittest
import numpy as np

from hd_inverse_column_control import schedule
from hd_nonquiet_word_control import first_nonquiet_words


class FirstNonquietTests(unittest.TestCase):
    def test_independent_difference_for_exact_matrices(self):
        # At exponent one, exact full=q+n gives a direct independent oracle.
        rows = tuple((i,) for i in range(16))
        quiet = (np.arange(64, dtype=float).reshape(16, 2, 2)+1)/128
        nonquiet = quiet[::-1]/2
        full = quiet+nonquiet
        for maps in ((tuple(range(16)),)*16, schedule()):
            result = first_nonquiet_words(np.log2(quiet), np.log2(nonquiet),
                                         np.log2(full), rows, maps)
            expected = []
            for i in range(16):
                q, total = np.array([1., 0.]), np.array([1., 0.])
                for p in maps:
                    q = q @ quiet[p[i]]
                    total = total @ full[p[i]]
                expected.append(np.log2(total.sum()-q.sum()))
            np.testing.assert_allclose(result, expected, rtol=1e-13, atol=1e-12)

    def test_zero_nonquiet_branch(self):
        zero = np.full((1, 2, 2), -np.inf)
        quiet = np.zeros((1, 2, 2))
        rows = (tuple(range(16)),)
        result = first_nonquiet_words(quiet, zero, quiet, rows, schedule())
        self.assertEqual(result[0], -np.inf)


if __name__ == '__main__':
    unittest.main()
