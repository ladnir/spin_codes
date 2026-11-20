"""Exact small Johnson-kernel normalization and projector checks."""
import itertools
import math
import unittest
from fractions import Fraction as F
from johnson_overlap import hahn


class JohnsonTests(unittest.TestCase):
    def test_orthogonality_and_projectors(self):
        for n, w in [(6, 2), (7, 3)]:
            subsets = [set(s) for s in itertools.combinations(range(n), w)]
            size = len(subsets)
            for degree in range(w+1):
                values = [hahn(n, w, degree, i) for i in range(w+1)]
                self.assertEqual(values[0], 1)
                mass = sum(math.comb(w, i)*math.comb(n-w, i)*values[i] for i in range(w+1))
                self.assertEqual(mass, size if degree == 0 else 0)
                multiplicity = math.comb(n, degree)-(math.comb(n, degree-1) if degree else 0)
                matrix = [[values[w-len(s & t)] for t in subsets] for s in subsets]
                for i, j in itertools.product(range(size), repeat=2):
                    square = sum(matrix[i][k]*matrix[k][j] for k in range(size))
                    self.assertEqual(square, F(size, multiplicity)*matrix[i][j])
                # Symmetric scaled idempotence proves PSD, exactly.
                selected = [i for i in range(size) if i % 3 != 0]
                self.assertGreaterEqual(sum(matrix[i][j] for i in selected for j in selected), 0)


if __name__ == '__main__':
    unittest.main()
