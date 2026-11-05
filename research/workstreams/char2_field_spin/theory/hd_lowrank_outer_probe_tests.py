"""Small exact controls for the prospective low-rank mixer probe."""
from itertools import combinations
import unittest

from gf8_packet_maps import MUL, matmul, rank
import hd_lowrank_outer_probe as target


class LowrankProbeTests(unittest.TestCase):
    def test_rank2_is_exact_reference(self):
        self.assertEqual(target.matrix(2), target.reference_matrix())

    def test_gram_zero_and_involution(self):
        identity = tuple(tuple(int(i == j) for j in range(16)) for i in range(16))
        for degree in (2, 3, 4):
            v, h = target.vandermonde(degree), target.matrix(degree)
            self.assertFalse(any(x for row in matmul(tuple(zip(*v)), v) for x in row))
            self.assertEqual(matmul(h, h), identity)
            self.assertEqual(rank(h), 16)

    def test_small_gram_nullity_matches_principal_matrix(self):
        for degree in (2, 3, 4):
            v, h = target.vandermonde(degree), target.matrix(degree)
            for arity in (1, 2, 3, 16):
                for active in combinations(range(16), arity):
                    gram = [[int(i == j) for j in range(degree)] for i in range(degree)]
                    for index in active:
                        for i in range(degree):
                            for j in range(degree):
                                gram[i][j] ^= MUL[v[index][i]][v[index][j]]
                    principal = [[h[i][j] for j in active] for i in active]
                    self.assertEqual(arity-rank(principal), degree-rank(gram))

    def test_rank3_affine_gf4_principal_counterexample(self):
        # These labels form 06 + GF4 in the existing AES-field basis.
        # Full H remains invertible, but its four-coordinate restriction
        # cannot use the frozen rank2 principal-inverse argument.
        labels = (0x06, 0x07, 0xba, 0xbb)
        indices = tuple(target.LABELS.index(u) for u in labels)
        v, h = target.vandermonde(3), target.matrix(3)
        gram = [[int(i == j) for j in range(3)] for i in range(3)]
        for index in indices:
            for i in range(3):
                for j in range(3):
                    gram[i][j] ^= MUL[v[index][i]][v[index][j]]
        self.assertEqual(gram, [[1, 0, 0], [0, 1, 1], [0, 1, 1]])
        self.assertEqual(rank([[h[i][j] for j in indices] for i in indices]), 3)


if __name__ == '__main__':
    unittest.main()
