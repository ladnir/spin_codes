"""Independent exact checks for the prospective two-sum column mixer."""
from collections import Counter
from itertools import combinations
import random
import unittest

from gf8_packet_maps import MUL, matvec, matmul, xor_sum
from hd_two_sum_outer import (LABELS, matrix, direct, circuit, uniform_mask_counts,
                              nonzero_mask_counts, silent_moment_count)


class TwoSumOuterTests(unittest.TestCase):
    def test_moments_and_involution(self):
        self.assertEqual(xor_sum(LABELS), 0)
        self.assertEqual(xor_sum(MUL[i][i] for i in LABELS), 0)
        identity = tuple(tuple(int(i == j) for j in LABELS) for i in LABELS)
        self.assertEqual(matmul(matrix(), matrix()), identity)

    def test_explicit_circuit_every_binary_basis_and_random_bytes(self):
        examples = []
        for coordinate in LABELS:
            for bit in LABELS:
                x = [0]*8
                x[coordinate] = 1 << bit
                examples.append(x)
        rng = random.Random(1132026)
        examples += [[rng.randrange(256) for _ in LABELS] for _ in range(100)]
        for x in examples:
            self.assertEqual(circuit(x), matvec(matrix(), x))
            self.assertEqual(direct(x), circuit(x))
            self.assertEqual(circuit(circuit(x)), tuple(x))

    def test_all_one_and_two_support_projective_vectors(self):
        for i in LABELS:
            actual = Counter()
            x = [0]*8
            x[i] = 1
            actual[sum(bool(v) << j for j, v in enumerate(direct(x)))] = 255
            self.assertEqual(tuple(actual[t] for t in range(256)), nonzero_mask_counts(1 << i))
        for i, j in combinations(LABELS, 2):
            actual = Counter()
            for ratio in range(1, 256):
                x = [0]*8
                x[i], x[j] = 1, ratio
                actual[sum(bool(v) << k for k, v in enumerate(direct(x)))] += 255
            self.assertEqual(tuple(actual[t] for t in range(256)), nonzero_mask_counts((1 << i) | (1 << j)))
            self.assertEqual(sum(actual.values()), 255**2)

    def test_explicit_branch_six_witness_and_silent_moments(self):
        # For three distinct labels i,j,k, (j+k,i+k,i+j) is nonzero
        # and annihilates both moments. H fixes this weight-three vector.
        for i, j, k in combinations(LABELS, 3):
            x = [0]*8
            x[i], x[j], x[k] = j ^ k, i ^ k, i ^ j
            self.assertEqual(direct(x), tuple(x))
            self.assertEqual(sum(bool(v) for v in x), 3)
        self.assertEqual([silent_moment_count(a) for a in range(4)], [1, 0, 0, 255])

    def test_uniform_kernel_boundary_and_invalid_inputs(self):
        self.assertEqual(uniform_mask_counts(0), (1,)+(0,)*255)
        self.assertEqual(uniform_mask_counts(255), tuple(255**t.bit_count() for t in range(256)))
        for function in (uniform_mask_counts, nonzero_mask_counts):
            with self.assertRaises(ValueError):
                function(256)
        with self.assertRaises(ValueError):
            circuit([0]*7)

    def test_complete_kernel_involution_symmetry_and_branch(self):
        kernel = [nonzero_mask_counts(a) for a in range(256)]
        branch = 17
        for active, counts in enumerate(kernel):
            self.assertEqual(sum(counts), 255**active.bit_count())
            for target, count in enumerate(counts):
                self.assertEqual(count, kernel[target][active])
                if active and count:
                    branch = min(branch, active.bit_count()+target.bit_count())
            if active.bit_count() <= 6:
                self.assertEqual(counts[active], silent_moment_count(active.bit_count()))
        self.assertEqual(branch, 6)

    def test_pair_law_has_only_two_fixed_rowset_types(self):
        for i, j in combinations(LABELS, 2):
            active = (1 << i) | (1 << j)
            counts = nonzero_mask_counts(active)
            missing = range(1, 8) if i == 0 else range(8)
            expected = [0]*256
            for row in missing:
                expected[255 ^ (1 << row)] = 255
            expected[255] = 255*(255-len(missing))
            self.assertEqual(tuple(expected), counts)


if __name__ == '__main__':
    unittest.main()
