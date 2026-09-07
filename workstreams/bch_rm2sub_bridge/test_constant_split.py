from fractions import Fraction as F
import itertools
import math
import unittest

import numpy as np
from flint import arb

import certify_constant_split as split


class ConstantSplitTests(unittest.TestCase):
    def test_all_band_sequences_and_constant_offsets(self):
        q = 5
        region = [tuple(arb((j+1)*(k+1))/32 for k in range(16)) for j in range(q+1)]
        mantissas, exponents = split.density.initial(region)
        ps, roots = [F(1, 4), F(3, 4)], [F(2), F(3, 2)]
        left = np.array([float(r*(1-p)) for r, p in zip(roots, ps)])
        right = np.array([float(r*p) for r, p in zip(roots, ps)])
        rows = list(split.constant_matrices(mantissas, exponents, left, right))
        self.assertEqual([d for d, _, _ in rows], list(range(q+1)))
        for d, matrix, exponent in rows:
            h = q-d
            upper = [F(float(v))*F(2)**exponent for v in matrix.flat]
            for labels in itertools.product(range(2), repeat=d):
                law = [F(1)]
                scale = F(1)
                for label in labels:
                    p = ps[label]
                    law = [sum((law[j-bit]*(p if bit else 1-p)
                        for bit in (0, 1) if 0 <= j-bit < len(law)), F(0))
                        for j in range(len(law)+1)]
                    scale *= roots[label]
                for k in range(16):
                    exact = scale*sum((mass*F((h+j+1)*(k+1), 32) for j, mass in enumerate(law)), F(0))
                    self.assertGreaterEqual(upper[k], exact)

    def test_case_partition_count(self):
        for length in range(2, 9):
            for q in range(1, length+1):
                count = sum(math.comb(length, q)*math.comb(q, h)*12**(q-h) for h in range(q+1))
                self.assertEqual(count, math.comb(length, q)*13**q)


if __name__ == '__main__':
    unittest.main()
