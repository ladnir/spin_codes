import unittest
from itertools import product
from math import comb

from extension_moments import krawtchouk, support_cap


class ExtensionMomentTests(unittest.TestCase):
    def test_polynomials_and_orthogonality(self):
        for q in (2, 4, 16):
            for n in range(1, 7):
                rows = [krawtchouk(n, q, x, n) for x in range(n+1)]
                for x, row in enumerate(rows):
                    for j, value in enumerate(row):
                        exact = sum((-1)**a*(q-1)**(j-a)*comb(x, a)*comb(n-x, j-a)
                                    for a in range(max(0, j-(n-x)), min(j, x)+1))
                        self.assertEqual(value, exact)
                for j in range(n+1):
                    for k in range(n+1):
                        moment = sum(comb(n, x)*(q-1)**x*rows[x][j]*rows[x][k] for x in range(n+1))
                        expected = q**n*comb(n, j)*(q-1)**j if j == k else 0
                        self.assertEqual(moment, expected)

    def test_extended_even_parity_codes(self):
        # Addition in GF(2^a) is XOR. The binary single-parity-check code
        # extends to this q-ary code; its dual repetition code has distance n.
        for q, n in ((2, 7), (4, 5), (16, 3)):
            shells = [0]*(n+1)
            for prefix in product(range(q), repeat=n-1):
                parity = 0
                for value in prefix:
                    parity ^= value
                weight = sum(bool(x) for x in prefix)+bool(parity)
                shells[weight] += 1
            for u in range(2, n+1):
                self.assertLessEqual(sum(shells[1:u+1]), support_cap(n, q, q**(n-1), n, 2, u))

    def test_degree_guard(self):
        with self.assertRaises(ValueError):
            support_cap(7, 2, 64, 7, 2, 4, degree=4)


if __name__ == '__main__':
    unittest.main()
