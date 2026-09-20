"""Exact one-column conditioning and Fourier parity identities on toy slices."""
from collections import Counter
from fractions import Fraction as F
import itertools
import math
import unittest
from certify_tail_parity_mixing import kraw_values


class ConditionalParityTests(unittest.TestCase):
    def test_exact_conditioning(self):
        for n, w, q in [(4, 2, 2), (4, 2, 3), (6, 2, 3)]:
            words = [sum(1 << j for j in support)
                     for support in itertools.combinations(range(n), w)]
            counts, good = Counter(), Counter()
            for rows in itertools.product(words, repeat=q):
                j = sum(row & 1 for row in rows)
                parity = 0
                for row in rows:
                    parity ^= row
                counts[j] += 1
                good[j] += parity == 0
            one, zero = kraw_values(n-1, w-1), kraw_values(n-1, w)
            for j in range(q+1):
                self.assertEqual(F(counts[j], len(words)**q),
                                 math.comb(q, j)*F(w, n)**j*(1-F(w, n))**(q-j))
                conditional = sum(math.comb(n-1, h)*
                                  F(one[h], math.comb(n-1, h))**j*
                                  F(zero[h], math.comb(n-1, h))**(q-j)
                                  for h in range(n))/2**(n-1)
                if j % 2:
                    self.assertEqual(conditional, 0)
                self.assertEqual(F(good[j], counts[j]), conditional)


if __name__ == '__main__':
    unittest.main()
