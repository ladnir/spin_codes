"""Independent exact arithmetic checks for the GF16 outward engine."""

import copy
from decimal import Decimal
from fractions import Fraction as F
import math
import unittest

from gf16_certificate import (CONFIG, class_upper, interval_local_envelope,
                              upper_polynomial_product, verify_sparse)
from q1_batched_screen import fixed_support_multiplicity
from subfield_certificate import Interval


def rational_local(z):
    q, t, s = 16, 16, 4
    f = z+(1-z)/15
    result = []
    for w in range(17):
        p0 = F(sum((-1)**i*math.comb(w, i)*q**max(0, w-s-i) for i in range(w+1)), 15**w)
        g = F(0)
        for zeros in range(4):
            count = fixed_support_multiplicity(t, s, t-zeros, q)
            for active in range(max(0, zeros-t+w), min(w, zeros)+1):
                g += (F(count*math.comb(w, active)*math.comb(t-w, zeros-active), q**s-1)
                      *z**(t-w-zeros+2*active)*f**(w-active))
        ret = (1-p0)/(q**s-1)
        lo, hi = ret*z**16, ret*z**max(0, 13-w)
        if w == 1:
            lo = hi = (13*z**15+z**16+2*z**14)/(16*(q**s-1))
        elif w >= 2:
            inverse = sum(F(math.comb(t, h), (q**s-1)*15**max(0, h-t+w))
                          *z**(t-h)*(1-z)**h for h in range(t+1))
            hi = min(hi, inverse)
        result.append([[z**w*p0, z**w*(1-p0)], [hi, g-lo]])
    return result


class GF16CertificateTests(unittest.TestCase):
    def test_local_intervals_enclose_exact_rationals(self):
        for z in [F(1, 7), F(17, 19), F(1)]:
            actual = interval_local_envelope(Interval.rational(z.numerator, z.denominator))
            expected = rational_local(z)
            for w in range(17):
                for row in range(2):
                    for col in range(2):
                        self.assertLessEqual(F(actual[w][row][col].lo), expected[w][row][col])
                        self.assertGreaterEqual(F(actual[w][row][col].hi), expected[w][row][col])
                        self.assertGreaterEqual(actual[w][row][col].hi, 0)

    def test_positive_convolution_is_rounded_up(self):
        matrices = [[F(1, 3), F(2, 7), F(4, 5), F(5, 11)],
                    [F(2, 9), F(1, 5), F(1, 13), F(2, 3)]]
        upper = [[Interval.rational(x.numerator, x.denominator).hi for x in a] for a in matrices]
        actual = upper_polynomial_product(upper, upper, 2)
        expected = [[F(0)]*4 for _ in range(3)]
        for i, a in enumerate(matrices):
            for j, b in enumerate(matrices):
                for row in range(2):
                    for col in range(2):
                        expected[i+j][2*row+col] += sum(a[2*row+k]*b[2*k+col] for k in range(2))
        for a, b in zip(actual, expected):
            for upper_value, exact in zip(a, b):
                self.assertGreaterEqual(F(upper_value), exact)

    @classmethod
    def setUpClass(cls):
        entry = class_upper(1, 4290114385, 2**32)
        cls.receipt = dict(format='gf16-gl4-sparse-v1', complete=False, config=CONFIG,
                           covered_active_rows=[1, 1], uncovered_active_rows=[2, 8192],
                           entries=[entry], failure_upper=entry['failure_upper'], integer_bits=entry['integer_bits'])

    def test_sparse_class_replay(self):
        self.assertTrue(verify_sparse(self.receipt))
        self.assertGreaterEqual(self.receipt['integer_bits'], 48)

    def test_reject_overstated_receipts(self):
        for field, value in [('complete', True), ('covered_active_rows', [1, 8192]),
                             ('failure_upper', '1e-100'), ('integer_bits', 80)]:
            record = copy.deepcopy(self.receipt)
            record[field] = value
            with self.assertRaises(ValueError):
                verify_sparse(record)
        record = copy.deepcopy(self.receipt)
        record['entries'][0]['failure_upper'] = '1e-100'
        with self.assertRaises(ValueError):
            verify_sparse(record)
        record = copy.deepcopy(self.receipt)
        record['entries'][0]['active_rows'] = 2
        with self.assertRaises(ValueError):
            verify_sparse(record)


if __name__ == '__main__':
    unittest.main()
