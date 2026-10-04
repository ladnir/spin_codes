import importlib.util
from fractions import Fraction as F
from math import comb
from pathlib import Path
import sys
import unittest

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('rs8_transfer', HERE / 'rs8_transfer.py')
transfer = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = transfer
spec.loader.exec_module(transfer)
outer = transfer.helpers.load_module('rs8_test_outer', HERE.parent / 'k16_design/rs_outer.py')


class Rs8TransferTests(unittest.TestCase):
    def test_alternate_mds_formula(self):
        q, n, k = 1 << 32, 8, 4
        d = n-k+1
        expected = [1] + [0] * n
        for w in range(d, n+1):
            expected[w] = comb(n, w)*(q-1)*sum(
                (-1)**j * comb(w-1, j) * q**(w-d-j) for j in range(w-d+1))
        self.assertEqual(list(outer.mds_symbol_weight_counts(q, n, k)), expected)
        self.assertEqual(sum(expected), 1 << 128)

    def test_exact_shells_and_minimum(self):
        counts = outer.expected_group_support_counts(n=8, k=4, packet_bits=4, packets_per_symbol=8)
        self.assertEqual([i for i, x in enumerate(counts) if x], list(range(5, 65)))
        self.assertEqual(sum(counts), (1 << 128)-1)
        self.assertEqual(counts[5], F(comb(8, 5) * 120**5, ((1 << 32)-1)**4))

    def test_every_fixed_support_has_uniform_envelope(self):
        q = 1 << 32
        counts = outer.mds_symbol_weight_counts(q, 8, 4)
        for h in range(5, 9):
            self.assertLessEqual(F(counts[h], comb(8, h)*(q-1)**h), F(1, (q-1)**4))
        self.assertLess(F(1 << 256, (q-1)**4), F(1 << 256, 65535**8))

    def test_scalar_maps_are_transitive(self):
        # Exhaustive small analogue, including the adjoint orientation.
        def mul(a, b):
            result = 0
            while b:
                if b & 1:
                    result ^= a
                a <<= 1
                if a & 16:
                    a ^= 0x13
                b >>= 1
            return result
        for x in range(1, 16):
            self.assertEqual({mul(c, x) for c in range(1, 16)}, set(range(1, 16)))
            self.assertEqual({sum(((mul(c, 1 << bit) & x).bit_count() & 1) << bit
                                  for bit in range(4)) for c in range(1, 16)}, set(range(1, 16)))


if __name__ == '__main__':
    unittest.main()
