"""Exact/directed regression tests for the repaired eight-symbol candidate."""
from fractions import Fraction as F
import math
import unittest

import numpy as np

from gf8_packet_middle_certificate import (CONSTRUCTION, RETURN_MINIMA, WITNESSES,
                                          cell_upper, certify, certify_dense,
                                          exact_transfer, metadata, prefactor)
from gf8_packet_width_screen import packet, source_row_log, bound
from return_cover_certificate import Witness, interval


class RepairedPacketCertificateTests(unittest.TestCase):
    def test_fresh_geometry_and_map_metadata(self):
        g = metadata()
        self.assertEqual(g['construction'], CONSTRUCTION)
        self.assertEqual((g['outer'], g['width'], g['state']), ([128, 83], 8, 4))
        self.assertEqual((g['rows'], g['n'], g['padded_k'], g['output_threshold']),
                         (12672, 1622016, 1051776, 162201))
        self.assertEqual(g['requested_rate'], '64/99')
        self.assertEqual(g['return_minima'], list(RETURN_MINIMA))

    def test_exact_matrix_agrees_with_floating_packet(self):
        for _, w in (WITNESSES[0], WITNESSES[4], WITNESSES[-1]):
            expected = np.asarray(exact_transfer(w), dtype=float)
            actual = np.exp(packet(math.log(w.x), math.log(w.y), math.log(w.tau),
                                   8, 4, RETURN_MINIMA))
            np.testing.assert_allclose(actual, expected, rtol=2e-12, atol=1e-300)

    def test_single_cell_matches_relaxed_source_float(self):
        g = metadata()
        _, w = WITNESSES[0]
        weight = 48661
        direct = cell_upper(weight, weight, g, w, prefactor(w, g))
        point = [math.log(getattr(w, name)) for name in ('x', 'y', 'tau', 'tilt')]
        screened = bound(point, weight, dimension=83, width=8, state=4,
                         minima=RETURN_MINIMA)
        # Directed replay conservatively retains the old relaxed source.
        alpha = weight/g['n']
        correction = g['rows']*(source_row_log(alpha, point[2], point[3], 128, 83, False)
                                -source_row_log(alpha, point[2], point[3], 128, 83, True))/math.log(2)
        self.assertAlmostEqual(float(direct.hi.ln())/math.log(2), screened+correction, places=7)

    def test_interval_dominates_exact_weight_partition(self):
        g = metadata()
        _, w = WITNESSES[4]
        lo, hi = 162201, 162207
        upper = cell_upper(lo, hi, g, w, interval(1))
        individual = sum((cell_upper(j, j, g, w, interval(1)) for j in range(lo, hi+1)), interval(0))
        self.assertGreaterEqual(upper.lo, individual.hi)

    def test_complete_middle_and_dense_components(self):
        middle, dense = certify(), certify_dense()
        for receipt in (middle, dense):
            self.assertTrue(receipt['complete'] and receipt['rigorous'])
            self.assertFalse(receipt['whole_code'])
            for key, value in metadata().items():
                self.assertEqual(receipt[key], value)
        self.assertEqual(middle['weights'], [48661, 389283])
        self.assertEqual(dense['weights'], [389284, 1622016])
        self.assertEqual(middle['weights'][1]+1, dense['weights'][0])
        self.assertEqual(middle['margin_bits_lower'], 3303)
        self.assertEqual(dense['margin_bits_lower'], 212483)
        self.assertEqual(middle['partition_sha256'],
                         'a6bb9147e34372abfd664a6c75d5ea2bcbdec8c6fb3dbcf7bd99d37828216c9f')
        self.assertEqual(dense['witness_sha256'],
                         'efcdb9e1ded74c57ee496378e5665d96820755cf5b3aa2079aab01e5a2628b5b')

    def test_reject_inexact_or_invalid_witness(self):
        with self.assertRaises(ValueError):
            exact_transfer(Witness(.5, F(1, 2), F(7), F(1)))
        with self.assertRaises(ValueError):
            exact_transfer(Witness(F(1), F(1), F(256), F(1)))
        with self.assertRaises(ValueError):
            certify(witnesses=((.03, WITNESSES[0][1]),))
        with self.assertRaises(ValueError):
            certify(cell_bits=39)
        with self.assertRaises(ValueError):
            certify(first=F(1, 4), last=F(1, 5))


if __name__ == '__main__':
    unittest.main()
