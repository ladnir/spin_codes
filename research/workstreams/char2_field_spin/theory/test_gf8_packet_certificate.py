from fractions import Fraction as F
import unittest

from flint import ctx, fmpq

from gf8_packet_certificate import combine
from gf8_packet_middle_certificate import metadata
from gf8_packet_sparse_certificate import polynomial_power, certify_zone, ANCHORS
from gf8_packet_width_screen import local_matrix


class PacketCertificateTests(unittest.TestCase):
    def receipts(self):
        base = dict(metadata(), complete=True, rigorous=True)
        sparse = dict(base, active_range=[1, 1057], margin_bits_lower=65,
                      failure_upper_dyadic_hex=[hex(3), hex(1 << 67)])
        middle = dict(base, weights=[48661, 389283], margin_bits_lower=3303)
        dense = dict(base, weights=[389284, 1622016], margin_bits_lower=212483)
        return [sparse, middle, dense]

    def test_complete_join_and_union(self):
        result = combine(*self.receipts())
        self.assertEqual(result['margin_bits_lower'], 65)
        self.assertEqual(result['minimum_distance'], 162202)
        self.assertEqual(result['first_remaining_outer_weight'], 48668)
        self.assertEqual(result['requested_rate'], '64/99')

    def test_geometry_and_completion_rejected(self):
        for index in range(3):
            for key, value in (('width', 16), ('eta', 3), ('injection_constant', 1)):
                data = self.receipts()
                data[index][key] = value
                with self.assertRaises(ValueError):
                    combine(*data)
            data = self.receipts()
            data[index]['rigorous'] = False
            with self.assertRaises(ValueError):
                combine(*data)

    def test_join_gaps_rejected(self):
        for index, key, value in ((0, 'active_range', [1, 1056]),
                                  (0, 'active_range', [2, 1057]),
                                  (2, 'weights', [389285, 1622016]),
                                  (2, 'weights', [389284, 1622015])):
            data = self.receipts()
            data[index][key] = value
            with self.assertRaises(ValueError):
                combine(*data)

    def test_bad_probability_rejected(self):
        data = self.receipts()
        data[0]['failure_upper_dyadic_hex'] = ['0x1', '0x2']
        with self.assertRaises(ValueError):
            combine(*data)

    def test_single_packet_directed_polynomial(self):
        local = [local_matrix(w, F(9, 10), F(7), 8, 4) for w in range(9)]
        actual = polynomial_power(local, 1, 8)
        from math import comb
        for index, (row, col) in enumerate(((0, 0), (0, 1), (1, 0), (1, 1))):
            for w in range(9):
                expected = local[w][row][col]*comb(8, w)
                # Compare exact rational endpoints: contains(fmpq) first
                # rounds the rational to a ball, whose radius may be wider.
                exact = fmpq(expected.numerator, expected.denominator)
                self.assertLessEqual(actual[index][w].lower().fmpq(), exact)
                self.assertLessEqual(exact, actual[index][w].upper().fmpq())

    def test_directed_one_row_and_precision_restoration(self):
        precision = ctx.prec
        receipt = certify_zone(1, 1, 1, ANCHORS[1], precision=512)
        self.assertGreaterEqual(receipt['certified_integer_bits'], 230)
        self.assertEqual(ctx.prec, precision)
        p, q = (int(v, 16) for v in receipt['failure_upper_dyadic_hex'])
        self.assertLessEqual(F(p, q), F(1, 1 << receipt['certified_integer_bits']))


if __name__ == '__main__':
    unittest.main()
