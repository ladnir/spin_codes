"""Small exact checks for new geometry and count-folding glue, not a certificate."""
from contextlib import redirect_stdout
from fractions import Fraction as Q
from io import StringIO
import unittest
from unittest.mock import patch

from flint import arb, arb_mat
import packet_q1 as screen


class PacketQ1Tests(unittest.TestCase):
    def test_geometry(self):
        bch = screen.Geometry(128, 256, 512)
        rs = screen.Geometry(512, 64, 128)
        self.assertEqual((bch.K, bch.N, bch.macros_per_region), (65536, 131072, 4))
        self.assertEqual((rs.K, rs.N, rs.macros_per_region), (65536, 131072, 16))
        for arguments in ((127, 256, 512), (512, 64, 512), (True, 64, 128), (0, 64, 128)):
            with self.assertRaises(ValueError):
                screen.Geometry(*arguments)

    def test_cdf_needs_suffix_majorant(self):
        # With nonmonotone weights, summing CDF differences as shells is invalid.
        upper = screen.fold_cdf((Q(0), Q(2), Q(3)), (arb(0), arb(1) / 4, arb(1) / 2))
        self.assertEqual(upper, arb(3) / 2)

    def test_geometry_and_total_are_enforced_before_kernel(self):
        common = dict(group_count=32, regions=2, epochs_per_region=1,
                      group_dimension=4, count_kind='shells', tilts=['.01'])
        with self.assertRaises(ValueError):
            screen.evaluate_q1((0, 1, 1), **common)
        with self.assertRaises(ValueError):
            screen.evaluate_q1((0, 5, 10), **dict(common, epochs_per_region=4))

    def test_identity_kernel_counts_all_q1_messages(self):
        # A kernel which never attenuates gives exactly L*(2^groupK-1).
        local = [arb_mat([[1]]), arb_mat([[1]])]
        common = dict(group_count=32, regions=2, epochs_per_region=1,
                      group_dimension=4, tilts=['.01'], data={}, map_record={})
        with patch.object(screen.kernel_t64, 'local_operators', return_value=local), redirect_stdout(StringIO()):
            shells = screen.evaluate_q1((0, 5, 10), count_kind='shells', **common)
            cdf = screen.evaluate_q1((0, 5, 15), count_kind='cdf', **common)
        mantissa, exponent = shells['q1_upper']
        self.assertEqual(Q(mantissa) * Q(2) ** exponent, 32 * 15)
        self.assertEqual(shells['q1_upper'], cdf['q1_upper'])
        self.assertFalse(shells['whole_code_certificate'])
        self.assertEqual(shells['occupancies_not_covered'], [2, 32])


if __name__ == '__main__':
    unittest.main()
