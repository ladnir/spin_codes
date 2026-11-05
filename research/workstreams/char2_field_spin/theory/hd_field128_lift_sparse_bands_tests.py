"""Tiny replay and unfinished-coverage tests; no extension-field distance gate."""

from fractions import Fraction as F
import unittest

from flint import arb, ctx

import hd_field128_lift_certificate as cert
import hd_field128_lift_sparse_bands as bands
import hd_field128_lift_screen as floating


class SparseBandTests(unittest.TestCase):
    def test_default_schedule_and_exact_witnesses(self):
        schedule = bands.scheduled_bands()
        self.assertEqual(len(schedule), 9)
        self.assertEqual([r for _, lo, hi in schedule for r in range(lo, hi+1)], list(range(1, 65)))
        self.assertEqual(bands.scheduled_bands(3, 12), ((4, 3, 5), (8, 6, 10), (16, 11, 12)))
        for center, alpha in ((1, F(1, 16)), (8, F(1, 16)), (16, F(1, 40)),
                              (24, F(1, 80)), (32, F(1, 80)), (48, F(1, 160)), (64, F(1, 160))):
            z, actual = bands.default_witness(center)
            self.assertEqual(z, (1-F(center, 1200 if center <= 17 else 1400))**16)
            self.assertEqual(actual, alpha)
        for schedule in (((1, 1, 1), (2, 3, 3)), ((1, 1, 2), (2, 2, 3))):
            with self.assertRaises(ValueError):
                bands.scheduled_bands(1, 3, schedule)

    def test_shared_directed_table_matches_three_single_weights(self):
        z, alpha = F(7, 8), F(1, 2)
        result = bands.batch_receipt(1, 3, z, alpha, maximum=5, q=256, max_seconds=20)
        self.assertEqual(result['regional_tables_built'], 1)
        self.assertEqual(result['coefficient_maximum'], 5)
        self.assertEqual(result['coefficient_marker'], '5/4096')
        self.assertEqual([row['active_groups'] for row in result['weight_receipts']], [1, 2, 3])
        for row in result['weight_receipts']:
            old = cert.weight_receipt(row['active_groups'], z, alpha, q=256, max_seconds=20)
            self.assertLessEqual(abs(F(row['log2_upper_rational'])-F(old['log2_upper_rational'])), F(1, 10**10))
            self.assertEqual(row['profiles'], old['profiles'])
            self.assertEqual(row['witness']['coefficient_maximum'], 5)
        with ctx.workprec(256):
            mass = sum((arb(2)**cert.rational(F(row['log2_upper_rational']))
                        for row in result['weight_receipts']), arb(0))
            expected = cert.rational_ceiling(mass.log()/arb(2).log())
            self.assertEqual(F(result['log2_upper_rational']), expected)
        self.assertTrue(result['covers_all_supports_in_range'])
        self.assertFalse(result['whole_code_certificate'])

    def test_shared_float_table_matches_selected_screens(self):
        z, alpha = F(7, 8), F(1, 2)
        result = bands.float_band_screen(1, 3, z, alpha, maximum=5, q=256)
        for row in result['rows']:
            old = floating.weight_screen(row['active_groups'], float(z), float(alpha), q=256)
            self.assertAlmostEqual(row['log2_upper_screen'], old['log2_upper_screen'], places=8)
        self.assertEqual(result['regional_tables_built'], 1)
        self.assertFalse(result['directed'])

    def test_serial_driver_retains_complete_and_pending_ranges(self):
        snapshots = []
        result = bands.run_schedule(first=1, last=3, q=256, max_seconds=20,
                                    witnesses={c: (F(7, 8), F(1)) for c in (1, 2, 4)},
                                    checkpoint=snapshots.append)
        self.assertEqual(result['completed_band_count'], 3)
        self.assertTrue(result['complete_requested_range'])
        self.assertEqual(result['completed_weight_range'], [1, 3])
        self.assertEqual(result['pending_weight_ranges'], [])
        self.assertTrue(result['source_pins_unchanged_at_exit'])
        self.assertFalse(result['whole_code_certificate'])
        self.assertEqual(snapshots[0]['pending_weight_ranges'], [[1, 3]])
        self.assertEqual(snapshots[1]['pending_weight_ranges'], [[2, 3]])
        self.assertEqual(snapshots[2]['pending_weight_ranges'], [[3, 3]])
        incomplete = bands.run_schedule(first=1, last=3, q=256, max_seconds=1e-12)
        self.assertFalse(incomplete['complete_requested_range'])
        self.assertEqual(incomplete['completed_bands'], [])
        self.assertEqual(incomplete['pending_weight_ranges'], [[1, 3]])
        self.assertEqual(incomplete['failure']['error_type'], 'TimeoutError')

    def test_source_closure_and_bounds(self):
        pins = bands.source_pins()
        for name in ('hd_field128_lift_sparse_bands.py', 'hd_field128_lift_certificate.py',
                     'hd_common_refresh_principal.py', 'hd_field128_lift_screen.py'):
            self.assertIn(name, pins)
        for first, last, maximum in ((0, 3, 3), (2, 1, 3), (1, 4, 3), (1, 65, 65)):
            with self.assertRaises(ValueError):
                bands.batch_receipt(first, last, F(7, 8), F(1), maximum=maximum)
        with self.assertRaises(ValueError):
            bands.batch_receipt(1, 1, .9, F(1))


if __name__ == '__main__':
    unittest.main()
