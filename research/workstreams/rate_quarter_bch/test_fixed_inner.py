"""Rate, cutoff, map-identity, and receipt checks for the quarter-rate screen."""
import json
import math
import unittest
from fractions import Fraction

import numpy as np
import evaluate_fixed_inner as screen
import verify_fixed_inner_results as replay


class FixedInnerTests(unittest.TestCase):
    def test_exact_dense_partition_check(self):
        replay.check_coverage([dict(lower=[0,0,0],upper=[5,7,7])], 7, 2)
        with self.assertRaisesRegex(ValueError, 'gap|incomplete'):
            replay.check_coverage([dict(lower=[0,0,0],upper=[4,7,7])], 7, 2)
        with self.assertRaisesRegex(ValueError, 'overlap'):
            replay.check_coverage([dict(lower=[0,0,0],upper=[5,7,7])]*2, 7, 2)

    def test_quarter_rate_geometry(self):
        g = screen.geometry(20, Fraction(19, 100))
        self.assertEqual(g['outer_rows'], 16384)
        self.assertEqual(g['output_bits'], 4 * (1 << 20))
        self.assertEqual(g['epochs_per_region'], 128)
        self.assertEqual(g['bad_weight'], 796917)
        with self.assertRaises(ValueError):
            screen.geometry(12, Fraction(1, 10))

    def test_exact_selected_inner(self):
        a, kernel, _ = screen.load_inner()
        self.assertEqual(sum(a.values()), (1 << 19)-1)
        self.assertEqual(min(a), 48)
        self.assertEqual(next(w for w in range(1,129) if kernel[w]), 6)
        self.assertEqual(sum(kernel), 1 << 109)

    def test_q1_agrees_with_existing_rate_half_entry_point(self):
        a, _, _ = screen.load_inner()
        spectrum = screen.outer.read_spectrum(screen.HERE/'BCH256_64.wd')
        counts = {w:n for w,n in enumerate(spectrum) if w and n}
        tilts = np.array([-6., -5., -4.])
        g = screen.geometry(16, Fraction(1,10))
        moments = screen.q1.coefficient_logs(*screen.q1.region_logs(
            *screen.q1.epoch_logs(128,19,a,np.exp(tilts)), g['epochs_per_region']),256)
        actual = screen.aggregate_q1(moments,tilts,counts,g)
        # The old entry point's geometry already accepts arbitrary dimension;
        # only its distance target is hardcoded to 1/10.
        reference = screen.q1.screen(128,19,a,256,64,16,
                                    {w:math.log(n) for w,n in counts.items()},tilts)
        self.assertAlmostEqual(actual['margin_bits'],reference['margin_bits'],places=10)
        harder = screen.aggregate_q1(moments,tilts,counts,screen.geometry(16,Fraction(19,100)))
        self.assertLessEqual(harder['margin_bits'],actual['margin_bits'])

    def test_retained_partial_receipt(self):
        receipt = json.loads((screen.HERE/'FIXED_INNER_SCREEN.json').read_text())
        self.assertFalse(receipt['inner_reoptimized'])
        for path, digest in receipt['source_sha256'].items():
            self.assertEqual(screen.sha(screen.ROOT/path), digest, path)
        for row in receipt['results']:
            maximum = row['covered_occupations'][1]
            self.assertEqual([x['occupation'] for x in row['higher_occupations']],list(range(2,maximum+1)))
            self.assertEqual(row['output_bits'], 4*row['message_bits'])
            logs = [-x['margin_bits']*math.log(2) for x in [row['q1']]+row['higher_occupations']]
            self.assertAlmostEqual(-float(np.logaddexp.reduce(logs))/math.log(2),row['partial_union_margin_bits'],places=9)
            self.assertEqual(row['full_occupation_coverage'], maximum==row['outer_rows'])


if __name__ == '__main__':
    unittest.main()
