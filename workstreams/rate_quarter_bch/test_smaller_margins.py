"""Exact algebra and retained-result checks for the smaller quarter-rate outer."""
import json
import math
import unittest

import numpy as np
import smaller_outer as small
import evaluate_smaller_margins as evaluator
import refine_smaller_dense as refinement


class SmallerOuterTests(unittest.TestCase):
    def test_prepared_epochs(self):
        a,kernel,_ = evaluator.fixed.load_inner()
        prepared = refinement.Epochs(a,kernel)
        for z in (-18.,-8.3,-1.1,2.):
            np.testing.assert_allclose(prepared.at(z),evaluator.fixed.general.epoch_logs(128,19,a,kernel,math.exp(z),128),atol=2e-11,rtol=2e-13)

    def test_closed_targets(self):
        path = small.r.HERE/'SMALLER_MARGIN_CLOSED.json'
        if not path.exists():
            self.skipTest('Run dense refinement and sparse extension first')
        receipt = json.loads(path.read_text())
        for name,digest in receipt['source_sha256'].items():
            self.assertEqual(evaluator.fixed.sha(evaluator.fixed.ROOT/name),digest,name)
        targets = {'33/200':40.,'19/100':30.}
        self.assertEqual({r['distance_target'] for r in receipt['results']},set(targets))
        for row in receipt['results']:
            self.assertGreaterEqual(row['combined_margin_bits'],targets[row['distance_target']])
            self.assertTrue(row['exact_integer_coverage_checked'])

    def test_construction(self):
        c = small.construction()
        self.assertEqual((c['length'],c['dimension'],c['distance_lower_bound']),(128,32,32))
        self.assertEqual(c['nonzero_coset_orbit'],127)
        self.assertEqual(c['selected_nonzero_cosets'],7)
        self.assertEqual(c['parent_generator_hex'],'0xd5306d6bfdbc8574719e70d')
        self.assertEqual(c['subcode_generator_hex'],'0x5106dae17a61e520c606a4e29')

    def test_spectrum(self):
        counts = small.spectrum()
        self.assertEqual(counts[32],588)
        self.assertEqual(counts[36],896)
        self.assertEqual(sum(counts),1<<32)
        self.assertEqual(counts,counts[::-1])
        self.assertEqual(small.r.audit_spectrum(counts,32,32)['nonzero_coefficients'],19)

    def test_retained_outer(self):
        receipt = json.loads((small.r.HERE/'SMALLER_OUTER_AUDIT.json').read_text())
        self.assertEqual(receipt['construction'],small.construction())
        self.assertEqual(receipt['spectrum'],{str(w):str(n) for w,n in enumerate(small.spectrum()) if n})

    def test_retained_margin_receipt(self):
        path = small.r.HERE/'SMALLER_MARGIN_SCREEN.json'
        if not path.exists():
            self.skipTest('Run the smaller margin producer first')
        receipt = json.loads(path.read_text())
        self.assertFalse(receipt['inner_reoptimized'])
        self.assertEqual(receipt['construction'],small.construction())
        for name,digest in receipt['source_sha256'].items():
            self.assertEqual(evaluator.fixed.sha(evaluator.fixed.ROOT/name),digest,name)
        for row in receipt['results']:
            lo,hi = row['covered_occupations']
            self.assertEqual(lo,1)
            self.assertEqual(len(row['occupation_margins_bits']),hi)
            sparse = float(np.logaddexp.reduce(-np.array(row['occupation_margins_bits'])*math.log(2)))
            self.assertAlmostEqual(-sparse/math.log(2),row['sparse_union_margin_bits'],places=10)
            self.assertEqual(row['output_bits'],4*row['message_bits'])
            if 'dense' in row:
                self.assertEqual(row['dense']['occupation_min'],hi+1)
                evaluator.check_coverage(row['dense']['selected_boxes'],row['outer_rows'],hi+1)
                total = float(np.logaddexp(sparse,row['dense']['log_union_upper']))
                self.assertAlmostEqual(-total/math.log(2),row['combined_margin_bits'],places=9)


if __name__ == '__main__':
    unittest.main()
