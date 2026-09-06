"""Exact score identities and regressions for the spectrum-curve assessment."""
from fractions import Fraction as F
import unittest
import curve_spectrum_calibration as study
import frontier_ledger as ledger
import bridge as base


class CurveSpectrumTest(unittest.TestCase):
    def test_reference_mass_and_algebraic_properties(self):
        for block,d,even,one,path,counts in study.outer_inputs():
            ref=study.reference_counts(block,block//2,d,even,one)
            self.assertEqual(sum(ref.values(),F(0)),1<<(block//2))
            self.assertEqual(ref[0],1)
            if even:self.assertTrue(all(w%2==0 for w in ref))
            if one:self.assertEqual(ref[block],1)
            if block==8:self.assertEqual(ref,counts)
            if block==64:self.assertGreater(ref[13],0)
        self.assertEqual(study.reference_counts(256,128,38,True,True),study.model.even_binomial())

    def test_exact_conditional_endpoint_reconstruction(self):
        data=base.read(base.HERE/'generated/curve_spectrum_calibration_v1.json')
        for name,digest in data['source_sha256'].items():self.assertEqual(base.sha(base.ROOT/name),digest)
        rebuilt=study.conditional_endpoint(base.HERE/'generated/curve_k28_full_retained_v1.json',
            base.HERE/'generated/frontier_k28_q1_v1.json')
        self.assertEqual(rebuilt,data['conditional_endpoint'])
        ref=base.decode(rebuilt['reference_q1_upper']);higher=base.decode(rebuilt['higher_upper'])
        for row in rebuilt['cases']:
            self.assertEqual(base.decode(row['conditional_upper']),ref*(1<<row['weighted_score_inflation_bits'])+higher)
        self.assertLess(rebuilt['cases'][0]['higher_occupancy_penalty_bits'],.054)

    def test_retained_tree_partitions_and_budget(self):
        for tag in ('6144_8191','8192_524287','524288_1572864','1572865_2097152'):
            data=base.read(base.HERE/f'generated/curve_k28_dense_{tag}_100_v1.json')
            old=base.read(base.ROOT/data['seed_certificate']);lo,hi=data['occupancy_range']
            self.assertEqual(ledger.validate_tree(data['tree'],lo,hi),data['leaves'])
            self.assertGreaterEqual(data['leaves'],old['leaves'])
            self.assertEqual(base.decode(data['range_upper']),F(hi-lo+1,1<<100))
            self.assertEqual(base.decode(old['range_upper']),base.decode(data['range_upper'])*(1<<20))

    def test_backtest_labels_and_error_identity(self):
        data=base.read(base.HERE/'generated/curve_spectrum_calibration_v1.json')
        known=[r for r in data['rows'] if r['exact_spectrum_q1_bits'] is not None]
        self.assertEqual(len(known),16)
        for row in known:
            self.assertAlmostEqual(row['weighted_inflation_bits'],row['model_q1_bits']-row['exact_spectrum_q1_bits'])
        self.assertLess(data['largest_observed_weighted_inflation_bits'],1)
        for row in data['rows']:
            if row['block_bits']==256:self.assertIsNone(row['exact_spectrum_q1_bits'])


if __name__=='__main__':unittest.main()
