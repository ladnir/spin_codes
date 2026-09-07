from fractions import Fraction as F
import unittest

from flint import arb, ctx

import fine_split_search as fine


class FineProposalTests(unittest.TestCase):
    def test_fractional_proposal_at_narrow_window(self):
        powers = [[-80]*(q+1) for q in range(231, 252)]
        for row in powers:
            row[0] = 70
        index, q, h = fine.propose(powers, 231, 251, [])
        self.assertEqual((index, q, h), (F(-5, 2), 241, 0))

    def test_complete_stops(self):
        self.assertIsNone(fine.propose([[-80]*242], 241, 241, []))

    def test_integer_bank_conversion_preserves_probabilities(self):
        row = dict(tilt=-3, anchor_occupation=200, anchor_all_one_rows=0, p=['witness'])
        converted = fine.normal_shards(dict(status='CONSTANT_SPLIT_GRID_OUTWARD_PRODUCER', shards=[row]))
        self.assertEqual(F(converted[0]['tilt_numerator'], converted[0]['tilt_denominator']), F(-3))
        self.assertEqual(converted[0]['p'], row['p'])

    def test_fractional_tilt_closes_actual_q241(self):
        ctx.prec = 256
        spec = fine.core.instance('t128_s19', 16)
        t, s, ac, kernel = fine.core.inputs.load('t128_s19')
        caps = fine.core.inputs.caps_module.caps()
        for index in (F(-3), F(-2), F(-5, 2)):
            tilt = arb(index.numerator)/index.denominator
            region = fine.density.engine.regions(t, s, ac, kernel, (-(tilt/10).exp()).exp(), 241, 512)
            ps = fine.density.choose(region, 241, caps)[:-1]
            powers = fine.grid.evaluate(region, ps, spec, 241, 241, tilt)
            bound = fine.grid.bounds(powers, 241, 241)[0]
            if index.denominator == 1:
                self.assertGreater(bound, F(2)**-40)
            else:
                self.assertLess(bound, F(2)**-70)


if __name__ == '__main__':
    unittest.main()
