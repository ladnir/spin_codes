import math
import unittest
import csv
import tempfile
from pathlib import Path

import numpy as np

import activation_occupation as general
import balanced_occupation as balanced
import balanced_occupation_v2 as batched


class BalancedTest(unittest.TestCase):
    def test_density_dominates_every_shell(self):
        counts = {2: 17, 4: 36, 6: 17, 8: 1}
        for scale in (0., .2, .75, 1., 2.):
            for count in (1, 2, 4):
                bands, roots, lp, ln = balanced.density_roots(counts, 8, scale, count)
                for band, r, p, n in zip(bands, roots, lp, ln):
                    for w in band:
                        bound = 8*r+w*p+(8-w)*n if w < 8 else 8*r
                        actual = math.log(counts[w])-math.log(math.comb(8, w))
                        self.assertGreaterEqual(bound+1e-12, actual)

    def test_envelope_matches_all_measures_including_zero_entries(self):
        rng = np.random.default_rng(179)
        for size in (1, 3, 20, 100):
            roots = rng.uniform(-3, 3, size)
            ps = rng.uniform(0, 1, size); ps[-1] = 1.
            lp = np.log(ps)
            with np.errstate(divide='ignore'): ln = np.log1p(-ps)
            left = rng.uniform(-1000, 1000, (71, 3, 3))
            right = rng.uniform(-1000, 1000, left.shape)
            left[0] = -np.inf; right[1] = -np.inf
            left[2] = right[2] = -np.inf
            expected = np.full_like(left, -np.inf)
            for r, p, n in zip(roots, lp, ln):
                np.maximum(expected, r+np.logaddexp(n+left, p+right), out=expected)
            np.testing.assert_allclose(balanced.Envelope(roots, lp, ln).apply(left, right), expected, atol=5e-13)

    def test_recurrence_matches_unpruned_adaptive_bound(self):
        epoch = general.epoch_logs(8, 4, {4: 14, 8: 1}, [1,0,0,0,14,0,0,0,1], .7, 8)
        regions = general.region_logs(epoch, 8, 16, 12)
        counts = {2: 3, 4: 7, 8: 1}
        for scale in (0., .2, 1.):
            bands, roots, lp, ln = balanced.density_roots(counts, 8, scale)
            expected = []
            for q, matrix in general.adaptive_logs(regions, roots, lp, ln):
                value = math.log(math.comb(16, q))+q*math.log(len(bands))+12*.7+general.terminal_log(matrix, 8)
                expected.append(min(value, math.log(math.comb(16, q))+q*math.log(sum(counts.values()))))
            np.testing.assert_allclose(balanced.occupation_bounds(regions, counts, 8, 16, 12, .7, scale), expected, atol=1e-12)

    def test_batched_witnesses_and_region_ladder_preserve_bounds(self):
        epoch = general.epoch_logs(8,4,{4:14,8:1},[1,0,0,0,14,0,0,0,1],.7,8)
        regions = batched.region_ladder(epoch,8,[8,16,64],12)
        counts = {2:3,4:7,8:1}; scales = [0.,.2,1.]
        model = batched.PreparedModel(counts,8,scales)
        for length,region in regions.items():
            expected_region = general.region_logs(epoch,8,length,min(12,length))
            np.testing.assert_allclose(region,expected_region,atol=1e-12)
            values,witnesses = model.bounds(region,length,12,.7)
            expected = np.stack([balanced.occupation_bounds(region,counts,8,length,12,.7,s) for s in scales])
            np.testing.assert_allclose(values,expected.min(axis=0),atol=1e-12)
            np.testing.assert_allclose(witnesses,np.array(scales)[expected.argmin(axis=0)])

    def test_producer_replays_scalar_receipt_values(self):
        import run_balanced_occupation_grid as scalar
        import run_balanced_occupation_grid_v2 as vector
        _,observations,_ = scalar.coverage.snapshot()
        row = next(r for r in observations.values() if int(r['block_bits'])==512
                   and int(r['step_bits'])==64 and int(r['state_bits'])==12
                   and int(r['message_exponent'])==16 and r['outer_model']!='random-ensemble-average')
        arguments = dict(maximum_occupation=4,log_surprisals=[-5.,-3.],
                         probability_scales=[0.,.2],bands=0,setup_failure_bits=60)
        with tempfile.TemporaryDirectory(dir=scalar.grid.HERE) as temporary:
            root = Path(temporary)
            scalar.run_batch(root/'scalar',[row],arguments,{})
            vector.run_batch(root/'vector',[row],arguments,{})
            with (root/'scalar/occupations.csv').open() as handle:
                a = list(csv.DictReader(handle))
            with (root/'vector/occupations.csv').open() as handle:
                b = list(csv.DictReader(handle))
            for first,second in zip(a,b,strict=True):
                self.assertEqual(first['occupation'],second['occupation'])
                self.assertAlmostEqual(float(first['margin_bits']),float(second['margin_bits']),places=9)
                self.assertEqual(first['dominant_log_surprisal'],second['dominant_log_surprisal'])
                self.assertEqual(first['witness_probability_scale'],second['witness_probability_scale'])
            self.assertEqual(vector.verify(root/'vector',arguments)['row_count'],3)


if __name__ == '__main__': unittest.main()
