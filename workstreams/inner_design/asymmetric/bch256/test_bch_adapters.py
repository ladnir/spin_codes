"""Regression checks for length parameters, directed folds, and cover mapping."""
from fractions import Fraction as F
import unittest
from flint import arb, ctx
import bch_model as model
import sparse_ranges as ranges
import certify_dense as dense


class AdapterTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        ctx.prec = 256
        cls.engine = model.Engine(16)

    def test_regions_zero_one(self):
        engine = self.engine
        tilt = F(-8)
        zero, one = engine.epoch(tilt)[:2]
        reference = model.independent.single.regions(zero, one, engine.n, engine.length//128)
        actual = engine.region(tilt, 1)
        for a, b in zip(actual, reference):
            for x, y in zip(a, b):
                self.assertLess(abs(x-y), arb(2)**-200)

    def test_linear_epoch_q1(self):
        ordinary = self.engine.q1(F(-8))
        linear = self.engine.q1(F(-8), linear=True)
        self.assertEqual(len(ordinary), 257)
        for a, b in zip(ordinary, linear):
            self.assertLess(abs(a-b), arb(2)**-180)

    def test_directed_folds_dominate_arb(self):
        engine = self.engine
        tilt = F(-4)
        ps = [F(sum(band), len(band)*256) for band in model.BANDS]
        region = engine.region(tilt, 6)
        rows = ranges.evaluate(engine, region, ps, tilt, 2, 6)
        for row in rows:
            q = row['occupation']
            actual = engine.adaptive(region[:q+1], ps, q)*(engine.cutoff*model.number(tilt).exp()).exp()
            self.assertLessEqual(actual, arb(2)**row['power'])

    def test_geometry_rescaling(self):
        seed = dict(leaves={'':dense.geometry.box(512, 8192, F(0), F(1), None)}, splits={})
        result = dense.mapped_leaves(seed, 18, 128)
        self.assertEqual(result, [('', 128, 2048, F(0), F(1), None)])
        truncated = dense.mapped_leaves(seed, 16, 256)
        self.assertEqual(truncated[0][1:3], (256, 512))

    def test_missing_partition_rejected(self):
        seed = dict(leaves={'0':dense.geometry.box(512, 4352, F(0), F(1), None)}, splits={'':'q'})
        with self.assertRaises(ValueError):
            dense.mapped_leaves(seed, 20, 512)


if __name__ == '__main__':
    unittest.main()
