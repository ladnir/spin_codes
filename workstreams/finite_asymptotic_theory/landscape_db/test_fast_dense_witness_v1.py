import math
import unittest

import numpy as np

import fast_dense_witness_v1 as fast


class FastWitnessTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = fast.FastCoverRefiner({2: 2, 4: 1}, 4, 4, 2, {2: 2, 4: 1},
                                          [1, 0, 2, 0, 1], 8, [[2], [4]])

    def test_table_derivatives_inside_cell(self):
        x, y, h = -2.113, -1.271, 1e-6
        value, dx, dy = self.model.lookup(x, y)
        self.assertAlmostEqual(dx, (self.model.lookup(x+h, y)[0]-self.model.lookup(x-h, y)[0])/(2*h), places=8)
        self.assertAlmostEqual(dy, (self.model.lookup(x, y+h)[0]-self.model.lookup(x, y-h)[0])/(2*h), places=8)

    def test_retained_witness_uses_direct_transfer(self):
        p = np.array([0., .5, 1.]); proposal = np.array([.5, .3, .2])
        witness = dict(log_surprisal=-1., proposal=proposal.tolist(), probabilities=p.tolist(),
                       log_density_costs=self.model.costs(p).tolist())
        lower, upper = [0, 0, 0], [5, 8, 8]
        value = self.model.direct(lower, upper, witness)
        box = dict(lower=lower, upper=upper, own_log_bound=value, witness=witness)
        result = self.model.refine_box(box, -20*math.log(2))
        expected = self.model.direct(lower, upper, result['witness'])
        self.assertAlmostEqual(expected, result['own_log_bound'], places=11)
        self.assertLessEqual(expected, value+1e-11)
        self.assertLess(fast.FastCoverRefiner.__mro__.index(fast.FastRefiner),
                        fast.FastCoverRefiner.__mro__.index(fast.base.Refiner))


if __name__ == '__main__':
    unittest.main()
