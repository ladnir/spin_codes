from itertools import permutations,product
from math import comb,log2
import unittest
import numpy as np

from hd_local_shuffle_control import shuffled_coefficients,path_measure,column_moment,first_active_words
from hd_post_refresh_quiet_rows import log_matrix_power,moment_atoms

MEASURE = dict(log2_probability=-3.,log2_atom00=-2.,log2_atom01=-1.,log2_atom10=-.5,
               log2_scale01=-2.,log2_scale10=-4.,log2_scale11=-1.)


class LocalShuffleTests(unittest.TestCase):
    def test_exact_binary_toy_domination(self):
        # H=I+11^T in dimension4, A={0,1}; outside rank1, zero fiber d1.
        coeff = shuffled_coefficients((1,0,1),(1,0,1),1,size=4,q=2)
        for t in ((0.,.2,.8,1.),(.1,)*4,(1.,)*4):
            true = 0.
            for x0,x1 in product(range(2),repeat=2):
                x = (x0,x1,0,0)
                y = tuple(v^x0^x1 for v in x)
                true += sum(np.prod([t[p[i]] for i,v in enumerate(y) if v])
                            for p in permutations(range(4)))/24
            upper = sum(float(coeff[w])*sum(np.prod([t[i] for i in range(4) if mask>>i&1])
                for mask in range(16) if mask.bit_count()==w) for w in range(5))
            self.assertLessEqual(true,upper+1e-12)

    def test_count_dp_mass(self):
        for length in (1,2,4,16):
            paths,elementary,empty = path_measure(MEASURE,length=length)
            self.assertEqual(len(paths),comb(length+4,4))
            self.assertEqual(elementary.shape,(comb(length+4,4),length+1))
            np.testing.assert_allclose(np.logaddexp2.reduce(paths,axis=0),empty,atol=1e-11)
            np.testing.assert_allclose(empty,log_matrix_power(moment_atoms(MEASURE)[1],length),atol=1e-11)

    def test_constant_kernel_is_empty_measure(self):
        paths,elementary,empty = path_measure(MEASURE,length=4)
        from fractions import Fraction
        coeff = tuple(map(Fraction,(1,0,0,0,0)))
        np.testing.assert_allclose(column_moment(coeff,paths,elementary,18),empty,atol=1e-11)
        self.assertTrue(np.isneginf(column_moment(coeff,paths,elementary,18,nonzero=True)).all())

    def test_first_activity_at_power_one(self):
        def matrix(value):
            answer = np.full((1,2,2),-np.inf)
            answer[:,0,0] = answer[:,1,1] = log2(value)
            return answer
        active = .7
        got = first_active_words(matrix(1+active),matrix(active),matrix(1),1,columns=4,minimum=2)[0]
        expected = log2(sum(comb(4,i)*active**i for i in range(2,5)))
        self.assertAlmostEqual(got,expected,places=12)

    def test_monotonicity_scope_guard(self):
        with self.assertRaises(ValueError):
            path_measure(dict(MEASURE,log2_atom10=.1))


if __name__=='__main__':
    unittest.main()
