import io
import unittest
from contextlib import redirect_stdout
from math import comb
from fractions import Fraction as Q
from unittest.mock import patch
import numpy as np
from flint import arb,arb_mat,ctx
import sparse_cover
import occupancy_cdf_cover as cover


class OutwardFoldCacheTests(unittest.TestCase):
    def test_grouped_outward_masses_contain_exact_coefficients(self):
        ctx.prec=256;D=cover.DENOMINATOR
        for nums in ([],[0]*20,[D]*20,[D//2]*96,[D//4]*60+[3*D//4]*36,
                     [1,D-1,D//3,2*D//3]*5):
            expected=[Q(1)]
            for numerator in nums:
                p=Q(numerator,D);next_=[Q(0)]*(len(expected)+1)
                for j,m in enumerate(expected):next_[j]+=m*(1-p);next_[j+1]+=m*p
                expected=next_
            actual=cover.bernoulli_masses(nums)
            self.assertEqual(len(actual),len(expected))
            for value,exact in zip(actual,expected):
                lo,e=value.lower().man_exp();lower=Q(int(lo))*Q(2)**int(e)
                hi,e=value.upper().man_exp();upper=Q(int(hi))*Q(2)**int(e)
                self.assertLessEqual(lower,exact);self.assertGreaterEqual(upper,exact)
        for nums in ([-1],[D+1],[.5]):
            with self.assertRaises(ValueError):cover.bernoulli_masses(nums)

    def test_identical_group_folds_are_evaluated_once(self):
        ctx.prec=256
        args=sparse_cover.build_args(2,['.0032'],256,0,52,None)
        args.joint_witness=False
        counts=[max(0,u-37) for u in range(257)]
        exact=[arb_mat([[arb(2)**-100]]) for _ in range(3)]
        floating=[np.array([[2.**-100]]) for _ in range(3)]
        with patch.object(cover,'fold_arb',wraps=cover.fold_arb) as folds,redirect_stdout(io.StringIO()):
            upper=cover.cover(args,{('.0032','1'):(exact,floating)},{'1':counts},[1],cutoff=0)
        self.assertIsNotNone(upper)
        self.assertEqual(folds.call_count,1)
        _,lo,hi,p=folds.call_args.args
        expected=arb(comb(2048,2))*arb(2)**-25600*cover.fold_arb(counts,lo,hi,p)**2
        self.assertLess(abs(float(upper/expected)-1),1e-14)


if __name__=='__main__':unittest.main()
