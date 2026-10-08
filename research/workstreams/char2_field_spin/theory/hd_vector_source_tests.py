from fractions import Fraction as F
from math import comb
from itertools import product
import unittest

import numpy as np
from scipy.special import logsumexp

from hd_vector_source import (exact_row_degree,log_row_degree,exact_all_row_filter,
                              directed_ceiling,return_rank_charge)
from hd_vector_region import SCREEN_WITNESSES
from hd_vector_region_tests import rank,mul


class VectorSourceTests(unittest.TestCase):
    def test_distinguished_row_against_complete_support_enumeration(self):
        for active,length,rows in ((1,4,3),(2,4,3),(3,3,4)):
            for tau in (F(1),F(2),F(4)):
                total,good,degrees=exact_all_row_filter(active,length,rows,F(2,5),4,2,tau,length//2+1)
                self.assertEqual(degrees,exact_row_degree(active,length,rows,F(2,5),4,2,tau))
                bad=sum(degrees[:length//2+1])
                self.assertGreaterEqual(good,max(F(0),total-active*bad))
                self.assertLessEqual(good,total-bad)

    def test_logarithmic_polynomial_matches_exact(self):
        expected=exact_row_degree(2,4,3,F(2,5),4,2,F(2))
        actual=log_row_degree(2,4,3,.4,4,2,2.)
        np.testing.assert_allclose(np.exp(actual),list(map(float,expected)),rtol=2e-13,atol=1e-13)

    def test_mds_full_field_value_count_loss_bound(self):
        for b,k,q in ((4,2,4),(8,4,16),(16,8,256)):
            d=b-k+1
            for w in range(d,b+1):
                actual=(q-1)*sum((-1)**j*comb(w-1,j)*q**(w-d-j) for j in range(w-d+1))
                self.assertLessEqual(actual,q**(w-b+k))
                self.assertGreaterEqual(actual,F(q**(w-b+k))*(1-F(w,q)))

    def test_directed_row_filter_ceiling(self):
        for b,k,a,is_rank,point in SCREEN_WITNESSES:
            if is_rank and b==128 and a in (16,64) and len(point)==2:
                result=directed_ceiling(a,point)
                self.assertTrue(result['rigorous'])
                self.assertFalse(result['actual_failure_probability_claim'])
                self.assertGreaterEqual(result['filtered_proof_expression_lower_log2'],500 if a==16 else 5000)

    def test_conditional_return_rank_formula(self):
        u=(1,0,0);v=(0,1,0);x=(0,0,1);z=(0,0,0)
        for s in ((z,z),(x,x),(u,x),(u,v),(u,v,(1,1,0))):
            h=rank(s+(x,))-1
            for tau in (F(1),F(2),F(4)):
                total=F(0)
                for d in product(range(4),repeat=len(s)):
                    equations=tuple(tuple(a^mul(b,c) for a,c in zip(row,x)) for row,b in zip(s,d))
                    rr=rank(equations)
                    if rank(equations+(x,))>rr:total+=tau**(-rr)
                self.assertEqual(total/4**len(s),return_rank_charge(h,4,len(s),tau))


if __name__=='__main__':unittest.main()
