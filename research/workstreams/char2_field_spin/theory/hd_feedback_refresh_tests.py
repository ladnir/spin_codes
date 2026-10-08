"""Exact small-field check of the new feedback-refresh density transfer."""

from fractions import Fraction as F
from itertools import combinations,product
import unittest

from hd_blockwise_shared_transfer import multiply


class FeedbackRefreshTests(unittest.TestCase):
    def test_small_field_pointwise_density(self):
        # F=GF4, state dimension2. F_packet=I, A=[I;I], C=[I I].
        q,d,z=4,15,F(2,3)
        live=tuple(h for h in product(range(q),repeat=2) if any(h))
        all_states=((0,0),)+live
        a,b,c=F(1,5),F(1,20),F(3,4)
        mixer={}
        for h in live:
            line={tuple(multiply(scale,value,q) for value in h) for scale in range(1,q)}
            mixer[h]={target:a*(target==h)+b*F(target in line,q-1)+c/d for target in live}
            self.assertEqual(sum(mixer[h].values()),1)
        for size in range(5):
            for support in combinations(range(4),size):
                zero_law={target:F(0) for target in all_states}
                flat_law={target:F(0) for target in all_states}
                quiet=free=flat=flat_quiet=F(0)
                for values in product(range(q),repeat=size):
                    x=[0]*4
                    for coordinate,value in zip(support,values): x[coordinate]=value
                    syndrome=(x[0]^x[2],x[1]^x[3])
                    weight=z**sum(value!=0 for value in x)/q**size
                    free+=weight
                    if any(syndrome):
                        for target in live: zero_law[target]+=weight/d
                    else:
                        quiet+=weight
                        zero_law[(0,0)]+=weight
                    for h in live:
                        y=[x[i]^h[i%2] for i in range(4)]
                        tilted=z**sum(value!=0 for value in y)/q**size/d
                        flat+=tilted
                        if any(syndrome):
                            for target in all_states:
                                flat_law[target]+=tilted*(1-mixer[h].get(target,F(0)))/d
                        else:
                            flat_quiet+=tilted
                            for target in live: flat_law[target]+=tilted*mixer[h][target]
                self.assertEqual(zero_law[(0,0)],quiet)
                for target in live:
                    self.assertEqual(zero_law[target],(free-quiet)/d)
                # The affine-coset quiet column is bounded by the zero coset.
                bulk=flat-F(1,q)*flat_quiet+F(1,q)*quiet
                self.assertEqual(flat_law[(0,0)],(flat-flat_quiet)/d)
                for target in live:
                    self.assertLessEqual(flat_law[target],bulk/d)


if __name__=='__main__':
    unittest.main()
