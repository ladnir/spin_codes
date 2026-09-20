import itertools
import math
import unittest
from fractions import Fraction as F
import activation_bridge as q1
from test_adaptive_range import exact_adaptive


def power(matrix,n):
    result=tuple(F(int(i==j)) for i in range(3) for j in range(3))
    for _ in range(n):result=q1.positive_mul(result,matrix)
    return result


class ConstantSplitTests(unittest.TestCase):
    def test_shifted_group_domination(self):
        ps=[F(1,3),F(2,3)];roots=[F(3,2),F(5,4)];n=4
        regions=[tuple(F((j+2)*(k+1)%13+1,17) for k in range(9)) for j in range(7)]
        for d in range(4):
            for h in range(4):
                shifted=regions[h:h+d+1]
                upper=exact_adaptive(shifted,ps,roots) if d else shifted[0]
                bound=sum(power(upper,n)[:3]);total=F(0)
                for labels in itertools.product(range(2),repeat=d):
                    matrix=[F(0)]*9
                    for bits in itertools.product((0,1),repeat=d):
                        mass=math.prod(ps[g] if bit else 1-ps[g] for g,bit in zip(labels,bits))
                        matrix=[a+mass*b for a,b in zip(matrix,regions[h+sum(bits)])]
                    value=math.prod(roots[g]**n for g in labels)*sum(power(matrix,n)[:3])
                    self.assertLessEqual(value,bound);total+=value
                self.assertLessEqual(total,2**d*bound)

    def test_complete_three_type_count(self):
        ordinary=14
        for length in range(1,8):
            total=sum(math.comb(length,d)*math.comb(length-d,h)*ordinary**d
                      for d in range(length+1) for h in range(length-d+1) if d+h)
            self.assertEqual(total,16**length-1)

    def test_constant_heavy_region_floor(self):
        length=16;covered=7;limit=3
        eligible=[(d,h) for d in range(limit+1) for h in range(length-d+1) if d+h>covered]
        self.assertTrue(all(h>=covered+1-limit for d,h in eligible))
        count=sum(math.comb(length,d)*math.comb(length-d,h)*5**d for d,h in eligible)
        upper=sum(math.comb(length,d)*2**(length-d)*5**d for d in range(limit+1))
        self.assertLessEqual(count,upper)


if __name__=='__main__':unittest.main()
