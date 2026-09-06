import itertools
import math
import unittest
from fractions import Fraction as F
import activation_bridge as q1
import exponential_mode_bound as mode


def product_matrices(matrices):
    value=tuple(F(int(i==j)) for i in range(3) for j in range(3))
    for matrix in matrices:value=q1.positive_mul(value,matrix)
    return value


class ModeTests(unittest.TestCase):
    def test_maclaurin_all_weights(self):
        for values in itertools.product((F(0),F(1,3),F(1),F(5,3)),repeat=4):
            mean=sum(values)/4
            for w in range(5):
                exact=sum((math.prod(values[i] for i in subset) for subset in itertools.combinations(range(4),w)),F(0))/math.comb(4,w)
                self.assertLessEqual(exact,mean**w)

    def test_noncommutative_pattern_coefficients(self):
        modes=[tuple(F((i+2)*(k+1)%7+1,20) for k in range(9)) for i in range(3)]
        coefficients=mode.polynomial_power(modes,4)
        expected=[[F(0)]*9 for _ in range(9)]
        for labels in itertools.product(range(3),repeat=4):
            matrix=product_matrices([modes[i] for i in labels])
            expected[sum(labels)]=[a+b for a,b in zip(expected[sum(labels)],matrix)]
        self.assertEqual(coefficients,[tuple(row) for row in expected])

    def test_full_fixed_weight_row_count(self):
        # A [4,2] code has two weight-2 words and one weight-4 word.
        # After row permutation, each weight-2 support has counting mass 1/3.
        n=4;q=2;caps={2:2,4:1};a=F(1,10);delta=F(3,10)
        modes=[tuple(F((i+2)*(k+1)%7+1,30) for k in range(9)) for i in range(3)]
        regions=[tuple(sum((modes[i][k]*(a+i*delta)**j for i in range(3)),F(0)) for k in range(9)) for j in range(q+1)]
        supports=[(mask,F(caps[mask.bit_count()],math.comb(n,mask.bit_count()))) for mask in range(1<<n) if mask.bit_count() in caps]
        exact=F(0)
        for rows in itertools.product(supports,repeat=q):
            counts=[sum((mask>>j)&1 for mask,_ in rows) for j in range(n)]
            moment=sum(product_matrices([regions[j] for j in counts])[:3])
            exact+=math.prod(mass for _,mass in rows)*moment
        upper=mode.numerator(modes,a,delta,n,q,caps)
        self.assertLessEqual(exact,upper)
        # Before replacing a fixed-weight symmetric mean by an arithmetic
        # mean power, the expansion agrees with exhaustive row enumeration.
        expanded=F(0)
        for labels in itertools.product(range(3),repeat=n):
            rates=[a+i*delta for i in labels]
            row_cost=sum((count*sum((math.prod(rates[i] for i in subset) for subset in itertools.combinations(range(n),w)),F(0))/math.comb(n,w)
                          for w,count in caps.items()),F(0))
            expanded+=sum(product_matrices([modes[i] for i in labels])[:3])*row_cost**q
        self.assertEqual(exact,expanded)


if __name__=='__main__':unittest.main()
