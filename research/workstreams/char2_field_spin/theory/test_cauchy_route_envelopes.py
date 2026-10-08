"""Exact checks for the two exploratory outer-route envelopes."""

from fractions import Fraction
from itertools import product
from math import comb, prod
import unittest


def expectation(probabilities, function):
    return sum((function(sum(bits))*
                prod((p if bit else 1-p for bit,p in zip(bits,probabilities)),
                     start=Fraction(1))
                for bits in product((0,1), repeat=len(probabilities))), Fraction(0))


class RouteEnvelopes(unittest.TestCase):
    def test_count_function_discrete_convexity(self):
        for parity in (1,2,10):
            for tau in (Fraction(1), Fraction(14), Fraction(255)):
                for r in (Fraction(2), Fraction(3), Fraction(511,2)):
                    f = [Fraction(1)]+[Fraction(0)]*parity+[
                        r**w/tau**parity for w in range(parity+1,33)]
                    self.assertTrue(all(f[w]+f[w+2]>=2*f[w+1]
                                        for w in range(len(f)-2)))

    def test_pairwise_balancing_increases_expectation(self):
        def f(w):
            return Fraction(1) if w==0 else Fraction(0) if w<=2 else Fraction(3**w,49)
        for probabilities in product((Fraction(1,10),Fraction(1,3),Fraction(3,4)),repeat=4):
            mean=(probabilities[0]+probabilities[1])/2
            balanced=(mean,mean)+probabilities[2:]
            self.assertLessEqual(expectation(probabilities,f),expectation(balanced,f))

    def test_exact_selected_row_probability(self):
        for rows in range(1,21):
            for active in range(1,rows+1):
                for weight in range(rows+1):
                    exact=Fraction(comb(active,weight) if weight<=active else 0,
                                   comb(rows,weight))
                    self.assertLessEqual(exact,Fraction(active,rows)**weight)

    def test_conditioning_is_not_needed_for_selected_row_event(self):
        # Each column draws an independent fixed-size subset.  The event
        # that all slots are routed from selected rows factors exactly.
        from itertools import combinations
        rows, active, counts=5,3,(1,2,3)
        possibilities=[list(combinations(range(rows),c)) for c in counts]
        total=success=0
        for columns in product(*possibilities):
            total+=1
            success+=all(all(row<active for row in column) for column in columns)
        formula=Fraction(1)
        for count in counts:
            formula*=Fraction(comb(active,count),comb(rows,count))
        self.assertEqual(Fraction(success,total),formula)


if __name__=='__main__':
    unittest.main()
