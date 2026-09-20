"""High-precision spot checks supplement, but do not replace, the error proof."""
import unittest
from fractions import Fraction as F
import numpy as np
from flint import arb,acb,ctx
from screen_pair_type_bound import load_pairs,SIGNS
from certify_pair_type_fourier import roots
from pair_fourier_error import evaluate


class ErrorTests(unittest.TestCase):
    def test_against_true_roots_and_high_precision_evaluation(self):
        ctx.prec=512;powers,counts=load_pairs();grid=[128,128,64]
        p=[F(805,1000),F(92,1000),F(92,1000),F(11,1000)]
        # Use exactly representable positive masses summing to one.
        den=1<<40;nums=[round(float(v)*den) for v in p[1:]];nums=[den-sum(nums)]+nums
        ps=[F(v,den) for v in nums];floats=np.array([float(v) for v in ps])
        angles=[(0,0,0),(1,0,0),(0,1,0),(0,0,1),(64,0,32),(0,64,32),(64,64,0),(1,1,1)]
        rng=np.random.default_rng(4112);angles += [tuple(map(int,v)) for v in rng.integers([0,0,0],grid,size=(40,3))]
        tables=[roots(n) for n in grid]
        values=np.array([[floats[0]]*len(angles)]+[[floats[k+1]*tables[k][angle[k]] for angle in angles] for k in range(3)])
        computed,error=evaluate(values,powers,counts)
        for index,angle in enumerate(angles):
            exact=[acb(arb(ps[0].numerator)/ps[0].denominator)]
            for k in range(3):
                argument=2*arb.pi()*angle[k]/grid[k]
                exact.append((arb(ps[k+1].numerator)/ps[k+1].denominator)*acb(argument.cos(),argument.sin()))
            linear=[sum((int(sign)*value for sign,value in zip(row,exact)),acb(0)) for row in SIGNS]
            result=acb(0)
            for exponents,count in zip(powers,counts):
                term=acb(float(count))
                for value,e in zip(linear,exponents):term*=value**(4*int(e))
                result+=term
            distance=abs(result-acb(float(computed[index].real),float(computed[index].imag))).upper()
            self.assertLessEqual(distance,arb(float(error[index])))


if __name__=='__main__':unittest.main()
