import itertools
import unittest
from fractions import Fraction as F
import numpy as np
from flint import arb, ctx
import scaled_adaptive as scaled
from certify_larger_state_range import costs_for
from test_adaptive_range import exact_adaptive
import bridge as base
import larger_state_maps as maps


class ExactConstantBandTests(unittest.TestCase):
    def test_selected_generators_are_quadratic_evaluations(self):
        for name in maps.NAMES:
            selected = base.read(maps.DIRECTORY/f'{name}_selection.json')
            t = selected['parameters']['step_bits']
            variables = t.bit_length()-1
            self.assertEqual(t, 1 << variables)
            for encoded in selected['selected']['A_generator_words_hex']:
                word = int(encoded, 16)
                coefficients = [(word >> i) & 1 for i in range(t)]
                for bit in range(variables):
                    for mask in range(t):
                        if mask & (1 << bit):
                            coefficients[mask] ^= coefficients[mask ^ (1 << bit)]
                self.assertTrue(all(mask.bit_count() <= 2 for mask,c in enumerate(coefficients) if c))

    def test_cost_one_only_for_all_one_word(self):
        self.assertEqual(costs_for([(256,)], [F(1)], {256:1}), [F(1)])
        with self.assertRaises(AssertionError):
            costs_for([(254,256)], [F(1)], {254:1,256:1})
        with self.assertRaises(AssertionError):
            costs_for([(256,)], [F(1)], {256:2})

    def test_adaptive_dominates_every_fixed_band_sequence(self):
        ctx.prec = 256
        region = [tuple(F((j+2)*(k+1), 71) for k in range(9)) for j in range(6)]
        ps = [F(1,4), F(2,3), F(1)]
        roots = [F(5,4), F(7,6), F(1)]
        for q in range(1,6):
            expected = exact_adaptive(region[:q+1], ps, roots)
            for labels in itertools.product(range(3), repeat=q):
                law = [F(1)]
                cost = F(1)
                for label in labels:
                    p = ps[label]
                    new = [F(0)]*(len(law)+1)
                    for j,v in enumerate(law):
                        new[j] += (1-p)*v
                        new[j+1] += p*v
                    law = new
                    cost *= roots[label]
                actual = [cost*sum((v*region[j][k] for j,v in enumerate(law)), F(0)) for k in range(9)]
                self.assertTrue(all(a <= b for a,b in zip(actual, expected)))
        upper = [tuple((arb(v.numerator)/v.denominator).upper() for v in row) for row in region]
        mantissas, exponents = scaled.initial(upper)
        left = np.nextafter(np.array([float(r*(1-p)) for r,p in zip(roots,ps)]), np.inf)
        right = np.nextafter(np.array([float(r*p) for r,p in zip(roots,ps)]), np.inf)
        for q,matrix,exponent in scaled.matrices(mantissas,exponents,left,right):
            expected = exact_adaptive(region[:q+1], ps, roots)
            self.assertTrue(all(F.from_float(float(v))*F(2)**exponent >= w for v,w in zip(matrix.flat,expected)))


if __name__ == '__main__':
    unittest.main()
