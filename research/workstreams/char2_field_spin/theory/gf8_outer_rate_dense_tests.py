import unittest
from fractions import Fraction as F

from flint import arb,ctx

import gf8_outer_rate_dense as dense
from field128_dense_certificate import interval_base,TILTS


class DenseRateTests(unittest.TestCase):
    def test_geometry(self):
        g=dense.geometry(88,20)
        self.assertEqual((g['rows'],g['padded_k'],g['n'],g['output_threshold']),
                         (11968,1053184,1531904,153190))
        self.assertLess(g['padded_k']-g['requested_k'],64*88)
        self.assertEqual({dense.geometry(d,16)['n'] for d in range(87,94)},{98304})

    def test_frozen_bracket_and_monotonicity(self):
        for index in (48,67,130,199):
            tilt=TILTS[index-48]
            for dimension in (87,88,90,96):
                bound=dense.bracket(dimension,F(index,200),F(index+1,200),tilt)
                self.assertEqual(bound,interval_base(128,dimension,index,tilt))
                for alpha in (F(index,200),F(2*index+1,400),F(index+1,200)):
                    self.assertLessEqual(dense.bracket(dimension,alpha,alpha,tilt),bound)
            self.assertLess(dense.bracket(88,F(index,200),F(index+1,200),tilt),
                            dense.bracket(90,F(index,200),F(index+1,200),tilt))

    def test_one_cell_exact_integer_comparison(self):
        # Independently compare the directed arithmetic with an exact
        # integer-power inequality in a one-row instance of the dense lemma.
        # Production padding is tested separately; no receipt uses this toy.
        g=dict(outer=[128,88],rows=1,n=128)
        weight=g['n']//2
        tilt=F(2)
        base=dense.bracket(88,F(1,2),F(1,2),tilt)
        old=ctx.prec
        try:
            ctx.prec=192
            bits=dense.lower_integer(dense.interval_margin(g,weight,weight,tilt))
            self.assertLessEqual((g['rows']+1)**(128*100)*base**(g['rows']*100)
                                 *2**(127*g['n']+100*bits),tilt**(weight*100))
        finally:ctx.prec=old

    def test_dense_replay_and_adjacent_coverage(self):
        receipt=dense.certify(88,20,include_witnesses=True)
        self.assertTrue(receipt['complete'] and receipt['rigorous'])
        self.assertFalse(receipt['whole_code_certified'])
        self.assertEqual(receipt['weights'],[367657,1531904])
        self.assertGreaterEqual(receipt['margin_bits_lower'],64)
        entries=receipt['witnesses']
        self.assertTrue(all(a[1]+1==b[0] for a,b in zip(entries,entries[1:])))
        self.assertEqual(sum(hi-lo+1 for lo,hi,_,_ in entries),receipt['covered_weights'])

    def test_failed_class_and_invalid_arguments(self):
        with self.assertRaises(ArithmeticError):dense.certify(96,20,F(6,25))
        with self.assertRaises(ValueError):dense.geometry(128)
        with self.assertRaises(ValueError):dense.bracket(88,F(1,2),F(1,3),F(2))
        with self.assertRaises(ValueError):dense.bracket(88,F(1,3),F(1,2),F(1,2))


if __name__=='__main__':unittest.main()
