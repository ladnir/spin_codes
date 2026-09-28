"""Exact trajectory tests for persistent-type Holder envelopes."""
from fractions import Fraction as Q
from itertools import product
from math import prod
import unittest

from flint import ctx

from holder_exact import coordinate_sums,dyadic_root_bound,holder_envelope


class HolderTests(unittest.TestCase):
    def test_all_tiny_trajectories(self):
        families=[
            [(Q(2,3),(Q(1,2),Q(1,2),Q(0),Q(0),Q(0))),
             (Q(7,4),(Q(0),Q(1,3),Q(1,3),Q(0),Q(1,3)))],
            [(Q(5),(Q(1,5),)*5),
             (Q(2),(Q(1,7),Q(2,7),Q(3,7),Q(1,7),Q(0))),
             (Q(1,11),(Q(0),Q(0),Q(1),Q(0),Q(0)))],
            [(Q(1),(Q(1),Q(0),Q(0),Q(0),Q(0))),
             (Q(3),(Q(0),Q(0),Q(0),Q(0),Q(1)))]
        ]
        for family in families:
            for m in (2,3):
                z,phi=holder_envelope(family,m=m,precision=96)
                self.assertEqual(sum(phi),1)
                for trajectory in product(range(5),repeat=m):
                    actual=sum(h*prod(theta[b] for b in trajectory) for h,theta in family)
                    upper=z**m*prod(phi[b] for b in trajectory)
                    self.assertLessEqual(actual,upper)
                for b in range(5):
                    constant=sum(h*theta[b]**m for h,theta in family)
                    self.assertLessEqual(constant,(z*phi[b])**m)

    def test_dyadic_roots_have_exact_minimal_grid_certificates(self):
        values=[Q(0),Q(1),Q(2),Q(4,9),Q(123,789),Q(1<<5000,7),Q(7,1<<5000)]
        for value in values:
            for m in (1,2,3,256):
                upper,quantum=dyadic_root_bound(value,m,precision=96)
                self.assertGreaterEqual(upper**m,value)
                if value:
                    self.assertGreater(quantum,0)
                    self.assertLess((upper-quantum)**m,value)
                    self.assertEqual((upper/quantum).denominator,1)
                    self.assertEqual(upper.denominator&(upper.denominator-1),0)
                    self.assertEqual(quantum.numerator&(quantum.numerator-1),0)
                    self.assertEqual(quantum.denominator&(quantum.denominator-1),0)
                else:
                    self.assertEqual((upper,quantum),(0,0))

    def test_zero_and_degenerate_families(self):
        expected=(Q(0),(Q(1),Q(0),Q(0),Q(0),Q(0)))
        self.assertEqual(holder_envelope([],m=3),expected)
        self.assertEqual(holder_envelope([(Q(0),(Q(1,5),)*5)],m=3),expected)
        for b in range(5):
            point=tuple(Q(int(i==b)) for i in range(5))
            z,phi=holder_envelope([(Q(9),point)],m=2)
            self.assertEqual(z,3)
            self.assertEqual(phi,point)
        theta=(Q(1,2),Q(1,4),Q(0),Q(1,8),Q(1,8))
        z,phi=holder_envelope([(Q(1),theta)],m=256)
        self.assertEqual(z,1)
        self.assertEqual(phi,theta)

    def test_coordinate_sums_are_exact_and_type_splitting_invariant(self):
        theta=(Q(1,2),Q(1,3),Q(1,6),Q(0),Q(0))
        merged=[(Q(5,7),theta)]
        split=[(Q(2,7),theta),(Q(3,7),theta)]
        self.assertEqual(coordinate_sums(merged,m=3),tuple(Q(5,7)*p**3 for p in theta))
        self.assertEqual(holder_envelope(merged,m=3),holder_envelope(split,m=3))

    def test_precision_restored_and_no_float_input(self):
        previous=ctx.prec
        try:
            ctx.prec=64
            dyadic_root_bound(Q(2),3,precision=192)
            self.assertEqual(ctx.prec,64)
        finally:
            ctx.prec=previous
        with self.assertRaises(TypeError):
            holder_envelope([(1.,(Q(1),Q(0),Q(0),Q(0),Q(0)))])
        with self.assertRaises(TypeError):
            holder_envelope([(Q(1),(1.,Q(0),Q(0),Q(0),Q(0)))])
        with self.assertRaises(TypeError):
            dyadic_root_bound(2.,2)

    def test_invalid_parameters(self):
        for m in (0,-1,Q(3,2)):
            with self.assertRaises(ValueError):
                holder_envelope([],m=m)
        for family in ([(Q(-1),(Q(1),Q(0),Q(0),Q(0),Q(0)))],
                       [(Q(1),(Q(1),Q(0)))],
                       [(Q(1),(Q(1),Q(1),Q(0),Q(0),Q(0)))],
                       [(Q(1),(Q(2),Q(-1),Q(0),Q(0),Q(0)))]):
            with self.assertRaises(ValueError):
                holder_envelope(family)
        with self.assertRaises(ValueError):
            dyadic_root_bound(Q(-1),2)
        with self.assertRaises(ValueError):
            holder_envelope([],precision=8)


if __name__=='__main__':
    unittest.main()
