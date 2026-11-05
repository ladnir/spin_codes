"""Small exact and directed replay tests; no full-code numerical gate."""

from fractions import Fraction as F
import unittest

from flint import arb,ctx

import hd_field128_lift_certificate as lift
from hd_common_refresh_packets import local_moments,matrix
from hd_grouped_quiet import quiet_cap
from hd_shifted_swap_sparse import sparse_geometry,weight_receipt
from hd_shifted_good_pair_screen import good_geometry
from hd_grouped_packet_screen import OCC


class DirectedLiftTests(unittest.TestCase):
    def test_exact_packet_formulas(self):
        for q,z in ((256,F(7,8)),(1<<128,F(7,8)),(1<<128,F(1,1<<200))):
            local = local_moments(z,q=q,h=4)
            density = F(q**16-1,(q-1)**16)
            for occ in ((0,0,0,0),(0,0,0,1),(1,2,3,4),(4,4,4,4),(16,16,16,16)):
                k,b,r,v = matrix(occ,z,local,q=q,h=4)
                live = density*v
                if not sum(occ):
                    live = min(live,quiet_cap(z,q=q,state=16,width=64,h=4))
                self.assertEqual(lift.packet(occ,z,q),(k,density*b,r,live))
                if sum(occ):
                    self.assertTrue(all(value>0 for value in lift.packet(occ,z,q)))

    def test_geometry_and_swap_counts(self):
        for maximum in (1,2,4):
            self.assertEqual(lift.sparse_geometry(maximum),sparse_geometry(maximum))
            s,f,n,index = good_geometry()
            expected = tuple((int(ss),int(ff),int(nn),tuple(map(int,OCC[ii])))
                             for ss,ff,nn,ii in zip(s,f,n,index) if ss+2*ff<=maximum)
            self.assertEqual(lift.good_geometry(maximum),expected)
            for ss in range(maximum+1):
                for ff in range((maximum-ss)//2+1):
                    count = sum(F(n,2**s) for s,f,n,_ in lift.sparse_geometry(maximum)
                                if s==ss and f==ff)
                    from math import comb
                    self.assertEqual(count,comb(32,ss)*comb(32-ss,ff))

    def test_upward_rounding_retains_tiny_values(self):
        value = F(7,(1<<2400)+1)
        rounded = lift.round_up(value,80)
        self.assertGreaterEqual(rounded,value)
        self.assertLess(rounded,value*(1+F(1,1<<78)))
        with ctx.workprec(256):
            self.assertTrue(lift.rational(rounded)>0)

    def test_carry_free_truncation(self):
        # A discarded singleton degree 4 cannot become full-pair degree 1.
        with ctx.workprec(256):
            maximum,budget = 2,lift.Budget(10)
            source = {(2,0):(arb(1),arb(0),arb(0),arb(1))}
            squared = lift.polynomial_product(lift.encode(source,maximum,F(1)),
                                              lift.encode(source,maximum,F(1)),maximum,budget)
            self.assertTrue(all(poly[index].is_zero() for poly in squared for index in range(6)))

    def test_positive_polynomial_product(self):
        with ctx.workprec(256):
            tiny = F(1,1<<2048)
            source = {(0,0):tuple(map(arb,(1,0,0,1))),
                      (1,0):tuple(lift.rational(x) for x in (F(1),tiny,tiny*tiny,F(1,1<<512)))}
            encoded = lift.encode(source,2,F(1))
            squared = lift.polynomial_product(encoded,encoded,2,lift.Budget(10))
            self.assertTrue(squared[1][1]>=lift.rational(2*tiny))
            self.assertTrue(squared[2][1]>=lift.rational(2*tiny*tiny))
            self.assertTrue(squared[1][1]>0)

    def test_q256_tiny_directed_replay(self):
        old = weight_receipt(1,F(999,1000),F(1),precision=256,bits=80)
        new = lift.weight_receipt(1,F(999,1000),F(1),q=256,precision=256,bits=80,
                                  outer_measure='inherited-diagonal',max_seconds=20)
        self.assertLessEqual(abs(F(old['log2_upper_rational'])-F(new['log2_upper_rational'])),F(1,10**10))
        self.assertEqual(new['profiles'],old['profiles'])
        self.assertTrue(new['directed'])
        self.assertFalse(new['whole_code_certificate'])
        self.assertEqual(new['projective_divisor'],'255')

    def test_outer_density_and_source_closure(self):
        for q in (256,1<<128):
            self.assertEqual(lift.outer_density(q),F(q**160,(q-1)**32))
            self.assertLessEqual(lift.outer_density(q),lift.outer_density(q,'inherited-diagonal'))
        pins = lift.source_pins()
        for name in ('hd_field128_lift_certificate.py','hd_common_refresh_principal.py',
                     'hd_bucket_feedback.py','hd_systematic_mds.py','hd_blockwise_shared_transfer.py',
                     'hd_grouped_quiet.py','hd_shifted_pairs_count.py'):
            self.assertIn(name,pins)
        self.assertNotIn('hd_field128_lift_screen.py',pins)
        self.assertTrue(all(len(digest)==64 for digest in pins.values()))

    def test_scope_validation_and_budget(self):
        for value in (.9,.1):
            with self.assertRaises(ValueError):
                lift.exact(value)
        for q in (16,32,1024,1<<136):
            with self.assertRaises(ValueError):
                lift.checked_q(q)
        with self.assertRaises(ValueError):
            lift.regional_directed(F(7,8),F(1,16),65)
        with self.assertRaises(TimeoutError):
            lift.Budget(1e-12).check()


if __name__=='__main__':
    unittest.main()
