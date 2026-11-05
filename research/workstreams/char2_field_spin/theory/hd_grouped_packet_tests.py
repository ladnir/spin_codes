"""Exact combinatorial and local-operator tests for the packet-route screen."""
from fractions import Fraction as F
from itertools import product
from math import comb, factorial
import unittest

import numpy as np

from hd_common_refresh_packets import local_moments, matrix
from hd_grouped_packet_screen import OCC, bound, coefficients, local_logs, return_bernstein, size_profiles
from hd_grouped_quiet import quiet_cap


class PacketRouteTests(unittest.TestCase):
    def test_gale_ryser_matches_all_four_position_mask_sequences(self):
        actual = {}
        for masks in product(range(16),repeat=4):
            counts = tuple(sum(mask.bit_count()==s for mask in masks) for s in range(1,5))
            occupancies = tuple(sorted(sum((mask>>lane)&1 for mask in masks) for lane in range(4)))
            actual.setdefault(counts,set()).add(occupancies)
        for counts,_,indices in size_profiles((1,2,3,4)):
            if sum(counts)<=4:
                self.assertEqual({tuple(OCC[i]) for i in indices},actual[counts])

    def test_homogeneous_size_four_has_no_envelope_loss(self):
        for (count,),multiplicity,indices in size_profiles((4,)):
            self.assertEqual(multiplicity,comb(16,count))
            self.assertEqual([tuple(OCC[i]) for i in indices],[(count,)*4])

    def test_size_profile_multiplicity(self):
        for sizes in ((1,),(4,),(1,4),(1,2,3,4)):
            self.assertEqual(sum(n for _,n,_ in size_profiles(sizes)),(len(sizes)+1)**16)

    def test_float_local_entries_match_exact_majorants(self):
        z = F(3,4)
        exact_local = local_moments(z,q=256,h=4)
        beta = F(256**16-1,255**16)
        logs = local_logs(float(z))
        samples = ((0,0,0,0),(0,0,0,1),(1,1,1,1),(0,0,3,16),
                   (2,6,9,12),(8,8,8,8),(16,16,16,16))
        for occupancy in samples:
            k,birth,ret,live = matrix(occupancy,z,exact_local,q=256,h=4)
            expected = [k,beta*birth,ret,beta*live]
            if not any(occupancy):
                expected[3] = min(expected[3],quiet_cap(z,q=256,h=4))
            index = np.flatnonzero(np.all(OCC==occupancy,axis=1))[0]
            np.testing.assert_allclose(np.exp2(logs[index]).ravel(),
                                       np.array([float(v) for v in expected]),rtol=2e-12,atol=0)

    def test_one_group_matches_scalar_route_bound(self):
        from hd_grouped_subfield_screen import bound as scalar_bound,coefficients
        z,alpha = .999, .75
        old,_ = scalar_bound(1,z,alpha,coefficients(z,alpha,h=4,quiet_inverse=True),regional=True)
        new = bound((1,0,0,0),z,alpha,np.array([-11.]),regional=True)
        # The established FLOAT source rounds its local entries to 13 bits.
        self.assertLessEqual(new,old+1e-6)
        self.assertLess(abs(new-old),.1)

    def test_regional_markers_cancel(self):
        a = bound((0,0,0,4),.985,.5,np.array([-8.]),regional=True)
        b = bound((0,0,0,4),.985,.5,np.array([-6.]),regional=True)
        self.assertAlmostEqual(a,b,places=7)

    def test_paired_composition_is_equal_at_alpha_one(self):
        a = bound((0,0,0,4),.985,1.,np.array([-8.]),regional=True)
        b = bound((0,0,0,4),.985,1.,np.array([-8.]),regional=True,steps=2)
        self.assertAlmostEqual(a,b,places=7)

    def test_paired_composition_improves_fractional_bound(self):
        a = bound((0,0,0,16),.9566,.07,np.array([-7.]),regional=True)
        b = bound((0,0,0,16),.9566,.07,np.array([-7.]),regional=True,steps=2)
        self.assertLess(b,a)

    def test_larger_rank_empty_packet_uses_its_own_quiet_bound(self):
        for h in (8,16):
            z = F(9,10)
            occ = np.repeat((np.arange(17)*(h//4))[:,None],16//h,axis=1)
            logs = local_logs(.9,h,occ)
            flat = sum(F(comb(4*h,c),256**min(h,c))*z**(4*h-c)*(1-z)**c
                       for c in range(4*h+1))
            uniform = F(256**16-1,255**16)*(256**16*flat**(16//h)-1)/(256**16-1)
            expected = min(uniform,quiet_cap(z,q=256,h=h))
            self.assertAlmostEqual(2**logs[0,1,1],float(expected),places=13)
            with self.assertRaises(ValueError):
                coefficients(.9,.5,(1,),h)

    def test_two_group_route_pair_composition_at_alpha_one(self):
        for h in (4,8,16):
            a = bound((0,4,0,0),.985,1.,np.array([-9.]),regional=True,h=h,route_width=2)
            b = bound((0,4,0,0),.985,1.,np.array([-9.]),regional=True,h=h,route_width=2,steps=2)
            self.assertAlmostEqual(a,b,places=7)

    def test_two_group_route_total_combinatorial_mass(self):
        halves = tuple(product(range(17),repeat=2))
        for degree in range(33):
            self.assertEqual(sum(comb(16,a)*comb(16,b) for a,b in halves if a+b==degree),
                             comb(32,degree))

    def test_cached_return_polynomial_matches_original_helper(self):
        from hd_common_refresh_screen import local_return_cap
        for h in (4,8,16):
            length = 4*h
            degrees = np.arange(length+1)
            for z in (.49,.91):
                actual = return_bernstein(h)@(z**(length-degrees)*(1-z)**degrees)
                expected = [local_return_cap(l,z,h,length,256.) for l in degrees]
                np.testing.assert_allclose(actual,expected,rtol=5e-14,atol=0)

    def test_mixed_pair_global_composition_at_alpha_one(self):
        for steps in (1,2):
            value = bound((2,3,0,0),.99,1.,np.array([-8.,-7.]),h=16,route_width=2,steps=steps)
            if steps==1:
                first = value
            else:
                self.assertAlmostEqual(first,value,places=7)

    def test_exact_four_composition_at_alpha_one(self):
        a = bound((0,4,0,0),.985,1.,np.array([-9.]),regional=True,h=16,route_width=2)
        b = bound((0,4,0,0),.985,1.,np.array([-9.]),regional=True,h=16,route_width=2,steps=4,exact_four=True)
        self.assertAlmostEqual(a,b,places=7)

    def test_exact_four_retains_route_fluctuations(self):
        options = dict(regional=True,h=16,route_width=2)
        histogram,z,alpha,odds = (0,512,0,0),.48162188608540357,.08341661300705136,np.array([-2.807354922057604])
        paired = bound(histogram,z,alpha,odds,steps=2,**options)
        exact = bound(histogram,z,alpha,odds,steps=4,exact_four=True,**options)
        averaged = bound(histogram,z,alpha,odds,steps=4,**options)
        self.assertLess(exact,paired)
        self.assertLess(exact,averaged)
        self.assertLess(exact,-800)


if __name__=='__main__':
    unittest.main()
