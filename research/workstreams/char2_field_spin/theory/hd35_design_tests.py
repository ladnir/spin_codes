"""Regression checks for explicitly parameterized distance screens."""
import unittest
from math import comb, log2

import numpy as np

from hd35_design_screen import FORBIDDEN_WEIGHT, LENGTH
from hd_sumrank_fast_screen import bound
from hd_sumrank_profile_pairs import exact_region_screen, screen
from hd_systematic_mds_screen import packet_logs


class DistanceScreenTests(unittest.TestCase):
    def test_regional_single_composition_and_marker_cancellation(self):
        from fractions import Fraction as F
        from flint import arb,ctx
        from hd_sumrank_paired_certificate import receipt
        options=dict(feedthrough='common2',outer_length=256,precision=256,
                     forbidden_weight=FORBIDDEN_WEIGHT,threshold=1000000)
        with ctx.workprec(256):
            single=receipt(1,2,F(999,1000),F(1),F(1,1000),paired=False,**options)
            paired=receipt(1,2,F(999,1000),F(1),F(1,1000),paired=True,**options)
            self.assertEqual(single['certified_intervals'],[[1,2]])
            self.assertTrue(arb(single['union_log2_upper']).overlaps(arb(paired['union_log2_upper'])))
            a=receipt(1,3,F(99,100),F(1,2),F(1,1000),paired=False,**options)
            b=receipt(1,3,F(99,100),F(1,2),F(1,256),paired=False,**options)
            self.assertEqual(a['certified_intervals'],[[1,3]])
            self.assertTrue(arb(a['union_log2_upper']).overlaps(arb(b['union_log2_upper'])))

    def test_all_small_map_return_constraints(self):
        """Every input support and output-zero set, not sampled minors."""
        from fractions import Fraction as F
        from hd_bucket_feedback import maps
        from hd_blockwise_shared_transfer import mm
        from hd_common_refresh_enumerator import extend
        from hd_common_refresh_packets import principal_defects
        q=65536
        for h in (1,2):
            length=4*h
            nu=principal_defects(h)
            actual_nu=[0]*(length+1)
            _,full=maps(local_rank=h)
            c=tuple(row[:length] for row in full[:h])
            ac=mm(tuple(zip(*c)),c,q)
            p=tuple(tuple(value^int(i==j) for j,value in enumerate(row))
                    for i,row in enumerate(ac))
            for support in range(1<<length):
                u=tuple(i for i in range(length) if support>>i&1)
                l=len(u)
                rr=()
                for row in c:rr=extend(rr,tuple(row[j] for j in u))
                self.assertEqual(len(rr),min(h,l))
                constraints=[()]
                quiet_constraints=[rr]
                for zero_set in range(1<<length):
                    if zero_set:
                        bit=zero_set&-zero_set
                        previous=zero_set^bit
                        row=tuple(p[bit.bit_length()-1][j] for j in u)
                        constraints.append(extend(constraints[previous],row))
                        quiet_constraints.append(extend(quiet_constraints[previous],row))
                    actual=F(1,q**len(constraints[zero_set]))-F(1,q**len(quiet_constraints[zero_set]))
                    active=(zero_set&support).bit_count()
                    inactive=zero_set.bit_count()-active
                    if inactive>=h:
                        self.assertEqual(actual,0)
                    else:
                        rank=max(0,active-min(nu[active],max(0,h-l+active)))
                        count=F(q**min(min(h,l),h-inactive)-1,q**(l-max(0,l-active-h)))
                        cap=min(1-F(1,q**min(h,l)),F(1,q**rank),count)
                        self.assertLessEqual(actual,cap)
                actual_nu[l]=max(actual_nu[l],l-len(constraints[support]))
            self.assertEqual(tuple(actual_nu),nu)

    def test_small_local_maps_and_exact_caps(self):
        from fractions import Fraction as F
        from hd_common_refresh_packets import local_return_cap, grouped_packets
        from hd_common_refresh_screen import local_return_cap as floating_cap
        from hd_bucket_feedback import maps
        from hd_blockwise_shared_transfer import mm
        for h in (1,2):
            _,c=maps(local_rank=h)
            a=tuple(zip(*c))
            self.assertFalse(any(value for row in mm(c,a,65536) for value in row))
            for l in range(4*h+1):
                exact=local_return_cap(l,F(3,4),h=h)
                approximate=floating_cap(l,.75,h,4*h)
                if exact:self.assertAlmostEqual(float(exact)/approximate,1.,places=11)
                else:self.assertEqual(approximate,0.)
            groups=grouped_packets(F(3,4),maximum=8,h=h)
            for j in range(9):
                self.assertEqual(sum(count for k,count,_ in groups if k==j),comb(64,j))

    def test_single_profile_composition_matches_directed(self):
        from fractions import Fraction as F
        from flint import arb, ctx
        from hd_common_refresh_packets import grouped_packets
        from hd_sumrank_global_certificate import witness
        from hd_randomized_outer_state_profile import power,total
        z,a,p=F(3,4),F(1,8),F(1,16)
        groups=grouped_packets(z,h=2)
        matrix=np.zeros((2,2))
        for j,count,entries in groups:
            matrix+=count*float(p)**j*float(1-p)**(64-j)*np.array([float(x)**float(a) for x in entries]).reshape((2,2))
        run=total(power(tuple(tuple(row) for row in np.log2(matrix)),LENGTH//64)[0])
        with ctx.workprec(256):
            slope,intercept=witness(z,a,p,feedthrough='common2',outer_length=256,forbidden_weight=FORBIDDEN_WEIGHT)
            expected=run-float(a)*(log2(65535)+FORBIDDEN_WEIGHT*log2(float(z)))-LENGTH*log2(float(1-p))
            self.assertAlmostEqual(float(intercept),expected,places=6)

    def test_exact_common_local_caps(self):
        from fractions import Fraction as F
        from hd_common_refresh_packets import local_return_cap, local_moments, matrix, grouped_packets
        from hd_common_refresh_screen import local_return_cap as floating_cap, matrices
        for z in (F(3,4), F(1,4)):
            for l in range(17):
                exact = local_return_cap(l, z)
                approximate = floating_cap(l, float(z), 4, 16)
                if exact:
                    self.assertAlmostEqual(float(exact)/approximate, 1., places=11)
                else:
                    self.assertEqual(approximate, 0.)
            local = local_moments(z)
            for occ in ((0,0,0,0), (0,0,1,1), (0,0,8,8), (4,4,4,4), (16,16,16,16)):
                self.assertTrue(all(v >= 0 for v in matrix(occ,z,local)))
            groups = grouped_packets(z)
            degrees, weight, logs = matrices(float(z))
            # Grouping can split at a floating rounding boundary. Compare
            # total transition mass by degree, not the incidental grouping.
            expected = np.zeros((65,4))
            actual = np.zeros((65,4))
            for j, count, entries in groups:
                expected[j] += count*np.array([float(v) for v in entries])
            for j, count, entries in zip(degrees, np.exp2(weight), np.exp2(logs).reshape((-1,4))):
                actual[j] += count*entries
            np.testing.assert_allclose(actual, expected, rtol=1./2048, atol=0)

    def test_integer_boundary(self):
        self.assertEqual(FORBIDDEN_WEIGHT, 734003)
        self.assertLessEqual(20*FORBIDDEN_WEIGHT, 7*LENGTH)
        self.assertGreater(20*(FORBIDDEN_WEIGHT+1), 7*LENGTH)

    def test_single_packet_exact_distance_shift(self):
        z, alpha = .75, .125
        packets = alpha*packet_logs(z, 16)
        old = bound(128, z, alpha, -6., packets)
        explicit = bound(128, z, alpha, -6., packets, forbidden_weight=838860)
        new = bound(128, z, alpha, -6., packets, forbidden_weight=FORBIDDEN_WEIGHT)
        self.assertEqual(old, explicit)
        self.assertAlmostEqual(new-old, alpha*(838860-FORBIDDEN_WEIGHT)*log2(z), places=8)

    def test_region_normalization_and_distance_shift(self):
        width, z, alpha, rows = 128, .75, .5, 3
        transition = np.array(((.6, .4), (.2, .8)))
        coeff = np.array([log2(comb(width, j))+np.log2(transition) for j in range(width+1)])
        for outer in (64, 128, 256):
            length = LENGTH//outer
            expected = log2(comb(length, rows))+alpha*(
                rows*(16*outer-(outer//2)*log2(65535))-log2(65535)
                -FORBIDDEN_WEIGHT*log2(z))
            actual = exact_region_screen(rows, z, alpha, coeff, width=width,
                outer_length=outer, forbidden_weight=FORBIDDEN_WEIGHT)
            self.assertAlmostEqual(actual['log2_upper_screen'], expected, places=6)
            old = screen(rows, z, alpha, coeff, width=width, outer_length=outer)
            explicit = screen(rows, z, alpha, coeff, width=width,
                              outer_length=outer, forbidden_weight=838860)
            new = screen(rows, z, alpha, coeff, width=width,
                         outer_length=outer, forbidden_weight=FORBIDDEN_WEIGHT)
            self.assertEqual(old, explicit)
            self.assertAlmostEqual(new['log2_upper_screen']-old['log2_upper_screen'],
                alpha*(838860-FORBIDDEN_WEIGHT)*log2(z), places=6)

    def test_invalid_geometry_and_distance(self):
        coeff = np.zeros((129, 2, 2))
        for forbidden in (-1, LENGTH, .35*LENGTH):
            with self.assertRaises(ValueError):
                screen(3, .75, .5, coeff, forbidden_weight=forbidden)
        with self.assertRaises(ValueError):
            screen(3, .75, .5, coeff, width=127)

    def test_directed_certificate_distance_shift(self):
        from fractions import Fraction as F
        from flint import arb, ctx
        from hd_sumrank_paired_certificate import receipt
        from hd_sumrank_global_certificate import witness
        z, a, p = F(999,1000), F(1,2), F(1,1000)
        with ctx.workprec(256):
            # Build the exact rational Arb mark, not a binary float.
            expected = arb(838860-FORBIDDEN_WEIGHT)*(arb(999)/1000).log()/arb(2).log()/2
            for outer in (128,256):
                options = dict(feedthrough='systematic_mds',outer_length=outer,precision=256)
                old = receipt(1,1,z,a,p,threshold=10000000,**options)
                new = receipt(1,1,z,a,p,threshold=10000000,forbidden_weight=FORBIDDEN_WEIGHT,**options)
                difference = arb(new['union_log2_upper'])-arb(old['union_log2_upper'])
                self.assertTrue(difference.overlaps(expected))
                old_slope, old_intercept = witness(z,a,p,**options)
                new_slope, new_intercept = witness(z,a,p,forbidden_weight=FORBIDDEN_WEIGHT,**options)
                self.assertTrue((new_intercept-old_intercept).overlaps(expected))
                self.assertTrue(new_slope.overlaps(old_slope))


if __name__ == '__main__':
    unittest.main()
