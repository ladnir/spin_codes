"""Exact averaging, original-provider compatibility, and directed swap controls."""
from fractions import Fraction as F
from itertools import product
from math import comb
import unittest
from unittest.mock import patch

from flint import arb,ctx

from hd_shifted_alias_certificate import rational
from hd_shifted_cover_backend import CachedBound,CONSTRUCTION
from hd_shifted_endpoint import endpoint_packet_caps,fractional_average,swap_endpoint_receipt
from hd_shifted_general_certificate import general_coefficients
from hd_shifted_swap_certificate import (average_coefficients,SwapCachedBound,SWAP_CONSTRUCTION,
                                        boundary_receipt)
from hd_shifted_unoriented_certificate import collapse_coefficients


class SwapCertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with ctx.workprec(256):
            cls.oriented = general_coefficients(F(1,32),F(1,4),(F(1),)*3,bits=60)
            cls.averaged = average_coefficients(cls.oriented)
            cls.maximum = collapse_coefficients(cls.oriented)

    def test_small_labelled_arrangements_match_independent_swap_average(self):
        width = 4
        coefficients = [[arb(0) for _ in range(width**3)] for _ in range(4)]
        unsigned = {}
        for slots in product(range(4),repeat=3):
            u,v,f = slots.count(1),slots.count(2),slots.count(3)
            cap = tuple(1+sum((i+1)*(state+1)**(entry+1)
                             for i,state in enumerate(slots)) for entry in range(4))
            index = u*width**2+v*width+f
            for entry in range(4):
                coefficients[entry][index] += cap[entry]
            pattern = tuple(1 if state==2 else state for state in slots)
            unsigned.setdefault(pattern,[]).append(cap)
        result = average_coefficients(coefficients,width=width)
        for (s,f),values in result.items():
            exact = [F(0) for _ in range(4)]
            for pattern,caps in unsigned.items():
                if pattern.count(1)==s and pattern.count(3)==f:
                    self.assertEqual(len(caps),2**s)
                    for entry in range(4):
                        exact[entry] += F(sum(row[entry] for row in caps),len(caps))
            for actual,expected in zip(values,exact):
                self.assertTrue(actual==rational(expected))

    def test_actual_averages_are_dominated_by_orientation_maximum(self):
        with ctx.workprec(256):
            for degree,values in self.averaged.items():
                for value,maximum in zip(values,self.maximum[degree]):
                    # Exact convex domination; allow only the few final Arb
                    # endpoint ulps introduced by summing different balls.
                    self.assertTrue(value<=maximum*(1+arb(2)**-240))
            for (s,f),values in self.averaged.items():
                if s<=1:
                    for value,maximum in zip(values,self.maximum[s,f]):
                        self.assertTrue(abs(value-maximum)<=maximum*arb(2)**-240)

    def test_endpoint_coefficient_matches_direct_fair_swaps(self):
        with ctx.workprec(256):
            expected = fractional_average(((F(n,2**32),row)
                                            for _,n,row in endpoint_packet_caps(F(1,32),60)),F(1,4))
            for value,reference in zip(self.averaged[32,0],expected):
                self.assertTrue(value>=reference.lower())
                self.assertTrue(value<=reference.upper()*(1+arb(2)**-230))

    def test_endpoint_only_provider_reproduces_full_chronological_receipt(self):
        coefficient = self.averaged[32,0]
        class EndpointOnly(SwapCachedBound):
            @classmethod
            def unmarked_coefficients(cls,z,alpha,precision,bits):
                return {(32,0):coefficient},{}
        bound = EndpointOnly(F(1,32),F(1,4),((1,1),(1,1)),precision=256)
        receipt = bound.receipt(((4096,4096),(0,0),(0,0)))
        control = swap_endpoint_receipt(precision=256)
        self.assertEqual(receipt['log2_upper_rational'],control['log2_upper_rational'])
        self.assertEqual(receipt['construction'],SWAP_CONSTRUCTION)
        self.assertTrue(bound.witness()['pair_swaps'])
        self.assertEqual(bound.witness()['orientation_aggregation'],'independent-fair-swaps')
        self.assertFalse(receipt['whole_code_certificate'])
        for uniform in (False,True):
            alias = EndpointOnly(F(1,32),F(1,4),((1,1),(1,1)),precision=256,
                                 grid=64,uniform_aliases=uniform,run_markers=(1,1))
            result = alias.receipt(((4096,4096),(0,0),(0,0)))
            self.assertEqual(result['construction'],SWAP_CONSTRUCTION)
            self.assertTrue(alias.witness()['pair_swaps'])
            self.assertEqual(alias.witness()['uniform_aliases'],uniform)
            self.assertLess(abs(F(result['log2_upper_rational'])
                                -F(control['log2_upper_rational'])),F(1,1000))

    def test_default_provider_remains_the_original_cached_function(self):
        token = object()
        with patch('hd_shifted_cover_backend.unmarked',return_value=token) as original:
            self.assertIs(CachedBound.unmarked_coefficients(F(1,2),F(1,4),256,60),token)
            original.assert_called_once_with(F(1,2),F(1,4),256,60)
        self.assertEqual(CachedBound.construction,CONSTRUCTION)
        self.assertNotEqual(SwapCachedBound.construction,CONSTRUCTION)

    def test_full_boundary_provider_at_two_precisions(self):
        first = boundary_receipt(precision=256)
        second = boundary_receipt(precision=384)
        self.assertEqual(first['integer_profile_triples'],47377)
        self.assertEqual(first['log2_upper_rational'],second['log2_upper_rational'])
        self.assertEqual(F(first['log2_upper_rational']),F(-11303619398642231,200000000000))
        self.assertEqual(first['construction'],SWAP_CONSTRUCTION)
        self.assertFalse(first['whole_code_certificate'])

    def test_general_sparse_profile_with_interval_aliases(self):
        markers = (F('1321103641/494495897869'),F('1358293346/797067919595'))
        bound = SwapCachedBound(F('4447531136429335/4503599627370496'),
                                F('7230271873918145/72057594037927936'),
                                (markers,markers),precision=256,grid=64)
        result = bound.receipt(((16,16),(4,4),(4,4)))
        self.assertLess(abs(F(result['log2_upper_rational'])-F('-144.174076029899')),F(1,10**6))
        self.assertEqual(result['support_count_method'],'joint-profile-count')
        self.assertFalse(result['whole_code_certificate'])

    def test_invalid_array_shapes(self):
        for width,arrays in ((0,[[]]*4),(3,[[arb(0)]*27]*3),(3,[[arb(0)]*26]*4)):
            with self.assertRaises(ValueError):
                average_coefficients(arrays,width=width)


if __name__=='__main__':
    unittest.main()
