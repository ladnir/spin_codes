"""Bounded exact Walsh-width and direct density-envelope regressions."""
from fractions import Fraction as Q
from itertools import product
import unittest

import numpy as np
from flint import arb,ctx

from feedback_density import (INT64_MAX, class_convolution, class_transforms,
                              dyadic_powers, dyadic_precision, hadamard,
                              recover_counts, shape_bounds)


def rational(point):
    value = point.fmpq()
    return Q(int(value.numerator),int(value.denominator))


def traced_hadamard(values):
    """Python-integer reference records every partial butterfly magnitude."""
    current = list(map(int,values))
    peak = max(map(abs,current),default=0)
    width = 1
    while width < len(current):
        for start in range(0,len(current),2*width):
            for offset in range(width):
                left,right = current[start+offset],current[start+width+offset]
                current[start+offset],current[start+width+offset] = left+right,left-right
                peak = max(peak,abs(left+right),abs(left-right))
        width *= 2
    return current,peak


class FeedbackDensityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        ctx.prec = 192

    def test_all_small_partial_butterflies_obey_ND(self):
        for size in (2,4):
            for counts in product(range(3),repeat=size):
                denominator = sum(counts)
                if not denominator:
                    continue
                transformed,_ = traced_hadamard(counts)
                for indicator in product((0,1),repeat=size):
                    transformed_indicator,_ = traced_hadamard(indicator)
                    coefficients = [a*b for a,b in zip(transformed,transformed_indicator)]
                    convolution,peak = traced_hadamard(coefficients)
                    self.assertLessEqual(peak,size*denominator)
                    expected = [sum(counts[t ^ s]*indicator[s] for s in range(size)) for t in range(size)]
                    self.assertEqual(convolution,[size*x for x in expected])

    def test_exact_convolution_and_source_zero_exclusion(self):
        expansion = np.array([0,2,2,4,3,1,3,1],dtype=np.int64)
        prepared = class_transforms(expansion)
        for values in ((0,1,2,3,4,3,2,1),(17,0,0,0,0,0,0,0),(0,0,0,0,7,0,0,0)):
            counts = np.array(values,dtype=np.int64)
            denominator = int(counts.sum())
            transformed = hadamard(counts)
            np.testing.assert_array_equal(recover_counts(transformed,denominator),counts)
            total = np.zeros(len(counts),dtype=np.int64)
            for v,indicator in prepared[2].items():
                hits = class_convolution(transformed,counts,expansion,indicator,denominator,v)
                expected = [sum(int(counts[t ^ s]) for s in range(len(counts)) if expansion[s] == v)
                            for t in range(len(counts))]
                np.testing.assert_array_equal(hits,expected)
                total += hits
            np.testing.assert_array_equal(total,denominator-counts)

    def test_dyadic_powers_round_up(self):
        levels = (1,2,48,56,64,72,80)
        for bits in (5,23,44):
            for weight in (1,8,48):
                values = dyadic_powers(levels,weight,'.052',bits)
                ctx.prec = 384
                try:
                    for v in levels:
                        exact = (-arb(13)*abs(v-weight)/250).exp()
                        upper = arb(values[v])/(1 << bits)
                        self.assertTrue(upper >= exact or (upper-exact).contains(0))
                        self.assertLessEqual(values[v],1 << bits)
                finally:
                    ctx.prec = 192

    def test_direct_toy_output_density(self):
        # Two-bit state, four output bits, and one uniformly selected bit.
        # Actual output weights retain their correlation with feedback.
        images = (0,3,12,15)
        columns = (1,2,3,1)
        expansion = np.array([x.bit_count() for x in images],dtype=np.int64)
        counts = np.bincount(columns,minlength=4).astype(np.int64)
        prepared = class_transforms(expansion)
        result = shape_bounds(hadamard(counts),4,(1,),expansion,prepared,['.17'],maximum_bits=23)['.17']
        ctx.prec = 384
        try:
            for target in range(1,4):
                value = arb(0)
                classes = {v:arb(0) for v in prepared[0]}
                for state in range(1,4):
                    for bit,feedback in enumerate(columns):
                        if state ^ feedback != target:
                            continue
                        moment = (-arb(17)*(images[state] ^ (1 << bit)).bit_count()/100).exp()/16
                        value += moment
                        classes[expansion[state]] += moment/prepared[1][expansion[state]]
                self.assertTrue(value <= result['density'] or (value-result['density']).contains(0))
                for v,bound in result['uniform'].items():
                    self.assertTrue(classes[v] <= bound or (classes[v]-bound).contains(0))
        finally:
            ctx.prec = 192

    def test_integer_final_maximum_and_uniform_normalization(self):
        expansion = np.array([0,1,2,3,3,2,1,4],dtype=np.int64)
        counts = np.array([1,2,3,4,5,6,7,8],dtype=np.int64)
        denominator = int(counts.sum())
        prepared = class_transforms(expansion)
        result = shape_bounds(hadamard(counts),denominator,(1,2),expansion,prepared,['.03','.2'],maximum_bits=12)
        for tilt,record in result.items():
            powers = dyadic_powers(prepared[0],3,tilt,12)
            scalar_denominator = denominator*(1 << 12)*4
            exact = max(sum(int(counts[t ^ s])*powers[int(expansion[s])]
                            for s in range(1,len(counts))) for t in range(1,len(counts)))
            self.assertGreaterEqual(rational(record['density']),Q(exact,scalar_denominator))
            for v in prepared[0]:
                peak = max(sum(int(counts[t ^ s]) for s in range(len(counts)) if expansion[s] == v)
                           for t in range(1,len(counts)))
                expected = Q(peak*powers[v],scalar_denominator*prepared[1][v])
                self.assertGreaterEqual(rational(record['uniform'][v]),expected)

    def test_width_guards_and_largest_six_window_denominator(self):
        denominator = 751_631_892_480
        self.assertLess((1 << 19)*denominator,1 << 63)
        self.assertEqual(dyadic_precision(denominator),23)
        self.assertLessEqual(denominator*(1 << 23),INT64_MAX)
        self.assertGreater(denominator*(1 << 24),INT64_MAX)
        with self.assertRaises(ValueError):
            recover_counts(np.ones(8,dtype=np.int64),INT64_MAX//8+1)
        with self.assertRaises(ValueError):
            recover_counts(np.array([1,0,0,0],dtype=np.int64),1)
        with self.assertRaises(ValueError):
            dyadic_precision(INT64_MAX)
        with self.assertRaises(ValueError):
            dyadic_powers((1,2),1,'0',23)

    def test_nonintegral_or_overflowing_inputs_are_rejected(self):
        for values in ([1.,1.],np.array([1<<63,0],dtype=np.uint64),[1<<80,0]):
            with self.assertRaises(ValueError):
                recover_counts(values,1)
        with self.assertRaises(ValueError):
            class_transforms([0.,1.5])
        with self.assertRaises(ValueError):
            class_transforms([0,0,1,1])


if __name__ == '__main__':
    unittest.main()
