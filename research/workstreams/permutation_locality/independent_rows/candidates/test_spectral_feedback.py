"""Exact small-state tests of compressed Fourier sums and atom bounds."""
import unittest
import numpy as np

from spectral_feedback import summaries, build, compare_known, combined_records
from pair_tail import hadamard


class SpectralFeedback(unittest.TestCase):
    def check_counts(self,counts,expansion):
        transforms=np.array([hadamard(np.array(row,dtype=np.int64)) for row in counts])
        patterns,inverse,multiplicities=np.unique(transforms.T,axis=0,return_inverse=True,return_counts=True)
        result=summaries(patterns.T,[sum(row) for row in counts],multiplicities,inverse,np.array(expansion,dtype=np.int64))
        for actual,(zero,peak,den,classes) in zip(counts,result):
            self.assertEqual(zero,actual[0]); self.assertEqual(den,sum(actual))
            self.assertGreaterEqual(peak,max(actual)); self.assertLessEqual(peak,den)
            self.assertEqual(classes,{v:sum(n for n,w in zip(actual,expansion) if w==v) for v in set(expansion[1:])})

    def test_complete_small_distributions(self):
        from itertools import product
        rows=[x for x in product(range(3),repeat=4) if sum(x)>0]
        self.check_counts(rows,(0,2,2,4))
        self.check_counts([(1,2,0,3,1,0,4,2),(0,)*7+(17,),(7,)*8],(0,1,2,2,3,3,3,4))

    def test_large_exact_dot_products(self):
        # Fixed-width intermediate products would exceed 2^63 here.
        self.check_counts([tuple((1<<57)*x for x in (7,1,3,0,2,4,0,1))],(0,1,1,1,2,2,2,2))

    def test_rejects_inconsistent_inputs(self):
        with self.assertRaises(ValueError):
            summaries(np.array([[2,0]]),[2],np.array([1,3]),np.array([0,0,1,1]),np.array([0,1,1,2]))
        with self.assertRaises(ValueError):
            summaries(np.array([[2.,0.]]),[2],np.array([1,3]),np.array([0,1,1,1]),np.array([0,1,1,2]))
        for maximum in (0,11,True):
            with self.assertRaises(ValueError): build(maximum)

    def test_reference_census_may_include_empty_shape(self):
        known={(): (1,0,1,{1:0}),(1,):(0,2,4,{1:4})}
        result={(1,):(0,6,8,{1:8})}
        self.assertEqual(compare_known(result,known,1),1)
        result[(1,)]=(0,1,8,{1:8})
        with self.assertRaises(AssertionError): compare_known(result,known,1)

    def test_peak_intersection_preserves_exact_class_counts(self):
        from fractions import Fraction as Q
        from unittest.mock import patch
        known={(1,):(0,2,4,{1:4})}
        spectral={(1,):(0,3,4,{1:4}),(1,1):(1,8,17,{1:16})}
        with patch('spectral_feedback.conditioned_peaks',return_value={(1,):Q(1,2),(1,1):Q(1,3)}):
            result=combined_records(spectral,known,1)
        self.assertEqual(result[(1,)],known[(1,)])
        self.assertEqual(result[(1,1)],(1,6,17,{1:16}))
        self.assertEqual(spectral[(1,1)][1],8)


if __name__=='__main__':
    unittest.main()
