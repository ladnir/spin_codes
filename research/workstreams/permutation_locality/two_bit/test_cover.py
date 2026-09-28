from fractions import Fraction as Q
from math import comb
from unittest.mock import patch
import unittest

import numpy as np
from flint import arb,arb_mat,ctx

import model
from cover import complete_cover,outward_box,DENOMINATOR,check_partition,replay_cover
from shell_cover import IntervalFolds
from occupancy_adaptive import geometry_test,split
from full_cover import saved_threshold,ENSEMBLE_PREFIX,replay_saved_covers,validate_record


class CoverTests(unittest.TestCase):
    def setUp(self):ctx.prec=192

    def test_reused_partition_geometry(self):
        geometry_test()

    def test_saved_partition_rejects_gaps_overlaps_and_bad_labels(self):
        root=((38,256),)*2
        children=[dict(box=box,multiplicity=mult) for box,mult in split(root,1)]
        check_partition(children,2)
        with self.assertRaises(ValueError):check_partition(children[:-1],2)
        with self.assertRaises(ValueError):check_partition(children+[dict(box=root,multiplicity=1)],2)
        corrupted=[dict(item) for item in children];corrupted[0]['multiplicity']+=1
        with self.assertRaises(ValueError):check_partition(corrupted,2)

    def test_location_factor_and_fold(self):
        shells=[Q(int(u==128)) for u in range(257)]
        cdf=[Q(int(u>=128)) for u in range(257)]
        fold=IntervalFolds(cdf,shells)
        matrix=arb_mat([[int(i==j) for j in range(9)] for i in range(9)])
        args=([matrix]*3,fold,'0',((128,128),)*2,1,[DENOMINATOR//2]*2)
        current=outward_box(*args);old=outward_box(*args,groups=2048)
        expected=arb(comb(4096,2))*arb(2)**512/arb(comb(256,128))**2
        self.assertTrue(current>=expected.lower())
        self.assertTrue(abs(current/old-arb(comb(4096,2))/comb(2048,2))<arb(2)**-150)
        with self.assertRaises(ValueError):outward_box([matrix],*args[1:])

    def test_complete_root_replay(self):
        shells=[Q(int(u==128)) for u in range(257)]
        cdf=[Q(int(u>=128)) for u in range(257)]
        matrix=arb_mat([[arb(int(i==j))/16 for j in range(9)] for i in range(9)])
        floating=np.array([[[float(matrix[i,j]) for j in range(9)] for i in range(9)]]*3)
        certificate={};operators={'0':([matrix]*3,floating)}
        with patch('builtins.print'):
            bound=complete_cover(operators,cdf,shells,2,target_bits=55,max_splits=0,certificate=certificate)
            screen=complete_cover({'0':([matrix]*3,floating)},cdf,shells,2,target_bits=55,max_splits=0,screen_only=True)
        self.assertTrue(0<bound<arb(2)**-55);self.assertIsNone(screen)
        self.assertEqual(bound,replay_cover(operators,cdf,shells,2,certificate['leaves']))

    def test_output_cutoff_changes_the_chernoff_factor_only(self):
        shells=[Q(int(u==128)) for u in range(257)]
        cdf=[Q(int(u>=128)) for u in range(257)]
        fold=IntervalFolds(cdf,shells)
        matrix=arb_mat([[int(i==j) for j in range(9)] for i in range(9)])
        args=([matrix]*3,fold,'1/100',((128,128),)*2,1,[DENOMINATOR//2]*2)
        low=outward_box(*args,output_cutoff=100)
        high=outward_box(*args,output_cutoff=200)
        self.assertTrue(abs(high/low-arb(1).exp())<arb(2)**-150)
        with self.assertRaises(ValueError):outward_box(*args,output_cutoff=-1)

    def test_saved_sparse_cutoff_is_explicit_and_legacy_is_ten_percent(self):
        record=dict(schema='two-bit-support-cover-1',ensemble=[*ENSEMBLE_PREFIX,2])
        self.assertEqual(saved_threshold(record),209715)
        self.assertEqual(saved_threshold({**record,'threshold':199229}),199229)
        with self.assertRaises(ValueError):saved_threshold({**record,'threshold':2097152})

    def test_complete_cover_serializes_and_replays_changed_cutoff(self):
        shells=[Q(int(u==128)) for u in range(257)]
        cdf=[Q(int(u>=128)) for u in range(257)]
        matrix=arb_mat([[arb(int(i==j))/16 for j in range(9)] for i in range(9)])
        floating=np.array([[[float(matrix[i,j]) for j in range(9)] for i in range(9)]]*3)
        certificate={};operators={'.001':([matrix]*3,floating)}
        with patch('builtins.print'):
            bound=complete_cover(operators,cdf,shells,2,max_splits=0,certificate=certificate,output_cutoff=100)
        self.assertEqual(certificate['threshold'],100)
        self.assertEqual(bound,replay_cover(operators,cdf,shells,2,certificate['leaves'],output_cutoff=100))
        self.assertTrue(replay_cover(operators,cdf,shells,2,certificate['leaves'],output_cutoff=200)>bound)

    def test_sparse_batch_reconstructs_bounds_instead_of_trusting_saved_values(self):
        def record(q):
            return dict(schema='two-bit-support-cover-1',ensemble=[*ENSEMBLE_PREFIX,2],q=q,
                        threshold=199229,penalty='1',cutoff=2,output_degree=64,density_fold=True,
                        upper_dyadic=[1,100000],
                        leaves=[dict(box=[[38,256]]*q,multiplicity=1,tilt='.1',numerators=[DENOMINATOR//2]*q)])
        records=[record(2),record(3)]
        def replay(operators,cdf,shells,q,leaves,*,density_fold,output_cutoff):
            self.assertEqual(output_cutoff,199229);self.assertTrue(density_fold)
            return arb(2)**(-100*q)
        with patch('full_cover.build_operators',return_value={}) as build, \
             patch('full_cover.weighted_cdf_upper',return_value=[]), \
             patch('full_cover.weighted_union_shells',return_value=[1]), \
             patch('full_cover.replay_cover',side_effect=replay),patch('builtins.print'):
            results=replay_saved_covers(records,None)
        self.assertEqual(results,{2:arb(2)**-200,3:arb(2)**-300})
        build.assert_called_once_with(['.1'],3,2,'1',64,2,floating=False)
        with self.assertRaises(ValueError):replay_saved_covers([records[0],records[0]],None)
        with self.assertRaises(ValueError):validate_record({**records[0],'cutoff':0})
        with self.assertRaises(ValueError):validate_record({**records[0],'output_degree':65})

    def test_unpassed_cover_never_certifies(self):
        shells=[Q(int(u==128)) for u in range(257)]
        cdf=[Q(int(u>=128)) for u in range(257)]
        matrix=arb_mat([[int(i==j) for j in range(9)] for i in range(9)])
        floating=np.array([np.eye(9)]*3)
        with patch('builtins.print'):
            value=complete_cover({'0':([matrix]*3,floating)},cdf,shells,2,target_bits=55,max_splits=0)
        self.assertIsNone(value)


if __name__=='__main__':unittest.main()
