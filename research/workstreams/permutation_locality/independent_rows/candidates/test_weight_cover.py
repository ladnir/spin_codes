"""Guard the coupling of an inner weight witness and its outer measure."""
from fractions import Fraction as Q
from contextlib import redirect_stdout
from io import StringIO
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np
from flint import arb, arb_mat, ctx

import weight_cover as controller


class WeightCover(unittest.TestCase):
    def test_canonical_pair_keys_are_distinct(self):
        self.assertEqual(controller.measure_key('.9','.98'),controller.measure_key('9/10','49/50'))
        self.assertNotEqual(controller.measure_key('.9','1'),controller.measure_key('.9','.98'))
        self.assertNotEqual(controller.measure_key('.9','.98'),controller.measure_key('.98','.9'))
        for rho,a in (('0','1'),('1.1','1'),('1','0')):
            with self.assertRaises(ValueError): controller.measure_key(rho,a)

    def test_weighted_operators_and_counts_remain_paired(self):
        operators,counts,shells={},{},{}
        args=SimpleNamespace(weight_tilt='.98')
        built={('.056','.9'):'weighted operator'}
        with patch.object(controller,'build_operators',return_value=built) as builder, \
             patch.object(controller,'tilted_support_caps',return_value=([Q(0),Q(3,2)],[Q(0),Q(3,2)])) as outer:
            controller.append_grid(operators,counts,shells,[1],args)
            self.assertIsNot(builder.call_args.args[0],args)
            self.assertEqual(builder.call_args.args[0].weight_tilt,'.98')
            self.assertEqual(outer.call_args.kwargs,dict(full_weight=Q(10,9),input_weight=Q(49,50)))
        key=controller.measure_key('.9','.98')
        self.assertEqual(operators,{('.056',key):'weighted operator'})
        self.assertEqual(counts,{key:[0,2]})
        self.assertEqual(shells,{key:[0,Q(3,2)]})
        with patch.object(controller,'build_operators',return_value={('.056','.9'):'baseline operator'}), \
             patch.object(controller,'weighted_cdf_upper',return_value=[Q(0),Q(5,2)]), \
             patch.object(controller,'weighted_union_shells',return_value=[Q(1),Q(5,2)]):
            controller.append_grid(operators,counts,shells,[1],SimpleNamespace(weight_tilt='1'))
        baseline_key=controller.measure_key('.9','1')
        self.assertEqual(len(operators),2)
        self.assertEqual(counts[baseline_key],[0,3])
        self.assertEqual(shells[baseline_key],[0,Q(5,2)])
        self.assertEqual(counts[key],[0,2])
        self.assertEqual(operators[('.056',key)],'weighted operator')

    def test_duplicate_witness_is_rejected(self):
        key=controller.measure_key('.9','.98')
        with patch.object(controller,'build_operators',return_value={('.056','.9'):'new'}):
            with self.assertRaises(ValueError):
                controller.append_grid({('.056',key):'old'},{key:[0]},{key:[0]},[1],SimpleNamespace(weight_tilt='.98'))

    def test_opaque_measure_keys_reach_outward_cover(self):
        old_precision=ctx.prec
        try:
            ctx.prec=192
            labels=[controller.measure_key('.9',a) for a in ('1','.98')]
            operators={('.001',key):([arb_mat([[arb(2)**-16]])]*2,[np.array([[2.**-16]])]*2)
                       for key in labels}
            counts={key:[u+1 for u in range(257)] for key in labels}
            shells={key:[Q(1)]*257 for key in labels}
            args=SimpleNamespace(groups=1,target_bits=40,max_splits=0,retain_parents=True,
                                 joint_witness=False,joint_top=2,probe_supports=[],probe_vector=[],
                                 screen_only=False,precision=192)
            with redirect_stdout(StringIO()):
                result=controller.cover(args,operators,counts,[1],shells,prefix_rank=True)
            self.assertTrue(0<result<arb(2)**-40)
        finally:
            ctx.prec=old_precision


if __name__=='__main__':
    unittest.main()
