import unittest
from unittest.mock import patch
from contextlib import ExitStack
import numpy as np
from flint import arb,arb_mat
import sparse_cover


class SparseTransferScopeTests(unittest.TestCase):
    def test_bridge_cannot_import_weighted_lane_bounds(self):
        args=sparse_cover.build_args(64,['.0032'],256,100,52,None)
        self.assertEqual(args.penalties,['1'])
        self.assertEqual(args.weight_tilt,'1')
        self.assertIsNone(args.operator_cache)
        self.assertFalse(args.screen_only)
        self.assertTrue(args.retain_parents)
        self.assertEqual(args.check_interval,1)
        self.assertEqual(args.updates,2)
        self.assertEqual(sparse_cover.build_args(1,['.001'],192,1,48,None,updates=3).updates,3)
        self.assertEqual(sparse_cover.build_args(1,['.001'],192,1,48,None,updates=4).updates,4)

    def test_reject_invalid_scope_before_computing(self):
        for occupancies,threshold,tilts,precision in [
                ([],0,['.0032'],256),([1,1],0,['.0032'],256),
                ([0],0,['.0032'],256),([2049],0,['.0032'],256),
                ([1],1<<21,['.0032'],256),([1],0,['0'],256),
                ([1],0,['.0032'],64)]:
            with self.assertRaises(ValueError):
                sparse_cover.run(occupancies,threshold,tilts,precision,100,52)
        for degree in (0,15,2049,16.5):
            with self.assertRaises(ValueError):
                sparse_cover.run([16],0,['.0032'],256,100,52,operator_through=degree)
        for inner in ('shape-uniform','gf-shape-tilt'):
            with self.assertRaises(ValueError):
                sparse_cover.run([16],0,['.0032'],256,100,52,inner=inner,exact_feedback=True)
        for inner in ('gf-rank','gf-birth-classes'):
            with self.assertRaises(ValueError):
                sparse_cover.run([16],0,['.0032'],256,100,52,inner=inner)
        for inner in ('gf-occupancy','gf-density','gf-density-refined','gf-classes','gf-rank','gf-birth-classes'):
            with self.assertRaises(ValueError):
                sparse_cover.run([16],0,['.0032'],256,100,52,directory='old-cache',inner=inner,exact_feedback=True)
        for updates,inner in ((1,'gf-rank'),(5,'gf-birth-classes'),(True,'gf-rank'),
                              (3,'gf-occupancy'),(3,'gf-density'),(3,'shape-uniform'),
                              (4,'gf-occupancy'),(4,'gf-density'),(4,'shape-uniform')):
            with self.assertRaises(ValueError):
                sparse_cover.run([1],0,['.0032'],192,10,48,inner=inner,updates=updates,exact_feedback=True)
        for name,limit in (('joint_return_through',4),('lazy_density_through',32)):
            for value in (-1,limit+1,True,'1'):
                with self.assertRaises(ValueError):
                    sparse_cover.run([1],0,['.0032'],192,10,48,inner='gf-birth-classes',
                                     exact_feedback=True,**{name:value})
            with self.assertRaises(ValueError):
                sparse_cover.run([1],0,['.0032'],192,10,48,inner='gf-rank',exact_feedback=True,**{name:1})

    def test_birth_classes_use_all_mass_coordinates(self):
        import occupancy_birth_classes
        n=8
        matrix=arb_mat([[int(i==j) for j in range(n)] for i in range(n)])
        operators={('.0032','1'):([matrix],[np.eye(n)])}
        with ExitStack() as stack:
            for name in ('mixing_test','folding_test','geometry_test','retained_test','placement_prefix_test'):
                stack.enter_context(patch.object(sparse_cover.sparse,name))
            stack.enter_context(patch.object(sparse_cover.sparse,'authenticated_caps',return_value={}))
            stack.enter_context(patch.object(sparse_cover.sparse,'weighted_cdf_upper',return_value=[1]))
            stack.enter_context(patch.object(sparse_cover.sparse,'integer_cdf',return_value=[1]))
            build=stack.enter_context(patch.object(occupancy_birth_classes,'build_operators',return_value=operators))
            cover=stack.enter_context(patch.object(sparse_cover.sparse,'cover',return_value=arb(2)**-60))
            result=sparse_cover.run([1],12,['.0032'],192,10,52,inner='gf-birth-classes',exact_feedback=True)
            self.assertEqual(result,arb(2)**-60)
            self.assertTrue(build.call_args.args[0].exact_feedback)
            np.testing.assert_array_equal(cover.call_args.args[3],np.ones(n))
            self.assertEqual(cover.call_args.kwargs['cutoff'],12)
            sparse_cover.run([1],12,['.0032'],192,10,52,inner='gf-birth-classes',exact_feedback=True,updates=3)
            self.assertEqual(build.call_args.args[0].updates,3)
            sparse_cover.run([1],12,['.0032'],192,10,52,inner='gf-birth-classes',exact_feedback=True,updates=4)
            self.assertEqual(build.call_args.args[0].updates,4)
            sparse_cover.run([1],12,['.0032'],192,10,52,inner='gf-birth-classes',exact_feedback=True,
                             updates=3,joint_return_through=3,lazy_density_through=6)
            self.assertEqual(build.call_args.args[0].joint_return_through,3)
            self.assertEqual(build.call_args.args[0].lazy_density_through,6)
            sparse_cover.run([1],12,['.0032'],192,10,52,inner='gf-birth-classes',exact_feedback=True,
                             analytic_gradient=True)
            self.assertIs(cover.call_args.args[0].analytic_gradient,True)
            with self.assertRaises(ValueError):
                sparse_cover.run([1],12,['.0032'],192,10,52,analytic_gradient='yes')


if __name__=='__main__':unittest.main()
