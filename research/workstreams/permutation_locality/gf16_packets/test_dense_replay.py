from fractions import Fraction as Q
from types import SimpleNamespace
from unittest.mock import patch
import unittest
from flint import arb,ctx
import assemble


class DenseReplayTests(unittest.TestCase):
    def setUp(self):
        ctx.prec=192
        self.record=dict(schema=assemble.dense.SCHEMA,updates=4,threshold=207618,
            minimum_groups=49,kernel='birth-classes',unresolved={},
            leaves={'0':{'witness':{'fresh':0}},'1':{'witness':{'fresh':1}}})

    def test_sum_requires_exact_coverage_and_positive_dyadics(self):
        valid=[('0',(1,-50)),('1',(3,-52))]
        self.assertEqual(assemble.sum_dense_results(['0','1'],valid),arb(7)*arb(2)**-52)
        for results in (valid[:1],valid[::-1],valid+[valid[0]],
                [('0',(0,-50)),valid[1]],[('0',(1.,-50)),valid[1]],
                [('0',None),valid[1]],[('0',(True,-50)),valid[1]],
                [('0',(1,-50.)),valid[1]],[valid[0],valid[0]]):
            with self.assertRaises(ArithmeticError):assemble.sum_dense_results(['0','1'],results)
        with self.assertRaises(ArithmeticError):assemble.sum_dense_results(['0','0'],valid)

    def test_serial_and_parallel_replay_use_identical_result_order(self):
        class Model:
            root=(Q(0),Q(1))
            def outward(self,cell,witness):
                return arb(2)**-50 if witness['fresh']==0 else arb(3)*arb(2)**-52
        with patch('assemble.dense_model',return_value=Model()) as build:
            serial=assemble.replay_dense(self.record,192)
            self.assertEqual(build.call_args.kwargs,{'geometry_only':False})
            with patch('concurrent.futures.ProcessPoolExecutor') as executor:
                pool=executor.return_value.__enter__.return_value
                pool.map.return_value=[('0',(1,-50)),('1',(3,-52))]
                parallel=assemble.replay_dense(self.record,192,2)
                self.assertEqual(build.call_args.kwargs,{'geometry_only':True})
                self.assertEqual(pool.map.call_args.args,(assemble.dense_worker,['0','1']))
                self.assertIs(executor.call_args.kwargs['initializer'],assemble.initialize_dense_worker)
                self.assertEqual(executor.call_args.kwargs['initargs'],(self.record,192,None))
            self.assertEqual(serial,parallel)
        for workers in (0,5,True,1.0):
            with self.assertRaises(ValueError):assemble.replay_dense(self.record,192,workers)

    def test_worker_uses_fresh_witness_and_checks_precision(self):
        class Model:
            def outward(self,cell,witness):self.seen=(cell,witness);return arb(2)**-70
        model=Model()
        with patch.multiple(assemble,_dense_model=model,_dense_cells={'0':(Q(0),Q(1))},
                _dense_record=self.record,_dense_precision=192,_dense_log=None,create=True):
            self.assertEqual(assemble.dense_worker('0'),('0',(1,-70)))
            self.assertEqual(model.seen,((Q(0),Q(1)),{'fresh':0}))
            ctx.prec=128
            with self.assertRaises(ArithmeticError):assemble.dense_worker('0')
            ctx.prec=192
            def change_precision(cell,witness):ctx.prec=128;return arb(2)**-70
            with patch.object(model,'outward',side_effect=change_precision):
                with self.assertRaises(ArithmeticError):assemble.dense_worker('0')
        ctx.prec=192

    def test_initialization_builds_actual_maps_and_full_partition(self):
        model=SimpleNamespace(root=(Q(0),Q(1)))
        with patch('assemble.dense_model',return_value=model) as build:
            assemble.initialize_dense_worker(self.record,192,None)
            build.assert_called_once_with(self.record,192)
            self.assertEqual(assemble._dense_cells,{'0':(Q(0),Q(1,2)),'1':(Q(1,2),Q(1))})
            self.assertEqual(assemble._dense_precision,192)

    def test_model_rejects_invalid_claim_before_authentication(self):
        with patch('assemble.dense.actual_components') as authenticate:
            for precision in (127,True,192.0):
                with self.assertRaises(ValueError):assemble.dense_model(self.record,precision)
            for field,value in (('unresolved',{'0':{}}),('updates',5),('threshold',-1)):
                with self.assertRaises(ValueError):assemble.dense_model(dict(self.record,**{field:value}),192)
            authenticate.assert_not_called()


if __name__=='__main__':unittest.main()
