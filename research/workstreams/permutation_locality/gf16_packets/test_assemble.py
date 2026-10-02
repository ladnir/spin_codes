import unittest
from copy import deepcopy
from types import SimpleNamespace
from unittest.mock import patch
from flint import arb
import assemble


class AssemblyScopeTests(unittest.TestCase):
    def test_metadata_rejects_wrong_ensemble_or_incomplete_scope(self):
        valid=dict(schema=assemble.dense.SCHEMA,updates=2,minimum_groups=65,threshold=62914,
                   unresolved={},leaves={'':{}},kernel='refresh')
        assemble.validate(valid)
        assemble.validate(dict(valid,updates=3))
        assemble.validate(dict(valid,updates=4))
        for key,value in [('schema','old-lane-ensemble'),('updates',1),('updates',5),('updates','3'),('minimum_groups',0),
                          ('minimum_groups',2049),('threshold',-1),('threshold',1<<21),
                          ('unresolved',{'0':{}}),('unresolved',None),('leaves',{}),
                          ('screen_only',True),('kernel','unknown'),('base_tilt','0'),
                          ('central_bits',0),('central_bits',257),('central_scale','0'),
                          ('central_scale','-1/2'),('row_bias','0'),('row_bias','1/2'),('row_parity','true'),('row_parity',1),
                          ('variance_shuffle','true'),('variance_shuffle',1),
                          ('variance_bins',-1),('variance_bins',65),('variance_bins',True),('variance_bins',1),
                          ('regional_count',1),('regional_count','true'),('regional_count',True)]:
            record=deepcopy(valid);record[key]=value
            with self.assertRaises(ValueError):assemble.validate(record)

    def test_exact_single_group_and_remaining_ranges_are_disjoint(self):
        args=SimpleNamespace(single_group_exact=True,threshold=209715,sparse_tilts=['.00032'],
            precision=256,max_splits=100,sparse_inner='gf-rank',shape_penalties=['1'],
            shape_weight_tilt='1',exact_feedback=True)
        with patch('single_group.run',return_value=arb(2)**-42) as one, \
                patch('assemble.sparse_cover.run',return_value=arb(2)**-52) as rest:
            self.assertEqual(assemble.regenerate_sparse(1,args),0)
            one.assert_not_called();rest.assert_not_called()
            self.assertEqual(assemble.regenerate_sparse(2,args),arb(2)**-42)
            rest.assert_not_called()
            self.assertEqual(assemble.regenerate_sparse(5,args),arb(2)**-42+arb(2)**-52)
            self.assertEqual(rest.call_args.args[0],[2,3,4])
            rest.return_value=None
            self.assertIsNone(assemble.regenerate_sparse(5,args))
        args.single_group_exact=False
        with patch('single_group.run') as one, patch('assemble.sparse_cover.run',return_value=arb(2)**-52) as rest:
            self.assertEqual(assemble.regenerate_sparse(5,args),arb(2)**-52)
            self.assertEqual(rest.call_args.args[0],[1,2,3,4])
            one.assert_not_called()

    def test_separate_small_grid_still_covers_every_sparse_occupancy(self):
        args=SimpleNamespace(single_group_exact=True,threshold=188743,sparse_tilts=['.032'],
            sparse_small_tilts=['.001','.0024'],sparse_small_through=4,sparse_target_bits=48,
            precision=256,max_splits=150,sparse_inner='gf-birth-classes',shape_penalties=['1'],updates=3,
            shape_weight_tilt='1',exact_feedback=True)
        args.sparse_joint_return_through=3;args.sparse_lazy_density_through=6
        args.sparse_analytic_gradient=True
        with patch('single_group.run',return_value=arb(2)**-60) as one, \
                patch('assemble.sparse_cover.run',return_value=arb(2)**-50) as rest:
            upper=assemble.regenerate_sparse(8,args)
            self.assertEqual(upper,arb(2)**-60+arb(2)**-49)
            self.assertEqual(one.call_args.args[1],['.001','.0024','.032'])
            self.assertEqual(one.call_args.kwargs['updates'],3)
            calls=rest.call_args_list
            self.assertEqual([call.args[0] for call in calls],[[2,3,4],[5,6,7]])
            self.assertEqual([call.args[2] for call in calls],[['.001','.0024'],['.032']])
            self.assertEqual([call.args[5] for call in calls],[48,48])
            self.assertEqual([call.kwargs['updates'] for call in calls],[3,3])
            self.assertEqual([call.kwargs['joint_return_through'] for call in calls],[3,3])
            self.assertEqual([call.kwargs['lazy_density_through'] for call in calls],[6,6])
            self.assertEqual([call.kwargs['analytic_gradient'] for call in calls],[True,True])
            rest.reset_mock();args.updates=4
            self.assertEqual(assemble.regenerate_sparse(8,args),upper)
            self.assertEqual(one.call_args.kwargs['updates'],4)
            self.assertEqual([call.kwargs['updates'] for call in rest.call_args_list],[4,4])
            rest.reset_mock();rest.return_value=None
            self.assertIsNone(assemble.regenerate_sparse(8,args))
            self.assertEqual(rest.call_count,1)

    def test_small_grid_requires_a_range_and_positive_tilts(self):
        for through,tilts,bits in ((4,None,48),(0,['.001'],48),(4,['0'],48),(4,['.001'],39)):
            args=SimpleNamespace(sparse_small_through=through,sparse_small_tilts=tilts,sparse_target_bits=bits)
            with self.assertRaises(ValueError):assemble.regenerate_sparse(8,args)

    def test_parallel_chunks_partition_without_gaps_or_duplicates(self):
        for lo,hi in ((1,2),(2,4),(2,9),(9,97),(49,97),(2000,2049)):
            for workers in range(1,9):
                chunks=assemble.sparse_chunks(lo,hi,workers)
                self.assertEqual(len(chunks),min(workers,hi-lo))
                self.assertTrue(all(a<b for a,b in chunks))
                self.assertEqual([q for a,b in chunks for q in range(a,b)],list(range(lo,hi)))
        for lo,hi,workers in ((0,2,2),(2,2,2),(1,2050,2),(1,3,0),(1,3,9),(1,3,True)):
            with self.assertRaises(ValueError):assemble.sparse_chunks(lo,hi,workers)

    def test_worker_returns_only_fresh_outward_result(self):
        args=SimpleNamespace()
        with patch('assemble.sparse_range',return_value=arb(2)**-50) as fresh:
            task=(2,5,['.001'],args,None)
            result=assemble.sparse_worker(task)
            fresh.assert_called_once_with(2,5,['.001'],args)
            self.assertEqual(result[:2],(2,5))
            self.assertEqual(arb(result[2][0])*arb(2)**result[2][1],arb(2)**-50)
            fresh.return_value=None
            self.assertEqual(assemble.sparse_worker(task),(2,5,None))
            fresh.return_value=arb(0)
            with self.assertRaises(ArithmeticError):assemble.sparse_worker(task)

    def test_parallel_sum_requires_every_assigned_range(self):
        args=SimpleNamespace(sparse_log=None)
        ranges=[(2,4,['.001']),(4,6,['.032'])]
        def fresh_results(function,tasks):
            self.assertIs(function,assemble.sparse_worker)
            return [(lo,hi,(1,-50)) for lo,hi,*_ in tasks]
        with patch('concurrent.futures.ProcessPoolExecutor') as executor:
            pool=executor.return_value.__enter__.return_value
            pool.map.side_effect=fresh_results
            self.assertEqual(assemble.parallel_sparse(ranges,args,2),arb(2)**-48)
            tasks=pool.map.call_args.args[1]
            self.assertEqual([(t[0],t[1]) for t in tasks],[(2,3),(3,4),(4,5),(5,6)])
            self.assertEqual([t[2] for t in tasks],[['.001'],['.001'],['.032'],['.032']])
            for bad in ([(2,3,(1,-50))],[(2,4,(1,-50))],
                        [(2,3,(0,-50))],[(2,3,(1.0,-50))]):
                pool.map.side_effect=None;pool.map.return_value=bad
                with self.assertRaises(ArithmeticError):assemble.parallel_sparse(ranges,args,2)
            pool.map.return_value=[(2,3,(1,-50)),(3,4,None),(4,5,(1,-50)),(5,6,(1,-50))]
            self.assertIsNone(assemble.parallel_sparse(ranges,args,2))
        with patch('concurrent.futures.ProcessPoolExecutor') as executor:
            self.assertEqual(assemble.parallel_sparse([(2,2,[])],args,2),0)
            executor.assert_not_called()

    def test_parallel_assembly_preserves_exact_single_group_and_scope(self):
        args=SimpleNamespace(single_group_exact=True,threshold=193986,sparse_tilts=['.032'],
            sparse_small_tilts=['.001'],sparse_small_through=4,sparse_target_bits=48,
            precision=256,sparse_inner='gf-birth-classes',sparse_workers=4)
        with patch('single_group.run',return_value=arb(2)**-45), \
                patch('assemble.parallel_sparse',return_value=arb(2)**-46) as rest:
            self.assertEqual(assemble.regenerate_sparse(8,args),arb(3)*arb(2)**-46)
            self.assertEqual(rest.call_args.args,([(2,5,['.001']),(5,8,['.032'])],args,4))
            rest.return_value=None
            self.assertIsNone(assemble.regenerate_sparse(8,args))
        for workers in (0,9,True,1.0):
            args.sparse_workers=workers
            with self.assertRaises(ValueError):assemble.regenerate_sparse(8,args)

    def test_spawn_path_selects_this_assembler_and_restores_on_exception(self):
        original=list(assemble.sys.path)
        with self.assertRaisesRegex(RuntimeError,'test exception'):
            with assemble.sparse_spawn_path():
                self.assertEqual(assemble.sys.path[0],str(assemble.Path(assemble.__file__).resolve().parent))
                raise RuntimeError('test exception')
        self.assertEqual(assemble.sys.path,original)


if __name__=='__main__':unittest.main()
