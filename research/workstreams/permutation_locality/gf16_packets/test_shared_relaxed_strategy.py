import copy
import hashlib
import json
import unittest
from fractions import Fraction as Q
from types import SimpleNamespace
from unittest.mock import patch

from flint import arb, ctx
import shared_relaxed_strategy as s


class SharedAssemblyTests(unittest.TestCase):
    def record(self):
        return dict(schema=s.DENSE_SCHEMA, updates=2, K=1 << 20, N=1 << 21,
            minimum_groups=33, zero_bits=64, variance_bins=16,
            base_tilt='3/16', cost_tilt='1/4', results=[dict(distance='1/20',
                threshold=104857, root=['0', '1'], cover=dict(
                    leaves={'0':dict(witness={'fresh':0,'parameters':['1/100','0','0']}),
                            '1':dict(witness={'fresh':1,'parameters':['1/200','0','0']})},
                    unresolved={}))])

    def rows(self):
        return [dict(occupancy=q, updates=2, upper=[1,-40]) for q in range(1,33)]

    def test_schema_update_count_and_geometry(self):
        s.validate_record(self.record(), complete=True)
        for field,value in (('schema','pairwise-gf16-cover-1'),('updates',4),('updates',True),
                ('K',1 << 18),('N',1 << 20),('minimum_groups',True),('minimum_groups',0),
                ('zero_bits',-1),('variance_bins',0),('cost_tilt','0'),('base_tilt','0'),
                ('pruned',1),('pruned','true'),('pruned',None)):
            record=self.record();record[field]=value
            with self.assertRaises(ValueError):s.validate_record(record)

    def test_cutoff_matches_declared_strict_distance(self):
        for field,value in (('threshold',104856),('threshold',True),('distance','0'),
                ('distance','1/2'),('root',['0']),('root',['1','0'])):
            record=self.record();record['results'][0][field]=value
            with self.assertRaises(ValueError):s.validate_record(record)

    def test_dense_partition_holes_overlaps_and_empty_witness(self):
        for cover in (dict(leaves={'0':dict(witness={})},unresolved={}),
                dict(leaves={'':dict(witness={}), '0':dict(witness={})},unresolved={}),
                dict(leaves={'':dict(witness={})},unresolved={'1':{}}),
                dict(leaves={},unresolved={}),dict(leaves={'':{}},unresolved={})):
            record=self.record();record['results'][0]['cover']=cover
            with self.assertRaises(ValueError):s.validate_record(record,complete=True)

    def test_validate_before_expensive_model_build(self):
        with patch.object(s.shared_mixture,'actual_components') as rebuild:
            record=self.record();record['updates']=4
            with self.assertRaises(ValueError):s.build_model(record)
            for precision in (True,128.0,127):
                with self.assertRaises(ValueError):s.build_model(self.record(),precision)
            rebuild.assert_not_called()

    def test_pruned_exact_candidate_can_reduce_and_remove_components(self):
        caps=[0,2,1];original=[(Q(8),Q(1,2)),(Q(8),Q(3,4))]
        candidate=[dict(mass='4',activity='1/2')]
        record=dict(pruned=True,mixture=candidate,pruning=dict(success=False))
        self.assertEqual(s.checked_mixture(record,caps,original),[(Q(4),Q(1,2))])
        unchanged=[dict(mass=str(c),activity=str(p)) for c,p in original]
        self.assertEqual(s.checked_mixture(dict(mixture=unchanged),caps,original),original)
        self.assertEqual(s.checked_mixture(dict(pruned=False,mixture=unchanged),caps,original),original)
        with self.assertRaises(ValueError):s.checked_mixture(dict(mixture=candidate),caps,original)

    def test_pruned_invalid_or_enlarged_candidates_rejected(self):
        caps=[0,2,1];original=[(Q(8),Q(1,2)),(Q(8),Q(3,4))]
        for candidate in (None,[],[dict(mass='4')],[dict(mass='4',activity='1/2',extra=1)],
                [dict(mass=True,activity='1/2')],[dict(mass=4.0,activity='1/2')],
                [dict(mass='4',activity=.5)],[dict(mass='0',activity='1/2')],
                [dict(mass='-1',activity='1/2')],[dict(mass='9',activity='1/2')],
                [dict(mass='4',activity='1/3')],[dict(mass='1/0',activity='1/2')],
                [dict(mass='nan',activity='1/2')],
                [dict(mass='4',activity='1/2'),dict(mass='1',activity='2/4')]):
            with self.assertRaises(ValueError):
                s.checked_mixture(dict(pruned=True,mixture=candidate),caps,original)
        # Allowed component geometry does not suffice: every shell is checked.
        with self.assertRaises(ArithmeticError):
            s.checked_mixture(dict(pruned=True,mixture=[dict(mass='3',activity='1/2')]),caps,original)

    def test_pruned_model_rebuilds_caps_and_checks_saved_witness(self):
        caps=[0,2,1];original=[(Q(8),Q(1,2)),(Q(8),Q(3,4))]
        record=self.record();record.update(pruned=True,mixture=[dict(mass='4',activity='1/2')],
            cap_sha256=hashlib.sha256(json.dumps(list(map(str,caps))).encode()).hexdigest())
        data=dict(bits=19,windows=32,updates=2)
        model=SimpleNamespace(root=(Q(0),Q(1)))
        prior=ctx.prec
        try:
            with patch.object(s.shared_mixture,'actual_components',return_value=([],caps,original)) as fresh,\
                    patch.object(s.birth_classes,'actual',return_value=data),\
                    patch.object(s.sc,'Model',return_value=model) as construct:
                self.assertIs(s.build_model(record,384),model)
                fresh.assert_called_once_with(coupled=True,zero_bits=64,cost_tilt=Q(1,4))
                self.assertEqual(construct.call_args.args[0],s.shared_mixture.as_components([(Q(4),Q(1,2))]))
                construct.reset_mock()
                record['cap_sha256']='wrong'
                with self.assertRaises(ValueError):s.build_model(record,384)
                construct.assert_not_called()
        finally:ctx.prec=prior

    def test_sparse_coverage_and_routing_scope(self):
        self.assertEqual(s.sum_sparse(self.rows(),32),Q(32,1 << 40))
        self.assertEqual(s.sum_sparse([],0),0)
        for rows in (self.rows()[:-1],self.rows()[:-1]+[self.rows()[0]],
                self.rows()[:-1]+[dict(occupancy=32,updates=4,upper=[1,-40])],
                self.rows()[:-1]+[dict(occupancy=32,updates=2,upper=None)],
                self.rows()[:-1]+[dict(occupancy=32,updates=2,upper=[0,-40])]):
            with self.assertRaises(ValueError):s.sum_sparse(rows,32)

    def test_exact_sum_and_twenty_bit_strict_margin(self):
        dense=dict(complete=True,upper=[1,-40])
        result=s.combine(dense,self.rows(),32,20)
        self.assertEqual(s.dyadic(result['total_upper']),Q(33,1 << 40))
        equal=dict(dense,upper=s.endpoint(Q(2)**-20-Q(32,1 << 40)))
        with self.assertRaises(ValueError):s.combine(equal,self.rows(),32,20)
        for bits in (0,True,20.0):
            with self.assertRaises(ValueError):s.combine(dense,self.rows(),32,bits)
        with self.assertRaises(ValueError):s.endpoint(Q(1,3))

    def test_compact_aggregation_encloses_widely_separated_exponents(self):
        exact=Q(2)**-40+Q(2)**-100000
        compact=s.compact_endpoint(exact,128)
        self.assertGreaterEqual(s.dyadic(compact),exact)
        self.assertLessEqual(compact[0].bit_length(),128)
        self.assertLess(len(json.dumps(compact)),64)
        self.assertEqual(s.dyadic(s.compact_endpoint(Q(2)**-100000,128)),Q(2)**-100000)
        with self.assertRaises(ValueError):s.compact_endpoint(exact,True)
        with self.assertRaises(ValueError):s.compact_endpoint(exact,31)
        dense=dict(complete=True,upper=[1,-40])
        rows=[dict(occupancy=1,updates=2,upper=[1,-100000])]
        result=s.combine(dense,rows,1,20)
        self.assertGreaterEqual(s.dyadic(result['total_upper']),exact)
        self.assertGreaterEqual(s.dyadic(result['total_upper']),
            s.dyadic(result['dense_upper'])+s.dyadic(result['sparse_upper']))
        self.assertLess(len(json.dumps(result)),600)

    def test_saved_endpoints_and_acceptance_labels_not_trusted(self):
        class Model:
            root=(Q(0),Q(1))
            def outward(self,cell,witness):
                self.calls.append((cell,witness,ctx.prec))
                return arb(2)**-40
        model=Model();model.calls=[]
        record=self.record()
        for leaf in record['results'][0]['cover']['leaves'].values():
            leaf.update(upper=[1,-999999],accepted=True)
        prior=ctx.prec
        try:
            with patch.object(s,'build_model',return_value=model):
                result=s.replay_dense(record,384)
            self.assertEqual(s.dyadic(result['upper']),Q(2,1 << 40))
            self.assertEqual(model.calls,[((Q(0),Q(1,2)),{'fresh':0,'parameters':['1/100','0','0']},384),
                                          ((Q(1,2),Q(1)),{'fresh':1,'parameters':['1/200','0','0']},384)])
        finally:ctx.prec=prior

    def test_precision_mutation_during_outward_is_rejected(self):
        class Model:
            root=(Q(0),Q(1))
            def outward(self,cell,witness):
                ctx.prec=128
                return arb(2)**-40
        prior=ctx.prec
        try:
            with patch.object(s,'build_model',return_value=Model()):
                with self.assertRaises(ArithmeticError):s.replay_dense(self.record(),384)
        finally:ctx.prec=prior

    def test_verifier_freshly_regenerates_sparse_and_checks_sum(self):
        dense=dict(complete=True,upper=[1,-40],leaves=2,checked=[])
        prior=ctx.prec;ctx.prec=384
        try:
            with patch.object(s,'replay_dense',return_value=dense),patch.object(s.shared_sparse,'run',return_value=self.rows()) as run:
                result=s.verify(self.record(),precision=384,bits=20)
                self.assertEqual(run.call_args.args[0],[2])
                self.assertEqual(run.call_args.args[1],list(range(1,33)))
                self.assertEqual(run.call_args.args[2],104857)
                self.assertTrue(run.call_args.kwargs['coupled_counts'])
                self.assertTrue(result['complete'])
                self.assertEqual(result['ensemble'],s.ENSEMBLE)
                self.assertEqual(result['minimum_distance'],104858)
                self.assertEqual(result['dense_occupancies'],[33,2048])
                self.assertEqual(s.dyadic(result['total_upper']),Q(33,1 << 40))
        finally:ctx.prec=prior

    def test_failed_dense_budget_skips_sparse_work(self):
        prior=ctx.prec;ctx.prec=384
        try:
            with patch.object(s,'replay_dense',return_value=dict(complete=True,upper=[1,-20])),patch.object(s.shared_sparse,'run') as run:
                with self.assertRaises(ValueError):s.verify(self.record())
                run.assert_not_called()
        finally:ctx.prec=prior

    def test_full_verification_rejects_active_search_cache(self):
        import regional_count
        import shared_relaxed_replay
        def cached(*args):raise AssertionError('cache must not be used')
        cached._shared_search_region_cache=True
        with patch.object(regional_count,'outward',cached), patch.object(s,'build_model') as build:
            for action in (lambda:s.verify(self.record()), lambda:s.replay_dense(self.record()),
                           lambda:shared_relaxed_replay.replay_dense(self.record())):
                with self.assertRaises(ArithmeticError):action()
            build.assert_not_called()

    def test_optional_fresh_dense_replayer_preserves_assembly_contract(self):
        dense=dict(complete=True,upper=[1,-40],leaves=2,checked=[])
        calls=[]
        def fresh(record,precision,index):
            calls.append((record,precision,index));ctx.prec=precision
            return dense
        prior=ctx.prec
        try:
            record=self.record()
            with patch.object(s,'replay_dense') as serial,\
                    patch.object(s.shared_sparse,'run',return_value=self.rows()):
                result=s.verify(record,precision=384,dense_replayer=fresh)
                serial.assert_not_called()
                self.assertEqual(calls,[(record,384,0)])
                self.assertTrue(result['complete'])
                self.assertEqual(s.dyadic(result['total_upper']),Q(33,1 << 40))
            with self.assertRaises(ValueError):s.verify(record,dense_replayer=True)
        finally:ctx.prec=prior

    def test_fixed_witness_cutoff_scaling_and_frontier(self):
        dense=dict(complete=True, ensemble=s.ENSEMBLE, threshold=100, leaves=2,
            upper=[2,-40], checked=[dict(path='0',upper=[1,-40],output_tilt='1/100'),
                                   dict(path='1',upper=[1,-40],output_tilt='1/200')])
        prior=ctx.prec
        try:
            self.assertEqual(s.dyadic(s.fixed_witness_upper(dense,100)),Q(2,1 << 40))
            lower=s.dyadic(s.fixed_witness_upper(dense,50))
            original=s.dyadic(s.fixed_witness_upper(dense,100))
            higher=s.dyadic(s.fixed_witness_upper(dense,150))
            self.assertLess(lower,original);self.assertLess(original,higher)
            # Independent direct expression, using an interval rather than a float.
            expected=arb(2)**-40*((arb(1)/2).exp()+(arb(1)/4).exp())
            got=arb(higher.numerator)/arb(higher.denominator)
            self.assertTrue(got >= expected)
            frontier=s.fixed_witness_frontier(dense,20)
            self.assertLess(s.dyadic(frontier['upper']),Q(2)**-20)
            self.assertGreaterEqual(s.dyadic(s.fixed_witness_upper(dense,frontier['threshold']+1)),Q(2)**-20)
        finally:ctx.prec=prior

    def test_fixed_witness_scaling_rejects_scope_and_sum_mismatches(self):
        dense=dict(complete=True,ensemble=s.ENSEMBLE,threshold=100,leaves=1,
            upper=[1,-40],checked=[dict(path='',upper=[1,-40],output_tilt='1/100')])
        for changes in (dict(ensemble='independent4-gf16-r2'),dict(upper=[1,-41]),
                dict(leaves=2),dict(complete=False),
                dict(checked=[dict(path='0',upper=[1,-40],output_tilt='1/100')]),
                dict(checked=[dict(path='',upper=[1,-40],output_tilt='0')])):
            with self.assertRaises(ValueError):s.fixed_witness_upper(dict(dense,**changes),200)


if __name__ == '__main__':unittest.main()
