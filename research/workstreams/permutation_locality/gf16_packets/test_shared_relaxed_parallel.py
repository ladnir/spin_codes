from fractions import Fraction as Q
import copy
from types import SimpleNamespace
from unittest.mock import patch
import unittest

from flint import arb,ctx
import shared_relaxed_parallel as p


class SharedParallelTests(unittest.TestCase):
    def model(self):
        return SimpleNamespace(root=(Q(0),Q(1)),components=[(Q(1),Q(0),0),(Q(100),Q(1,2),1)])

    def record(self):
        return dict(schema=p.proof.DENSE_SCHEMA,updates=2,K=1 << 20,N=1 << 21,
            minimum_groups=33,zero_bits=64,variance_bins=16,base_tilt='3/16',cost_tilt='1/4',
            pruned=True,mixture=[dict(mass='100',activity='1/2')],cap_sha256='f'*64,
            results=[dict(distance='3/50',threshold=125829,root=['0','1'],
                probes=[dict(mean='1/8',upper=[1,-40])],rechecked_precision=256,
                cover=dict(leaves={'0':dict(witness={'old':1})},unresolved={'1':{}},visited=40)),
                dict(distance='1/20',threshold=104857,root=['0','1'],probes=[])])

    def test_retarget_lower_and_higher_preserves_source_and_comparison(self):
        original=self.record();saved=copy.deepcopy(original)
        for distance,cutoff in (('1/20',104857),('7/100',146800)):
            result=p.retarget_record(original,distance,0,'a'*64)
            self.assertEqual(original,saved)
            self.assertEqual(result['results'][0]['distance'],distance)
            self.assertEqual(result['results'][0]['threshold'],cutoff)
            self.assertEqual(result['results'][0]['cover'],original['results'][0]['cover'])
            self.assertEqual(result['results'][0]['root'],original['results'][0]['root'])
            self.assertEqual(result['results'][1],original['results'][1])
            for key in ('updates','K','N','minimum_groups','zero_bits','variance_bins',
                        'base_tilt','cost_tilt','pruned','mixture','cap_sha256'):
                self.assertEqual(result[key],original[key])
            self.assertEqual(result['results'][0]['probes'],[])
            self.assertNotIn('rechecked_precision',result['results'][0])
            source=result['results'][0]['retargeted_from']
            self.assertEqual(source,dict(source_sha256='a'*64,result_index=0,distance='3/50',
                threshold=125829,prior_probes=original['results'][0]['probes']))
            self.assertIn('Retargeted claim only',result['proof_status'])
            # Witnesses remain proposals; no prior acceptance is counted.
            model=self.model()
            result=p.search(model,result['results'][0]['cover'],self.mapping(model,False),
                            2,2,8,40,lambda _:None)
            self.assertEqual(result['leaves'],{})
            self.assertTrue(result['unresolved'])

    def test_retarget_rejects_invalid_distance_and_provenance(self):
        for distance in ('0','1/2','-1/20','1','1/0','nan',True,.07):
            with self.assertRaises(ValueError):p.retarget_record(self.record(),distance,0,'a'*64)
        for digest in ('','x'*64,'A'*64,None):
            with self.assertRaises(ValueError):p.retarget_record(self.record(),'7/100',0,digest)

    def test_retained_witness_is_rechecked_not_imported(self):
        class Model:
            def outward(self,cell,witness):
                self.seen=(cell,witness,ctx.prec)
                return arb(2)**-50
            def proposal(self,cell):raise AssertionError('unneeded proposal')
        m=Model();prior=ctx.prec
        try:
            row=p.evaluate(m,('0',(Q(0),Q(1,2)),40,{'keep':1}),256,'digest')
            self.assertEqual(m.seen,((Q(0),Q(1,2)),{'keep':1},256))
            self.assertEqual(row['upper'],[1,-50])
            self.assertEqual(row['component_sha256'],'digest')
        finally:ctx.prec=prior

    def test_proposal_precision_reset_and_outward_target(self):
        class Model:
            def proposal(self,cell):ctx.prec=192;return -50.,{'new':1}
            def outward(self,cell,witness):
                self.seen=ctx.prec
                return arb(2)**-45
        m=Model();prior=ctx.prec
        try:
            row=p.evaluate(m,('',(Q(0),Q(1)),40,None),256,'digest')
            self.assertEqual(m.seen,256)
            self.assertEqual(row['upper'],[1,-45])
            with patch.object(m,'outward',return_value=arb(2)**-39):
                self.assertIsNone(p.evaluate(m,('',(Q(0),Q(1)),40,None),256,'digest')['upper'])
        finally:ctx.prec=prior

    def test_mutation_during_fresh_bound_is_rejected(self):
        class Model:
            def outward(self,*args):ctx.prec=128;return arb(2)**-50
        prior=ctx.prec
        try:
            with self.assertRaises(ArithmeticError):
                p.evaluate(Model(),('',(Q(0),Q(1)),40,{}),256,'digest')
        finally:ctx.prec=prior

    def test_midpoint_success_only_requests_bisection(self):
        class Model:
            split_on_midpoint=True
            def __init__(self):self.calls=[]
            def outward(self,cell,witness):
                self.calls.append(cell)
                return arb(2)**(-80 if cell[0]==cell[1] else 10)
            def proposal(self,cell):raise AssertionError('midpoint-guided bisection needs no new search')
        model=Model();old=dict(parameters=['1/100','0','0'])
        prior=ctx.prec
        try:
            row=p.evaluate(model,('0',(Q(1,4),Q(1,2)),32,old),256,'digest')
            self.assertIsNone(row['upper'])
            self.assertEqual(row['proposal'],10.)
            self.assertEqual(row['witness'],old)
            self.assertEqual(model.calls,[(Q(1,4),Q(1,2)),(Q(3,8),Q(3,8))])
        finally:ctx.prec=prior

    def test_point_reuse_requires_fresh_finite_cell_success(self):
        class Model:
            point_probes=[dict(mean='1/3',witness=dict(parameters=['1/100','0','0'],point=True),
                               upper=[1,-99999])]
            def outward(self,cell,witness):
                self.seen=cell
                return arb(2)**(-60 if witness.get('point') else 100)
            def proposal(self,cell):raise AssertionError('fresh finite cell already passed')
        model=Model();prior=ctx.prec
        try:
            row=p.evaluate(model,('0',(Q(1,4),Q(1,2)),32,{}),256,'digest')
            self.assertEqual(model.seen,(Q(1,4),Q(1,2)))
            self.assertEqual(row['upper'],[1,-60])
            self.assertTrue(row['witness']['point'])
            self.assertNotIn('upper',row['witness'])
        finally:ctx.prec=prior

    def test_failed_midpoint_still_runs_proposal(self):
        class Model:
            split_on_midpoint=True
            def outward(self,cell,witness):return arb(2)**(-70 if witness.get('new') else 5)
            def proposal(self,cell):self.proposed=True;return -80.,dict(new=True)
        model=Model();prior=ctx.prec
        try:
            row=p.evaluate(model,('0',(Q(1,4),Q(1,2)),32,{}),256,'digest')
            self.assertTrue(model.proposed)
            self.assertEqual(row['upper'],[1,-70])
        finally:ctx.prec=prior

    def mapping(self,model,accept=True,corruption=None):
        def mapping(fn,jobs):
            rows=[]
            for path,cell,bits,old in jobs:
                row=dict(path=path,cell=list(map(str,cell)),witness={'rechecked':old},
                    component_sha256=p.component_digest(model.components),proposal=-50. if accept else 2.,
                    upper=[1,-50] if accept else None)
                if corruption:row.update(corruption)
                rows.append(row)
            return rows
        return mapping

    def test_every_checkpoint_is_a_complete_partition(self):
        model=self.model();snapshots=[]
        cover=dict(leaves={},unresolved={'':{}},visited=5)
        state=p.search(model,cover,self.mapping(model,False),2,4,8,40,snapshots.append)
        self.assertEqual(state['visited'],9)
        self.assertEqual(state['leaves'],{})
        self.assertTrue(state['unresolved'])
        for snapshot in snapshots:
            cells=p.proof.sc.partition(model,snapshot['leaves'],snapshot['unresolved'])
            self.assertTrue(cells)

    def test_prior_leaves_rechecked_and_worker_result_bound(self):
        model=self.model();snapshots=[]
        cover=dict(leaves={'0':dict(witness={'prior':1})},unresolved={'1':{}},visited=2)
        state=p.search(model,cover,self.mapping(model),2,2,8,40,snapshots.append)
        self.assertEqual(state['unresolved'],{})
        self.assertEqual(state['leaves']['0']['witness'],{'rechecked':{'prior':1}})
        self.assertEqual(state['leaves']['1']['witness'],{'rechecked':None})
        self.assertEqual(state['visited'],4)
        for corruption in (dict(component_sha256='wrong'),dict(path='wrong'),dict(upper=[1,-40]),
                           dict(cell=['0','1']),dict(witness=None)):
            with self.assertRaises(ArithmeticError):
                p.search(model,cover,self.mapping(model,corruption=corruption),2,2,8,40,lambda _:None)

    def test_missing_result_and_invalid_budget_rejected(self):
        model=self.model();cover=dict(leaves={},unresolved={'':{}})
        with self.assertRaises(ArithmeticError):
            p.search(model,cover,lambda fn,jobs:[],2,2,8,40,lambda _:None)
        for workers,budget,depth,bits in ((True,2,8,40),(0,2,8,40),(2,0,8,40),(2,2,0,40),(2,2,8,19)):
            with self.assertRaises(ValueError):
                p.search(model,cover,self.mapping(model),workers,budget,depth,bits,lambda _:None)

    def partition_witness(self):
        parts=[dict(interval=[lo,hi],dual=['2','3',str(i)]) for i,(lo,hi) in
               enumerate((('0','1/16'),('1/16','1/8'),('1/8','1/4')))]
        regional=[dict(copy.deepcopy(part),mgf_witnesses=[{'tilt':'1/8','dual':['3','1','2']}])
                  for part in parts]
        return dict(tilt='3/16',parameters=['1/100','2','3'],variance_dual=['4','5'],
                    variance_partition=parts,regional_count_parts=regional)

    def test_restriction_clips_full_partition_and_keeps_matching_duals(self):
        original=self.partition_witness();saved=copy.deepcopy(original)
        child=(Q(1,16),Q(3,32))
        result=p.restrict_witness((Q(1,32),Q(3,4)),child,original)
        self.assertEqual(original,saved)
        self.assertEqual([v['interval'] for v in result['variance_partition']],
                         [['0','1/16'],['1/16','3/32']])
        self.assertEqual([v['interval'] for v in result['regional_count_parts']],
                         [['0','1/16'],['1/16','3/32']])
        self.assertEqual([v['dual'] for v in result['variance_partition']],
                         [v['dual'] for v in result['regional_count_parts']])
        self.assertEqual(result['regional_count_parts'][1]['mgf_witnesses'],
                         original['regional_count_parts'][1]['mgf_witnesses'])
        result['regional_count_parts'][0]['mgf_witnesses'][0]['tilt']='9'
        self.assertEqual(original,saved)

    def test_restriction_rejects_bad_geometry_and_mismatched_parts(self):
        parent=(Q(1,32),Q(3,4));child=(Q(1,16),Q(3,32))
        for bad_child in ((Q(0),Q(1,16)),(Q(3,32),Q(1,16)),(Q(1,2),Q(1))):
            with self.assertRaises(ValueError):p.restrict_witness(parent,bad_child,{})
        for key in ('variance_partition','regional_count_parts'):
            bad=self.partition_witness();bad[key][1]['interval']=['1/15','1/8']
            with self.assertRaises(ValueError):p.restrict_witness(parent,child,bad)
        bad=self.partition_witness();bad['regional_count_parts'][1]['dual'][0]='4'
        with self.assertRaises(ValueError):p.restrict_witness(parent,child,bad)
        bad=self.partition_witness();bad.pop('variance_partition')
        with self.assertRaises(ValueError):p.restrict_witness(parent,child,bad)
        bad=self.partition_witness();bad['regional_count_parts'][0]['mgf_witnesses']=[0]
        with self.assertRaises(ValueError):p.restrict_witness(parent,child,bad)

    def test_restriction_copies_generic_witness(self):
        original=dict(parameters=['1/100','2','3'],tilt='3/16',variance_dual=['4','5'])
        result=p.restrict_witness((0,1),(0,Q(1,2)),original)
        self.assertEqual(result,original)
        result['parameters'][0]='1/50'
        self.assertNotEqual(result,original)

    def test_parent_reuse_is_opt_in_and_children_are_freshly_evaluated(self):
        class Model:
            root=(Q(0),Q(1));components=[(Q(1),Q(0),0)]
            def __init__(self):self.calls=[]
            def proposal(self,cell):
                self.calls.append(('proposal',cell))
                return 1.,dict(parameters=['1/100','2','3'],tilt='3/16')
            def outward(self,cell,witness):
                self.calls.append(('outward',cell))
                return arb(2)**(-50 if cell[1]-cell[0]<=Q(1,2) else 1)
        for reuse in (False,True):
            model=Model();cover=dict(leaves={},unresolved={'':{}});snapshots=[]
            def mapping(fn,jobs):
                return [p.evaluate(model,job,256,p.component_digest(model.components)) for job in jobs]
            state=p.search(model,cover,mapping,2,3,8,40,snapshots.append,reuse)
            calls=[name for name,_ in model.calls]
            self.assertEqual(calls,['proposal','outward','outward'] if reuse else ['proposal']*3)
            self.assertEqual(len(state['leaves']),2 if reuse else 0)
            if reuse:
                self.assertTrue(all('witness' in v for v in snapshots[0]['unresolved'].values()))
        with self.assertRaises(ValueError):
            p.search(model,cover,mapping,2,3,8,40,lambda _:None,1)


if __name__=='__main__':unittest.main()
