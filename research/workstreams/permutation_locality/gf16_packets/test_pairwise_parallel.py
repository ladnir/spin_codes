import copy
from fractions import Fraction as Q
import unittest
from unittest.mock import patch
from flint import arb,ctx
import pairwise_parallel as pp


class Model:
    root=(Q(0),Q(1));components=[(Q(1),Q(0),0),(Q(2),Q(1,2),1)]
    def proposal(self,cell):return (-100 if cell[1]-cell[0]<=Q(1,4) else 5),dict(fresh=True)
    def outward(self,cell,witness):return arb(2)**(-100 if witness.get('fresh') else -10)


class ParallelTests(unittest.TestCase):
    def test_cell_width_budget_sums_to_global_budget(self):
        for paths in ([''],['0','1'],['00','01','1'],['000','001','01','10','11']):
            self.assertEqual(sum((Q(2)**-pp.cell_target(p,64,44) for p in paths),Q(0)),Q(2)**-44)
        self.assertEqual(pp.cell_target('0101',64),64)
        for path,bits in (('2',44),('0',True),('0',39),('0','44')):
            with self.assertRaises(ValueError):pp.cell_target(path,64,bits)

    def test_adaptive_budget_checks_the_actual_cell_target(self):
        model=Model();record=self.record();record['unresolved']={'0':{},'1':{}}
        def mapping(unused,jobs):
            self.assertEqual([job[2] for job in jobs],[45,45])
            return [dict(path=path,cell=list(map(str,cell)),proposal=-47.,witness={},
                upper=[1,-46],component_sha256=pp.component_digest(model.components))
                for path,cell,bits,old in jobs]
        result=pp.search(model,record,mapping,2,2,4,64,lambda _:None,aggregate_bits=44)
        self.assertFalse(result['unresolved']);self.assertEqual(set(result['leaves']),{'0','1'})
        def wrong(unused,jobs):return [dict(r,upper=[1,-45]) for r in mapping(unused,jobs)]
        with self.assertRaises(ArithmeticError):
            pp.search(model,record,wrong,2,2,4,64,lambda _:None,aggregate_bits=44)

    def test_sufficient_proposal_still_requires_outward_check(self):
        model=Model()
        with patch.object(model,'proposal',return_value=(-67,{})),patch.object(model,'outward',return_value=arb(2)**-63) as outward:
            with self.assertRaises(ArithmeticError):
                pp.evaluate(model,('',model.root,64,None),ctx.prec,'digest')
            self.assertEqual(model.proposal_stop_bits,66)
            outward.assert_called_once()

    def record(self):
        return dict(schema=pp.cover.SCHEMA,ensemble=pp.cover.ENSEMBLE,
            parameters=copy.deepcopy(pp.cover.PARAMETERS),leaves={},unresolved={'':{}})

    def mapping(self,model):
        return lambda unused,jobs:[pp.evaluate(model,job,ctx.prec,pp.component_digest(model.components)) for job in jobs]

    def test_complete_and_budget_limited_partitions(self):
        model=Model();saved=[]
        result=pp.search(model,self.record(),self.mapping(model),4,20,4,64,saved.append)
        self.assertEqual(set(result['leaves']),{'00','01','10','11'})
        self.assertEqual(result['unresolved'],{})
        for row in saved:pp.cover.sc.partition(model,row['leaves'],row['unresolved'])
        result=pp.search(model,self.record(),self.mapping(model),4,3,4,64,lambda _:None)
        self.assertFalse(result['leaves']);self.assertEqual(len(result['unresolved']),4)

    def test_saved_accepted_labels_are_recomputed(self):
        record=self.record();record['leaves']={'':dict(witness={},proposal=-100000)};record['unresolved']={}
        model=Model();result=pp.search(model,record,self.mapping(model),2,1,4,64,lambda _:None)
        self.assertFalse(result['leaves']);self.assertEqual(set(result['unresolved']),{'0','1'})

    def test_wrong_worker_identity_comparison_and_endpoint_rejected(self):
        model=Model()
        for changed in (dict(path='111'),dict(component_sha256='wrong'),dict(upper=[-1,-100])):
            def mapping(unused,jobs):
                return [dict(pp.evaluate(model,job,ctx.prec,pp.component_digest(model.components)),**changed) for job in jobs]
            with self.assertRaises((ArithmeticError,ValueError)):pp.search(model,self.record(),mapping,1,1,4,64,lambda _:None)
        with self.assertRaises(ArithmeticError):
            pp.search(model,self.record(),lambda f,j:[],1,1,4,64,lambda _:None)

    def test_worker_preserves_exact_parent_components(self):
        prior=ctx.prec
        try:
            with patch.object(pp.cover.birth_classes,'actual',return_value=dict(bits=19,windows=32,updates=4)):
                model=pp.worker_model(Model.components,pp.cover.PARAMETERS,256)
                self.assertEqual(model.components,Model.components)
                self.assertEqual(ctx.prec,256)
                self.assertEqual(pp.component_digest(model.components),pp.component_digest(Model.components))
        finally:ctx.prec=prior

    def test_worker_rejects_precision_change(self):
        with self.assertRaises(ArithmeticError):
            pp.evaluate(Model(),('',(Q(0),Q(1)),64,None),ctx.prec+1,'digest')


if __name__=='__main__':unittest.main()
