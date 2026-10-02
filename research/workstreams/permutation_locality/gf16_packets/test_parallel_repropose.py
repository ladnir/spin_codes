from fractions import Fraction as Q
from types import SimpleNamespace
import unittest
from flint import arb,ctx
import parallel_repropose as parallel


class ParallelCoverTests(unittest.TestCase):
    def test_fresh_success_does_not_import_old_labels(self):
        ctx.prec=192
        class Model:
            def proposal(self,cell):return -100,{'fresh':True}
            def outward(self,cell,witness):
                self.checked=witness;return arb(2)**-99
        model=Model();result=parallel.check_cell(model,(Q(0),Q(1)),{'witness':{'old':True}},60)
        self.assertEqual(model.checked,{'fresh':True})
        self.assertEqual(result['upper'],[1,-99])

    def test_old_witness_is_freshly_rechecked_or_left_unresolved(self):
        ctx.prec=192
        class Model:
            def proposal(self,cell):return 1,{}
            def outward(self,cell,witness):self.checked=witness;return self.bound
        model=Model();old={'witness':{'old':True},'proposal':-1000}
        model.bound=arb(2)**-70
        self.assertEqual(parallel.check_cell(model,(Q(0),Q(1)),old,60)['upper'],[1,-70])
        self.assertEqual(model.checked,old['witness'])
        model.bound=arb(2)**-50
        self.assertIsNone(parallel.check_cell(model,(Q(0),Q(1)),old,60)['upper'])
        self.assertIsNone(parallel.check_cell(model,(Q(0),Q(1)),None,60)['upper'])

    def test_parent_rejects_missing_wrong_or_insufficient_evidence(self):
        model=SimpleNamespace(root=(Q(0),Q(1)))
        for path,row in [('1',{'upper':[1,-70],'witness':{}}),
                         ('0',{'upper':[0,-70],'witness':{}}),
                         ('0',{'upper':[1.,-70],'witness':{}}),
                         ('0',{'upper':[1,-50],'witness':{}}),
                         ('0',{'upper':[1,-70]})]:
            with self.assertRaises(ArithmeticError):
                parallel.accept_result(model,{}, {'0':{},'1':{}},'0',(path,row),60)
        leaves={};unresolved={'0':{},'1':{}}
        result=('0',dict(upper=[1,-70],witness={'fresh':True},proposal=-70))
        self.assertTrue(parallel.accept_result(model,leaves,unresolved,'0',result,60))
        self.assertEqual(set(leaves),{'0'});self.assertEqual(set(unresolved),{'1'})
        with self.assertRaises(ArithmeticError):parallel.accept_result(model,leaves,unresolved,'0',result,60)
        self.assertFalse(parallel.accept_result(model,leaves,unresolved,'1',('1',dict(upper=None,proposal=1)),60))
        cells=parallel.sc.partition(model,leaves,unresolved)
        self.assertEqual(sum(hi-lo for lo,hi in cells.values()),1)

    def test_bad_metadata_fails_before_authentication(self):
        base=dict(schema=parallel.sc.SCHEMA,kernel='birth-classes',updates=4)
        for record,precision in ((dict(base,updates=5),192),(dict(base,updates=True),192),
                (dict(base,kernel='old'),192),(dict(base,row_bias='1/2'),192),(base,127)):
            with self.assertRaises(ValueError):parallel.make_model(record,precision,geometry_only=True)


if __name__=='__main__':unittest.main()
