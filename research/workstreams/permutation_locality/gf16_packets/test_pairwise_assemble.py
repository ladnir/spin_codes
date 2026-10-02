import copy
import unittest
from unittest.mock import patch
from fractions import Fraction as Q
from types import SimpleNamespace
from flint import ctx
import pairwise_assemble as pa


class PairwiseAssembleTests(unittest.TestCase):
    def record(self):
        return dict(schema=pa.cover.SCHEMA,ensemble=pa.cover.ENSEMBLE,
            parameters=copy.deepcopy(pa.cover.PARAMETERS),leaves={'':dict(witness={})},unresolved={})

    def rows(self):return [dict(occupancy=q,upper=[1,-60]) for q in range(1,49)]

    def model(self):
        return SimpleNamespace(threshold=209715,q_min=49,data=dict(bits=19,windows=32,updates=4))

    def test_actual_geometry_is_fixed(self):
        pa.validate_geometry(self.model())
        with patch.object(pa.cover.sc,'G',1024):
            with self.assertRaises(ValueError):pa.validate_geometry(self.model())
        model=self.model();model.data['bits']=20
        with self.assertRaises(ValueError):pa.validate_geometry(model)

    def test_dense_scope_rejects_holes_and_other_ensembles(self):
        pa.validate_dense(self.record())
        for change in (dict(unresolved={'1':{}}),dict(ensemble='independent4-gf16-r4'),
                       dict(leaves={}),dict(leaves={'0':dict(witness={})}),
                       dict(leaves={'':dict(witness={}), '1':dict(witness={})}),
                       dict(screen_only=True),dict(leaves={'':{}})):
            record=self.record();record.update(change)
            with self.assertRaises(ValueError):pa.validate_dense(record)

    def test_sparse_scope_rejects_missing_duplicate_null_or_zero(self):
        for rows in (self.rows()[:-1],self.rows()[:-1]+[self.rows()[0]],
                     self.rows()[:-1]+[dict(occupancy=48,upper=None)],
                     self.rows()[:-1]+[dict(occupancy=48,upper=[0,-60])]):
            with self.assertRaises(ValueError):pa.sum_sparse(rows,48)
        self.assertEqual(pa.sum_sparse(list(reversed(self.rows())),48),Q(48,1<<60))

    def test_exact_sum_and_strict_margin(self):
        dense=dict(dense_complete=True,unresolved=0,upper=[1,-60])
        result=pa.combine(dense,self.rows(),48,40)
        self.assertEqual(pa.dyadic(result['total_upper']),Q(49,1<<60))
        with self.assertRaises(ValueError):pa.combine(dict(dense,unresolved=1),self.rows(),48,40)
        with self.assertRaises(ValueError):pa.combine(dict(dense,dense_complete=False),self.rows(),48,40)
        dense=dict(dense,upper=pa.exact_endpoint(Q(2)**-40-Q(48,1<<60)))
        with self.assertRaises(ValueError):pa.combine(dense,self.rows(),48,40)
        with self.assertRaises(ValueError):pa.exact_endpoint(Q(1,3))

    def test_verification_regenerates_sparse_instead_of_importing_labels(self):
        record=self.record();record['complete']=True;record['total_upper']=[1,-999999]
        dense=dict(dense_complete=True,unresolved=0,upper=[1,-60])
        prior=ctx.prec;ctx.prec=384
        try:
            with patch.object(pa.cover,'build_model',return_value=self.model()) as build,patch.object(pa.cover,'replay',return_value=dense) as replay,patch.object(pa.sparse,'run',return_value=self.rows()) as fresh:
                result=pa.verify(record)
                build.assert_called_once_with(384,record['parameters'])
                self.assertEqual(fresh.call_args.args[0],list(range(1,49)))
                self.assertTrue(fresh.call_args.kwargs['refined_counts'])
                self.assertTrue(result['complete'])
                self.assertEqual(result['minimum_distance'],209716)
                self.assertEqual(result['dense_occupancies'],[49,2048])
                self.assertEqual(pa.dyadic(result['total_upper']),Q(49,1<<60))
            record['unresolved']={'1':{}}
            with patch.object(pa.cover,'build_model') as build:
                with self.assertRaises(ValueError):pa.verify(record)
                build.assert_not_called()
        finally:ctx.prec=prior

    def test_failed_dense_replay_cannot_be_rescued_by_saved_sparse_numbers(self):
        record=self.record();record['sparse_upper']=[1,-10000]
        prior=ctx.prec;ctx.prec=384
        try:
            with patch.object(pa.cover,'build_model',return_value=self.model()),patch.object(pa.cover,'replay',return_value=dict(
                    dense_complete=False,unresolved=0,upper=[1,-39])),patch.object(pa.sparse,'run') as fresh:
                with self.assertRaises(ValueError):pa.verify(record)
                fresh.assert_not_called()
        finally:ctx.prec=prior

    def test_precision_change_is_rejected(self):
        prior=ctx.prec
        try:
            def replay(*args):
                ctx.prec=128
                return dict(dense_complete=True,unresolved=0,upper=[1,-60])
            with patch.object(pa.cover,'build_model',return_value=self.model()),patch.object(pa.cover,'replay',side_effect=replay),patch.object(pa.sparse,'run') as fresh:
                with self.assertRaises(ArithmeticError):pa.verify(self.record())
                fresh.assert_not_called()
        finally:ctx.prec=prior


if __name__=='__main__':unittest.main()
