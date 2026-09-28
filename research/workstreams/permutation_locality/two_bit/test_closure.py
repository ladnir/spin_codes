import unittest

from flint import arb,ctx

from closure import sparse_upper,dense_scope,SPARSE_LEMMAS,LIFTED_SCHEMA
from probe import aq
from full_cover import ENSEMBLE_PREFIX
from cover import DENOMINATOR


class ClosureTests(unittest.TestCase):
    @staticmethod
    def bridge_record(q,**changes):
        return dict(schema='two-bit-support-cover-1',ensemble=[*ENSEMBLE_PREFIX,2],q=q,
                    threshold=199229,penalty='1',cutoff=8,output_degree=64,density_fold=True,
                    leaves=[dict(box=[[38,256]]*q,multiplicity=1,tilt='.1',numerators=[DENOMINATOR//2]*q)])|changes

    def test_bridge_must_fill_every_gap_with_matching_scope(self):
        record=dict(schema=LIFTED_SCHEMA,minimum_groups=409,threshold=199229,unresolved={})
        bridge=[self.bridge_record(q) for q in range(401,409)]
        self.assertEqual(dense_scope(record,bridge),199229)
        self.assertEqual(dense_scope(record,list(reversed(bridge))),199229)
        for corrupt in (bridge[:-1],bridge+[bridge[-1]],bridge[1:],
                        [self.bridge_record(400),*bridge],
                        [{**bridge[0],'threshold':199228},*bridge[1:]],
                        [{**bridge[0],'ensemble':[*ENSEMBLE_PREFIX,3]},*bridge[1:]],
                        [{**bridge[0],'penalty':'1/2'},*bridge[1:]],
                        [{**bridge[0],'leaves':[]},*bridge[1:]]):
            with self.assertRaises(ValueError):dense_scope(record,corrupt)
        # Event inclusion permits a bridge proved at a higher cutoff.
        self.assertEqual(dense_scope(record,[{**r,'threshold':209715} for r in bridge]),199229)

    def test_conservative_aggregate_closes_49_bits(self):
        ctx.prec=256
        total=sparse_upper()+arb(2)**-312
        self.assertTrue(total<(-aq('49.11')*arb(2).log()).exp())
        self.assertTrue(total>arb(2)**-50)
        with self.assertRaises(ValueError):sparse_upper(SPARSE_LEMMAS[:-1])
        with self.assertRaises(ValueError):sparse_upper(SPARSE_LEMMAS[:2]+SPARSE_LEMMAS[3:])

    def test_dense_scope_cannot_leave_a_gap_or_raise_sparse_threshold(self):
        record=dict(schema=LIFTED_SCHEMA,minimum_groups=401,threshold=104857,unresolved={})
        self.assertEqual(dense_scope(record),104857)
        for changes in ({'minimum_groups':402},{'threshold':209716},
                        {'unresolved':{'0':{}}},{'screen_only':True},{'schema':'unknown'},{'updates':3}):
            with self.assertRaises(ValueError):dense_scope({**record,**changes})


if __name__=='__main__':unittest.main()
