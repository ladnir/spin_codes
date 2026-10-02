import hashlib
import json
import unittest
from copy import deepcopy

import audit_complete as audit


class ReceiptAuditTests(unittest.TestCase):
    def setUp(self):
        self.dense=dict(schema='gf16-packet-scalar-cover-1',updates=4,threshold=209715,
                        minimum_groups=49,unresolved={},leaves={p:{'witness':{}} for p in ('0','10','11')})
        self.report=dict(schema='gf16-complete-replay-1',precision=384,updates=4,threshold=209715,
                         minimum_distance=209716,output_length=1<<21,sparse_through=48,dense_through=2048,
                         dense_upper=[1,-51],sparse_upper=[1,-51],total_upper=[1,-50],
                         sparse_inner='gf-birth-classes',exact_feedback=True,single_group_exact=True)

    def check(self,report=None,dense=None,prefix=None,**options):
        report=deepcopy(self.report if report is None else report)
        raw=json.dumps(self.dense if dense is None else dense).encode()
        report['dense_sha256']=hashlib.sha256(raw).hexdigest()
        return audit.audit(report,raw,options.get('distance','1/10'),options.get('bits',49),prefix)

    def test_complete_scope_and_root_partition(self):
        self.assertEqual(self.check(),3)
        self.dense['leaves']={'':{'witness':{}}}
        self.assertEqual(self.check(),1)

    def test_hash_binding(self):
        self.report['dense_sha256']='wrong'
        with self.assertRaisesRegex(ValueError,'hash'):
            audit.audit(self.report,json.dumps(self.dense).encode(),'1/10',40)

    def test_schema_and_claim_mismatch(self):
        for key,value in (('schema','wrong'),('updates',2),('updates',True),('threshold',209714),
                          ('precision',True),('precision',64),('threshold',1<<21)):
            with self.subTest(key=key,value=value),self.assertRaises(ValueError):
                self.check(report=dict(self.report,**{key:value}))

    def test_full_occupancy_scope_required(self):
        for key,value in (('sparse_through',47),('dense_through',2047),('minimum_distance',209715),
                          ('output_length',1<<20)):
            with self.subTest(key=key),self.assertRaises(ValueError):
                self.check(report=dict(self.report,**{key:value}))
        for key,value in (('minimum_groups',0),('minimum_groups',2049),('minimum_groups',True),('minimum_groups','49'),
                          ('unresolved',{'0':{}}),('screen_only',True)):
            with self.subTest(key=key,value=value),self.assertRaises(ValueError):
                self.check(dense=dict(self.dense,**{key:value}))

    def test_dense_partition_gaps_overlaps_and_witnesses(self):
        for paths in ((),('0','10'),('','0'),('0','10','11','000'),('0','12')):
            with self.subTest(paths=paths),self.assertRaises(ValueError):
                self.check(dense=dict(self.dense,leaves={p:{'witness':{}} for p in paths}))
        with self.assertRaises(ValueError):self.check(dense=dict(self.dense,leaves={'':{}}))

    def test_sum_and_strict_margin(self):
        with self.assertRaises(ValueError):self.check(report=dict(self.report,total_upper=[1,-52]))
        with self.assertRaises(ValueError):self.check(bits=50)
        with self.assertRaises(ValueError):self.check(report=dict(self.report,dense_upper=[0,0]))
        with self.assertRaises(ValueError):self.check(report=dict(self.report,sparse_upper=[0,0]))

    def test_exact_requested_distance(self):
        with self.assertRaises(ValueError):self.check(distance='.100001')
        with self.assertRaises(ValueError):self.check(distance='209716/2097152')
        for distance,bits in (('0',40),('1',40),('.1',39),('.1',True)):
            with self.assertRaises(ValueError):self.check(distance=distance,bits=bits)

    def test_invalid_dyadic_endpoints(self):
        for pair in ((-1,-50),(1.0,-50),(True,-50),(1,False),(1,),None):
            with self.subTest(pair=pair),self.assertRaises(ValueError):audit.dyadic(pair)
        self.assertEqual(audit.dyadic([0,0]),0)

    def test_independent_prefix_exact_match_and_scope(self):
        prefix=dict(schema='gf16-fresh-sparse-prefix-1',parameters=deepcopy(self.report),
                    covered_occupancies=[1,48],upper=[1,-51])
        self.assertEqual(self.check(prefix=prefix),3)
        for key,value in (('covered_occupancies',[1,47]),('upper',[1,-52]),('schema','wrong')):
            with self.subTest(key=key),self.assertRaises(ValueError):
                self.check(prefix=dict(prefix,**{key:value}))
        for key,value in (('updates',2),('threshold',207618),('precision',256),('exact_feedback',False)):
            bad=deepcopy(prefix);bad['parameters'][key]=value
            with self.subTest(key=key),self.assertRaises(ValueError):self.check(prefix=bad)


if __name__=='__main__':unittest.main()
