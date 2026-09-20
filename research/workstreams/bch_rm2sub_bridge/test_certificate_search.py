import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import sys

import certificate_search_core as core
import search_certificates as search
import verify_certificate_search as verifier

F=core.F


class SearchTests(unittest.TestCase):
    def test_global_budget_not_eighty_bits_each(self):
        best={q:F(1,1<<43) for q in range(1,5)}
        self.assertEqual(core.summarize(best,4,40)['status'],'CERTIFIED')
        self.assertEqual(core.summarize({1:F(1,1<<100)},4,40)['status'],'UNRESOLVED')
        self.assertEqual(core.summarize({q:F(1,1<<40) for q in range(1,5)},4,40)['status'],'UNRESOLVED')

    def test_allocation(self):
        best={1:F(1,1<<50)};a=core.allocation(best,8192,40)
        self.assertLessEqual(sum(best.values())+8191*a,F(1,1<<40))
        self.assertGreater(a,F(1,1<<80))
        self.assertGreater(sum(best.values())+8191*2*a,F(1,1<<40))
        self.assertEqual(core.allocation({1:F(1)},2,40),0)

    def test_bounds_and_intervals(self):
        row=dict(occupation=1,upper=core.base.encode(F(1,8)))
        with self.assertRaises(ValueError):core.merge_rows({},[row,row],3)
        self.assertEqual(core.intervals([1,2,5,7,8]),[[1,2],[5,5],[7,8]])

    def test_proposal_skips_tried_and_probes_endpoint(self):
        spec=dict(rows=8192,message_exponent=20)
        best={1:F(1,1<<50)};task=search.propose(spec,best,set(),[])
        self.assertEqual(task['interval'],[2,9])
        following=search.propose(spec,best,{search.task_key(task)},[])
        self.assertEqual(following['interval'],[8192,8192])
        self.assertNotEqual(search.task_key(task),search.task_key(following))

    def test_timeout_reaps_child(self):
        with tempfile.TemporaryDirectory() as tmp:
            status=search.worker([sys.executable,'-c','import time; time.sleep(10)'],Path(tmp)/'log',.05)
            self.assertEqual(status,'TIMEOUT')

    def test_changed_outer_input_rejected(self):
        row=dict(path='bch_spectrum_work/nonexistent',sha256='0'*64)
        with patch.object(core.base,'read',return_value={'files':[row]}):
            with self.assertRaises(FileNotFoundError):core.outer_dependencies()

    def test_hash_authentication(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'input.json';core.base.write_new(path,{'value':1})
            data={'source_sha256':{'input.json':core.base.sha(path)}}
            core.authenticate(data,Path(tmp))
            data['source_sha256']['input.json']='0'*64
            with self.assertRaises(ValueError):core.authenticate(data,Path(tmp))

    def test_cross_instance_receipt_rejected(self):
        data={'task':{'instance':{'configuration':'wrong'}},'source_sha256':{}}
        with patch.object(core.base,'read',side_effect=[data,{'source_sha256':{}}]):
            with self.assertRaisesRegex(ValueError,'Cross-instance'):
                verifier.load_verified(Path('unused'),Path('unused'),{'configuration':'right'})

    def test_report_rejects_false_closure(self):
        spec=dict(configuration='test',message_exponent=20,rows=4)
        report=dict(instance=spec,reused_legacy_certificate=False,accepted_certificates=[],target_bits=40,
            coverage=core.summarize({},4,40))
        with tempfile.TemporaryDirectory() as tmp,patch.object(core,'instance',return_value=spec):
            path=Path(tmp)/'report.json';core.base.write_new(path,report)
            self.assertEqual(verifier.verify_report(path)['coverage']['status'],'UNRESOLVED')
            report['coverage']['status']='CERTIFIED'
            other=Path(tmp)/'bad.json';core.base.write_new(other,report)
            with self.assertRaises(ValueError):verifier.verify_report(other)

    def test_resume_no_jobs_preserves_coverage(self):
        spec=dict(configuration='test',message_exponent=20,rows=4)
        with tempfile.TemporaryDirectory() as tmp,patch.object(core,'instance',return_value=spec),\
                patch.object(search,'prior_anchors',return_value=[]):
            a=search.run(Path(tmp),'test',20,40,1,1,0,False)
            b=search.run(Path(tmp),'test',20,40,1,1,0,False)
            self.assertNotEqual(a,b)
            self.assertEqual(core.base.read(a)['coverage'],core.base.read(b)['coverage'])
            with self.assertRaises(ValueError):search.run(Path(tmp),'test',20,41,1,1,0,False)

    def test_resume_skips_timed_out_trial(self):
        spec=dict(configuration='test',message_exponent=20,rows=8)
        with tempfile.TemporaryDirectory() as tmp,patch.object(core,'instance',return_value=spec),\
                patch.object(search,'prior_anchors',return_value=[]),\
                patch.object(search,'worker',return_value='TIMEOUT') as worker:
            search.run(Path(tmp),'test',20,40,1,1,1,False)
            result=search.run(Path(tmp),'test',20,40,1,1,1,False)
            report=core.base.read(result)
            self.assertEqual(worker.call_count,2)
            self.assertEqual(len(report['attempts']),2)
            self.assertNotEqual(*(search.task_key(r['task']) for r in report['attempts']))
            self.assertEqual(report['coverage']['status'],'UNRESOLVED')


if __name__=='__main__':unittest.main()
