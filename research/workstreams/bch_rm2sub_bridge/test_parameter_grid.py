import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import parameter_grid as grid


def settings():return grid.policy(['t64_s16','t128_s19'],[20],40,2,.05,9,7,False)


def row(name,t,s,q1,mixed,m=20):
    return dict(configuration=name,message_exponent=m,screen_state='COMPLETED',step_bits=t,state_bits=s,
        epoch_updates=2097152//t,output_bits=2097152,q1_margin_bits=q1,mixed_min_margin_bits=mixed)


class GridTests(unittest.TestCase):
    def test_weak_reference_is_not_rejected(self):
        result=grid.classify([row('t64_s20',64,20,50,-85000),row('t128_s19',128,19,50,-180000)],settings())
        self.assertEqual(result[0]['priority'],'REFERENCE')
        self.assertEqual(result[1]['priority'],'SHORTLIST_HEURISTIC')
        self.assertNotIn('full_margin_bits',result[1])

    def test_watch_not_failure_and_cost_first_shortlist(self):
        rows=[row('t64_s20',64,20,50,-85000),row('t128_s19',128,19,50,-180000),
            row('t64_s16',64,16,50,-115000),row('t256_s14',256,14,50,-600000)]
        result=grid.classify(rows,settings())
        self.assertEqual(result[1]['rank'],1);self.assertEqual(result[2]['rank'],2)
        self.assertEqual(result[3]['priority'],'WATCH_HEURISTIC')
        chosen=grid.select_refinements(result,1,1)
        self.assertEqual([r['configuration'] for r in chosen],['t128_s19','t256_s14'])

    def test_missing_reference_is_explicit(self):
        result=grid.classify([row('t64_s16',64,16,50,-115000)],settings())
        self.assertEqual(result[0]['priority'],'NO_REFERENCE');self.assertIsNone(result[0]['rank'])

    def test_q1_headroom_is_a_policy_not_a_full_margin(self):
        rows=[row('t64_s20',64,20,50,-85000),row('t64_s16',64,16,41,-85000)]
        result=grid.classify(rows,settings())
        self.assertEqual(result[1]['priority'],'WATCH_HEURISTIC')

    def test_unknown_configuration_rejected(self):
        with self.assertRaises(ValueError):grid.policy(['t64_s17'],[20],40,2,.05,9,7,False)

    def test_task_keys_bind_instance_and_stage(self):
        self.assertNotEqual(grid.key({'stage':'screen','instance':1}),grid.key({'stage':'screen','instance':2}))

    def test_resume_and_job_cap(self):
        def instance(name,m):
            t,s=map(int,name.replace('t','').replace('_s',' ').split())
            return dict(configuration=name,message_exponent=m,step_bits=t,state_bits=s,rows=8192,
                fingerprint=name+str(m))
        def execute(directory,task,seconds):
            folder=directory/'jobs'/grid.key(task);folder.mkdir(parents=True)
            grid.base.write_new(folder/'task.json',task)
            grid.base.write_new(folder/'outcome.json',dict(status='COMPLETED'))
            grid.base.write_new(folder/'result.json',dict(schema='parameter-grid-job-v1',task=task,source_sha256={},
                result=dict(q1_margin_bits=50,mixed_min_margin_bits=-85000,full_margin_bits=None)))
        with tempfile.TemporaryDirectory() as tmp,patch.object(grid.core,'instance',side_effect=instance),\
                patch.object(grid,'execute',side_effect=execute) as work:
            directory=Path(tmp)
            first=grid.run(directory,settings(),10,3,3,0,0,1)
            self.assertEqual(grid.base.read(first)['successful_screens'],1)
            second=grid.run(directory,settings(),10,3,3,0,0,10)
            self.assertEqual(grid.base.read(second)['successful_screens'],3)
            self.assertEqual(work.call_count,3)
            third=grid.run(directory,settings(),10,3,3,0,0,10)
            self.assertEqual(grid.base.read(third)['jobs_executed'],0)
            self.assertEqual(grid.verify_report(third),'GRID_REPORT_VERIFIED_NOT_A_SECURITY_CERTIFICATE')
            bad=grid.base.read(third);bad['rows'][0]['full_margin_bits']=80
            path=directory/'bad.json';grid.base.write_new(path,bad)
            with self.assertRaises(ValueError):grid.verify_report(path)
            changed=settings()|{'target_bits':41}
            with self.assertRaises(ValueError):grid.run(directory,changed,10,3,3,0,0,10)

    def test_failed_jobs_are_terminal_but_not_screens(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory=Path(tmp);task={'stage':'screen','instance':1}
            folder=directory/'jobs'/grid.key(task);folder.mkdir(parents=True)
            grid.base.write_new(folder/'task.json',task)
            grid.base.write_new(folder/'outcome.json',{'status':'TIMEOUT'})
            self.assertEqual(grid.load_job(directory,task)['state'],'TIMEOUT')


if __name__=='__main__':unittest.main()
