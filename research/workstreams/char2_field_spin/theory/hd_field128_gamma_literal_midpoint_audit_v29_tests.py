from contextlib import ExitStack
from copy import deepcopy
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

import hd_field128_gamma_literal_midpoint_audit_v29 as mode
from hd_field128_gamma_portable_midpoint_directed_v29_tests import complete, encoded as encode, seed


class Fixture:
    def __init__(self,weight=914,target=54,bound='-59'):
        self.root=mode.HERE/'__literal_midpoint_v29_memory__'
        self.directory=self.root/'run';self.files={};self.loaded=[]
        self.seed=seed()
        self.put('seed.json',self.seed);self.seed_hash=self.hash('seed.json')
        self.link=dict(path='seed.json',sha256=self.seed_hash)
        self.execution={mode.runner.NAME:mode.RUNNER_SHA,'p.py':'1'*64,'b.py':'2'*64}
        self.pins=dict(self.execution,**{mode.NAME:'3'*64})
        first,last,rows,profiles=mode.runner.geometry(weight)
        witness=dict(self.seed['witness'],inner_width=64,gamma_extension_degree=16,gamma_base_field_order=256,
            joint_rank_return_cap=True,all_pair_before_alpha=True,precision=256,packet_mantissa_bits=80)
        self.bank=dict(mode.frozen.previous.EXPECTED,witness=witness,directed=False,source_sha256={'p.py':'1'*64},
            bank_source_sha256={'b.py':'2'*64},witness_source=self.link,runner_source_sha256=self.execution,
            runtime={'OPENBLAS_NUM_THREADS':'1'},fresh_bank_rebuilt=True,construction_unchanged=True)
        self.driver=mode.runner.driver_record(self.execution,self.link,self.seed,weight,4<<30,target)
        self.record=dict(active_groups=weight,requested_failure_bits=target,max_seconds=500.,
            count_provenance={'maximum_refined_singletons':64},witness_source=self.link)
        self.record.update(self.context())
        self.record.update(complete(weight))
        self.record.update(requested_failure_bits=target,source_pins_unchanged_at_exit=True,
            required_normal_process_exit=True,
            meets_requested_selected_weight_target=F(bound)<=-target,wrapper_admitted=F(bound)<=-target,
            literal_seed_only=True,historical_origin_verified=False,calibration_history_replayed=False,
            source_calibration_weight=913,total_job_seconds=510.,log2_upper_rational=bound,
            execution_driver=self.driver['driver'],execution_source_sha256=self.execution,
            rows=[dict(even_full=a,odd_full_range=[first,last],profile_pairs_evaluated=rows)
                  for a in range(first,last+1)])
        self.sync()

    def context(self):
        with patch.object(mode.frozen.producer,'source_pins',return_value={'p.py':'1'*64}):
            return mode.expected_context(self.record,self.bank)

    def put(self,name,value):self.files[self.root/name]=encode(value)
    def hash(self,name):return sha256(self.files[self.root/name]).hexdigest()

    def sync(self):
        context=self.context()
        digest=sha256(json.dumps(context,sort_keys=True,separators=(',',':')).encode()).hexdigest()
        self.record['context_sha256']=digest
        self.put('run/driver.json',self.driver);self.put('run/bank.json',self.bank)
        self.put('run/context.json',dict(context=context,context_sha256=digest,complete_selected_weight=False,
            covers_all_supports_at_weight=False,whole_code_certificate=False,scope='metadata only; no numerical bound'))
        self.put('run/summary.json',self.record)
        for row in self.record['rows']:
            self.put(f"run/row_{row['even_full']:04d}.json",dict(row,context_sha256=digest,
                active_groups=self.record['active_groups'],directed=True,whole_code_certificate=False,
                scope='completed even-profile row only; requires midpoint-count metadata context'))

    def audit(self,change_during_core=None):
        fixture=self
        class Loader:
            def __init__(self,distance):self.files={}
            def load(self,path,*,expected=None,role='input'):
                path=Path(path);fixture.loaded.append(path)
                if path not in fixture.files:raise FileNotFoundError(path)
                raw=fixture.files[path];digest=sha256(raw).hexdigest()
                mode.require(expected is None or digest==expected,'Input byte pin mismatch')
                pin=dict(path=str(path),sha256=digest,bytes=len(raw))
                self.files[str(path)]=dict(pin,roles=[role])
                return json.loads(raw),pin
        portable=mode.runner.baseline().pinned_module(mode.runner.PORTABLE,mode.runner.PORTABLE_SHA)
        prior=SimpleNamespace(validation=lambda:SimpleNamespace(digest=portable.digest_string),
            pinned_module=lambda *args:portable)
        oldstat,oldglob=Path.stat,Path.glob
        def stat(path,*args,**kw):
            if path in fixture.files:return SimpleNamespace(st_size=len(fixture.files[path]))
            return oldstat(path,*args,**kw)
        def glob(path,pattern,**kw):
            if path==fixture.directory:
                return (p for p in fixture.files if p.parent==path and p.match(pattern))
            return oldglob(path,pattern,**kw)
        def links(pins):
            for pin in pins:
                mode.require(sha256(fixture.files[Path(pin['path'])]).hexdigest()==pin['sha256'],'Changed input at exit')
        def numeric(record,*,pairs,precision):
            self.core_calls.append((deepcopy(record),pairs,precision))
            if change_during_core:change_during_core(self)
            return dict(active_groups=record['active_groups'],even_rows=record['expected_even_rows'],
                profile_pairs=record['expected_profile_pairs'],saved_log2_upper_rational=record['log2_upper_rational'],
                saved_row_reaggregation_log2_upper_rational=record['log2_upper_rational'],
                all_producer_midpoint_and_final_counts_match=True,saved_ceiling_dominates_row_enclosures=True)
        self.core_calls=[]
        with ExitStack() as stack:
            for obj,name,value in ((mode,'HERE',self.root),(mode,'source_pins',lambda:self.pins.copy()),
                (mode.runner,'baseline',lambda:prior),(mode.runner,'source_pins',lambda:self.execution.copy()),
                (mode.frozen.producer,'source_pins',lambda:{'p.py':'1'*64}),
                (mode.frozen.producer.provider,'source_pins',lambda:{'b.py':'2'*64}),
                (mode.frozen.producer,'check_source_pins',Mock()),(mode.frozen.index,'Validator',Loader),
                (mode.frozen,'audit_record',numeric),(mode.frozen.index,'union_ceiling',lambda values,prec:values[0]),
                (mode.frozen.previous.band,'_check_links',links),(Path,'stat',stat),(Path,'glob',glob)):
                stack.enter_context(patch.object(obj,name,value))
            return mode.audit_summary('run/summary.json',expected_summary_sha256=self.hash('run/summary.json'),
                expected_seed_sha256=self.seed_hash)


class AuditTests(unittest.TestCase):
    def test_source_identity_and_path_boundaries(self):
        for name in ('../x.json','C:/x.json','//server/x.json','/tmp/x.json','a\\b','a//b','./a'):
            with self.assertRaises(ValueError):mode.portable_path(name)
        with patch.object(Path,'resolve',lambda p:mode.HERE if p==mode.HERE else mode.HERE.parent/'outside.json'):
            with self.assertRaises(ValueError):mode.portable_path('symlink.json')

    def test_complete_914_and_later1145_with_literal_only_core_dispatch(self):
        for weight,target,bound in ((914,54,'-59'),(1145,48,'-49'),(8191,40,'-42')):
            fixture=Fixture(weight,target,bound)
            report=fixture.audit()
            self.assertEqual(report['requested_failure_bits'],target)
            self.assertEqual(fixture.core_calls[0][1:],(4096,384))
            self.assertFalse(report['historical_origin_verified'])

    def test_complete_target_miss_is_audited_not_relabelled_incomplete(self):
        fixture=Fixture(1145,48,'-47');result=fixture.audit()
        self.assertEqual(result['audit_result'],'PASS');self.assertTrue(result['complete_selected_weight'])
        self.assertFalse(result['meets_requested_selected_weight_target'])
        self.assertEqual(result['saved_log2_upper_rational'],'-47')

    def test_partial_and_false_admission_rejected_before_core(self):
        for key,value in (('status','incomplete'),('failure',{}),('completed_even_rows',457),
            ('profile_pairs_evaluated',1),('historical_origin_verified',True),
            ('requested_failure_bits',39),('meets_requested_selected_weight_target',False),('wrapper_admitted',False)):
            fixture=Fixture();fixture.record[key]=value;fixture.put('run/summary.json',fixture.record)
            with self.subTest(key=key),self.assertRaises(ValueError):fixture.audit()
            self.assertEqual(fixture.core_calls,[])

    def test_literal_bank_context_driver_and_checkpoint_mutations_fail(self):
        changes=(('run/bank.json',lambda v:v['witness']['markers'][0].__setitem__(0,'2')),
            ('run/driver.json',lambda v:v.__setitem__('requested_failure_bits',48)),
            ('run/context.json',lambda v:v['context'].__setitem__('active_groups',915)),
            ('run/row_0000.json',lambda v:v.__setitem__('active_groups',915)))
        for name,change in changes:
            fixture=Fixture();value=json.loads(fixture.files[fixture.root/name]);change(value);fixture.put(name,value)
            with self.subTest(name=name),self.assertRaises(ValueError):fixture.audit()
            self.assertEqual(fixture.core_calls,[])
        fixture=Fixture();fixture.put('run/row_9999.json',{})
        with self.assertRaises(ValueError):fixture.audit()
        fixture=Fixture();del fixture.files[fixture.directory/'row_0000.json']
        with self.assertRaises(FileNotFoundError):fixture.audit()

    def test_exit_byte_checks_and_unknown_fresh_source(self):
        fixture=Fixture()
        with self.assertRaisesRegex(ValueError,'Changed input'):
            fixture.audit(lambda f:f.files.__setitem__(f.directory/'row_0000.json',b'changed'))
        with self.assertRaises(ValueError):mode.fresh_source_maps({'source_sha256':{'evil.py':'x'}},{})

    def test_main_enters_reviewed_hard_limits_before_any_inputs(self):
        with patch.object(mode.argparse.ArgumentParser,'parse_args',return_value=None),\
                patch.object(mode.limits,'remote_limits',side_effect=RuntimeError('hard limits')):
            with self.assertRaisesRegex(RuntimeError,'hard limits'):mode.main()

    def test_audit_output_guard_without_writes(self):
        for free,value in ((20<<30,{}),(21<<30,{'x':'x'*(256<<10)})):
            with patch.object(mode.shutil,'disk_usage',return_value=SimpleNamespace(free=free)),\
                    patch.object(mode.frozen.producer,'write_new') as writer:
                with self.assertRaises(ValueError):mode.write_remote('output.json',value)
                writer.assert_not_called()


if __name__=='__main__':unittest.main()
