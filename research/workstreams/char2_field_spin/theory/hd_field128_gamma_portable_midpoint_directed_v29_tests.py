"""Memory-backed orchestration tests; no bank, coefficient or objective runs."""
from contextlib import ExitStack, nullcontext
from copy import deepcopy
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

import hd_field128_gamma_portable_midpoint_directed_v29 as mode
import hd_field128_gamma30_d40_r913_midpoint_directed_v25 as previous
import hd_field128_gamma30_d40_r913_directed_v25 as limits


def encoded(value):
    return (json.dumps(value,indent=2,allow_nan=False)+'\n').encode('utf-8')


class MemoryFiles:
    def __init__(self):
        self.root = mode.HERE/'__midpoint_v29_memory__'
        self.files, self.dirs = {}, {self.root}
        self.old_read, self.old_exists = Path.read_bytes, Path.exists

    def inside(self,path):
        return path == self.root or self.root in path.parents

    def read(self,path):
        if self.inside(path):
            if path not in self.files: raise FileNotFoundError(path)
            return self.files[path]
        return self.old_read(path)

    def exists(self,path):
        return path in self.files or path in self.dirs if self.inside(path) else self.old_exists(path)

    def mkdir(self,path,*args,**kwargs):
        if not self.inside(path): raise AssertionError('Unexpected real mkdir')
        if path in self.dirs and not kwargs.get('exist_ok'): raise FileExistsError(path)
        self.dirs.add(path)

    def write(self,path,value):
        path = Path(path)
        if not self.inside(path): raise AssertionError('Unexpected real write')
        if path in self.files: raise FileExistsError(path)
        self.files[path] = encoded(value)

    def __enter__(self):
        self.stack=ExitStack()
        self.stack.enter_context(patch.object(Path,'read_bytes',lambda p:self.read(p)))
        self.stack.enter_context(patch.object(Path,'exists',lambda p:self.exists(p)))
        self.stack.enter_context(patch.object(Path,'mkdir',lambda p,*a,**k:self.mkdir(p,*a,**k)))
        return self

    def __exit__(self,*args): return self.stack.__exit__(*args)


def seed():
    return dict(schema='gamma-portable-exact-witness-v1', active_groups=913,
        field_order=1<<128,state_dimension=30,target_distance='2/5',construction=mode.CONSTRUCTION,
        directed=False,whole_code_certificate=False,
        source={'path':'C:/unavailable/history.json','sha256':'f'*64},
        source_sha256={'missing_optimizer.py':'f'*64},
        witness=dict(z='1/4',alpha='1/8',grid=128,state_dimension=30,
            markers=[['1',str(k+1)] for k in range(83)]))


def complete(weight):
    first,last,rows,profiles=mode.geometry(weight)
    return dict(status='complete',failure=None,active_groups=weight,directed=True,
        complete_selected_weight=True,covers_all_supports_at_weight=True,covers_all_even_odd_splits=True,
        whole_code_certificate=False,completed_even_rows=rows,expected_even_rows=rows,
        completed_even_full_range=[first,last],pending_even_full_range=None,
        profile_pairs_evaluated=profiles,expected_profile_pairs=profiles,
        requested_failure_bits=54,meets_requested_selected_weight_target=True,log2_upper_rational='-59')


class PureTests(unittest.TestCase):
    def test_geometry_original_weight_dense_and_terminal_exclusion(self):
        expected={914:(0,457,458,209764),915:(0,457,458,209764),
            919:(0,459,460,211600),1145:(0,572,573,328329),1149:(0,574,575,330625),
            4097:(1,2048,2048,4194304),8191:(4095,4095,1,1)}
        for weight,geometry in expected.items(): self.assertEqual(mode.geometry(weight),geometry)
        for weight in (True,0,-1,8192,914.0,'914'):
            with self.assertRaises(ValueError): mode.geometry(weight)

    def test_admission_exact_coverage_and_actual_54_bit_ceiling(self):
        for weight in (914,919,1145,1149,8191): self.assertTrue(mode.admitted(complete(weight),weight))
        for key,value in [('status','incomplete'),('failure',{}),('active_groups',913),
            ('directed',1),('whole_code_certificate',True),('covers_all_supports_at_weight',False),
            ('covers_all_even_odd_splits',False),('complete_selected_weight',False),
            ('completed_even_rows',457),('expected_even_rows',457),
            ('completed_even_full_range',[1,457]),('pending_even_full_range',[457,457]),
            ('profile_pairs_evaluated',208849),('expected_profile_pairs',208849),
            ('requested_failure_bits',40),('meets_requested_selected_weight_target',False),
            ('log2_upper_rational','-53'),('log2_upper_rational',None),('log2_upper_rational','NaN')]:
            result=complete(914);result[key]=value
            self.assertFalse(mode.admitted(result,914),key)

    def test_explicit_relaxed_target_and_guards(self):
        result=complete(914);result.update(requested_failure_bits=48,log2_upper_rational='-49')
        self.assertTrue(mode.admitted(result,914,48));self.assertFalse(mode.admitted(result,914))
        self.assertEqual(result['status'],'complete')
        for value in (True,39,4097,48.0,'48'):
            with self.assertRaises(ValueError):mode.target_bits(value)
        for key in complete(914):
            result=complete(914);del result[key]
            self.assertFalse(mode.admitted(result,914),key)

    def test_frozen_identity_source_closure_and_no_seed_ancestry(self):
        self.assertIs(mode.baseline(),previous)
        with patch.object(mode,'BASE_SHA','0'*64),patch.object(mode.importlib,'import_module') as imp:
            with self.assertRaises(RuntimeError):mode.baseline()
            imp.assert_not_called()
        prior=SimpleNamespace(source_pins=Mock(return_value={'old.py':'a'*64}))
        portable=SimpleNamespace(source_pins=lambda:{'new.py':'b'*64})
        result=mode.source_pins(prior,portable,object())
        self.assertEqual(result['old.py'],'a'*64)
        self.assertEqual(result[mode.NAME],sha256((mode.HERE/mode.NAME).read_bytes()).hexdigest())
        self.assertEqual(set(prior.source_pins.call_args.kwargs),{'producer'})
        self.assertNotIn('missing_optimizer.py',result)
        prior.source_pins.return_value={'new.py':'c'*64}
        with self.assertRaisesRegex(RuntimeError,'conflict'):mode.source_pins(prior,portable,object())

    def test_canonical_pinned_seed_loader_does_not_open_history(self):
        portable=previous.pinned_module(mode.PORTABLE,mode.PORTABLE_SHA)
        with MemoryFiles() as mem,patch.object(portable,'HERE',mem.root):
            proposal=seed();mem.files[mem.root/'seed.json']=encoded(proposal)
            digest=sha256(mem.files[mem.root/'seed.json']).hexdigest()
            got,parameters,link=portable.load_seed('seed.json',digest)
            self.assertEqual(got,proposal);self.assertEqual(parameters[0:2],(F(1,4),F(1,8)))
            self.assertEqual(len(parameters[2]),83);self.assertEqual(link,dict(path='seed.json',sha256=digest))
            with self.assertRaises(ValueError):portable.load_seed('seed.json','0'*64)
            for mutate in (lambda s:s.update(state_dimension=31),
                lambda s:s['witness'].update(z='2/8'),lambda s:s['witness'].update(alpha='2'),
                lambda s:s['witness']['markers'].append(s['witness']['markers'][0])):
                bad=deepcopy(proposal);mutate(bad)
                with self.assertRaises(ValueError):portable.validate_seed(bad)
        for name in ('../seed.json','C:/seed.json','//server/seed.json','x\\seed.json'):
            with self.assertRaises(ValueError):limits.relative_name(name,json_file=True,basename=True)

    def test_resource_budget_exact_serialization_and_failed_charge(self):
        with MemoryFiles() as mem:
            value={'tiny':True};size=len(encoded(value));guard=Mock()
            with patch.object(previous,'MAX_OUTPUT_BYTES',size):
                budget=previous.OutputBudget(mem.root,guard,mem.write)
                budget.write_new(mem.root/'one.json',value)
                self.assertEqual(budget.used_bytes,size);guard.assert_called_once_with(size)
                with self.assertRaises(FileExistsError):budget.write_new(mem.root/'one.json',value)
                with self.assertRaises(ValueError):budget.write_new(mem.root/'child'/'two.json',value)
                with self.assertRaisesRegex(RuntimeError,'2MiB'):budget.write_new(mem.root/'two.json',value)
            broken=previous.OutputBudget(mem.root,guard,Mock(side_effect=OSError('disk')))
            with self.assertRaises(OSError):broken.write_new(mem.root/'broken.json',value)
            self.assertEqual(broken.used_bytes,size)
            denied=previous.OutputBudget(mem.root,Mock(side_effect=RuntimeError('reserve')),mem.write)
            with self.assertRaises(RuntimeError):denied.write_new(mem.root/'denied.json',value)
            self.assertNotIn(mem.root/'denied.json',mem.files)
        with patch.object(limits.shutil,'disk_usage',return_value=SimpleNamespace(free=20<<30)):
            with self.assertRaises(RuntimeError):limits.reserve_guard(1)
        with patch.object(limits.shutil,'disk_usage',return_value=SimpleNamespace(free=(20<<30)+1)):
            limits.reserve_guard(1)


class MockedRunTests(unittest.TestCase):
    def exercise(self,weight=914,action='complete',failure_bits=54):
        with MemoryFiles() as mem:
            events=[];captured={};proposal=seed();first,last,rows,profiles=mode.geometry(weight)
            bank=SimpleNamespace(check_pins=Mock(),metadata=lambda:{'mocked':True})
            def limit():events.append('limits');return nullcontext(4<<30)
            validation=SimpleNamespace(relative_name=limits.relative_name,digest=limits.digest,
                reserve_guard=Mock(),linux_limits=limit,check_pins=Mock())
            portable=SimpleNamespace(load_seed=Mock(return_value=(proposal,
                (F(1,4),F(1,8),tuple((F(x),F(y)) for x,y in proposal['witness']['markers'])),
                {'path':'seed.json','sha256':'a'*64})),runtime_metadata=lambda:{'mocked':True})
            def receipt(actual,original,**options):
                self.assertIs(actual,bank);self.assertEqual(original,weight);captured.update(options)
                options['context_checkpoint']({'context':{'active_groups':weight},'context_sha256':'c'*64})
                count=min(rows,2) if action=='partial' else rows
                for a in range(first,first+count):
                    options['checkpoint']({'even_full':a,'context_sha256':'c'*64,'complete_selected_weight':False})
                if action=='exception':raise MemoryError('mock stop')
                if action=='mutation':mem.files[mem.root/'out'/f'row_{first:04d}.json']=b'changed'
                result=complete(weight);result['requested_failure_bits']=failure_bits
                if failure_bits==48:result['log2_upper_rational']='-49'
                if action=='partial':result.update(status='incomplete',failure={'error_type':'TimeoutError'},
                    complete_selected_weight=False,covers_all_supports_at_weight=False,
                    covers_all_even_odd_splits=False,completed_even_rows=count,
                    profile_pairs_evaluated=count*rows,pending_even_full_range=[first+count,last],
                    log2_upper_rational=None,meets_requested_selected_weight_target=False)
                if action=='miss':result['log2_upper_rational']='-53'
                return result
            producer=SimpleNamespace(provider=SimpleNamespace(AliasBank=Mock(return_value=bank)),
                weight_receipt=receipt,write_new=mem.write)
            def helper(name,digest):
                events.append('load:'+name)
                return portable if name==mode.PORTABLE else producer
            prior=SimpleNamespace(validation=lambda:validation,pinned_module=helper,
                PRODUCER=previous.PRODUCER,PRODUCER_SHA256=previous.PRODUCER_SHA256,
                OutputBudget=previous.OutputBudget)
            with patch.object(mode,'HERE',mem.root),patch.object(mode,'baseline',return_value=prior),\
                    patch.object(mode,'source_pins',return_value={mode.NAME:'d'*64}):
                if action in ('exception','mutation'):
                    with self.assertRaises((MemoryError,RuntimeError)):mode.run('seed.json','a'*64,weight,'out',failure_bits=failure_bits)
                    self.assertGreaterEqual(bank.check_pins.call_count,2)
                    if action=='exception':self.assertNotIn(mem.root/'out'/'summary.json',mem.files)
                    return
                result=mode.run('seed.json','a'*64,weight,'out',failure_bits=failure_bits)
            self.assertEqual(events,['limits','load:'+mode.PORTABLE,'load:'+previous.PRODUCER])
            self.assertEqual(captured['failure_bits'],failure_bits);self.assertEqual(captured['max_singletons'],64)
            self.assertTrue(0<captured['max_seconds']<=595)
            self.assertEqual(producer.provider.AliasBank.call_args.args[0:2],(F(1,4),F(1,8)))
            self.assertEqual(len(producer.provider.AliasBank.call_args.args[2]),83)
            options=producer.provider.AliasBank.call_args.kwargs
            for key,value in dict(state=30,q=1<<128,distance=F(2,5),grid=128,
                precision=256,bits=80,max_table_bytes=64<<20,outer_measure='invertible-toeplitz').items():
                self.assertEqual(options[key],value)
            validation.reserve_guard.assert_any_call(2<<20)
            self.assertLessEqual(validation.check_pins.call_count,8)
            self.assertEqual(sum(p.name.startswith('row_') for p in mem.files),min(rows,2) if action=='partial' else rows)
            driver=json.loads(mem.files[mem.root/'out'/'driver.json'])
            self.assertEqual(driver['source_calibration_weight'],913)
            self.assertEqual(driver['active_groups'],weight)
            self.assertEqual(driver['requested_failure_bits'],failure_bits)
            self.assertEqual(driver['expected_profile_pairs'],profiles)
            self.assertFalse(driver['fresh_bank_rebuilt']);self.assertFalse(driver['historical_origin_verified'])
            self.assertEqual(driver['max_output_bytes'],2<<20)
            self.assertEqual(driver['memory_address_space_bytes'],4<<30)
            self.assertTrue(json.loads(mem.files[mem.root/'out'/'bank.json'])['fresh_bank_rebuilt'])
            self.assertFalse(result['whole_code_certificate']);self.assertFalse(result['calibration_history_replayed'])
            self.assertEqual(result['wrapper_admitted'],action=='complete')
            self.assertEqual(mode.admitted(result,weight,failure_bits),action=='complete')

    def test_complete_original914_and_later1149(self):
        self.exercise(914);self.exercise(1149,failure_bits=48)

    def test_dense_original_weight_not_encoder_complemented(self):self.exercise(8191)
    def test_partial_rows_no_bound_or_admission(self):self.exercise(action='partial')
    def test_complete_but_target_missed_not_admitted(self):self.exercise(action='miss')
    def test_resource_failure_preserves_rows_and_checks_exit(self):self.exercise(action='exception')
    def test_output_mutation_rejected_at_exit(self):self.exercise(action='mutation')


if __name__=='__main__':unittest.main()
