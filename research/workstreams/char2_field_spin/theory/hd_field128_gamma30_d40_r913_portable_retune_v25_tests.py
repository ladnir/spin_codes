"""Pure/mock tests only: no numerical objective, coefficient table, or bank."""
from contextlib import ExitStack, nullcontext
from copy import deepcopy
from hashlib import sha256
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

import hd_field128_gamma30_d40_r913_portable_retune_v25 as mode


def parent():
    return dict(active_groups=913, state_dimension=30, target_distance='2/5', construction=mode.CONSTRUCTION,
        source_sha256={}, source=dict(path='C:/unchanged/ancestor.json', sha256='a'*64),
        witness=dict(z='1/4', alpha='1/8', grid=128, state_dimension=30,
                     markers=[['1', str(i+1)] for i in range(79)]))


def normalize(value, root):
    if isinstance(value, list): return [normalize(v,root) for v in value]
    if not isinstance(value, dict): return value
    result={k:normalize(v,root) for k,v in value.items()}
    if isinstance(value.get('path'),str) and 'sha256' in value:
        try: result['path']=Path(value['path']).relative_to(root).as_posix()
        except ValueError: pass
    return result


def row(witness, score=-60.):
    return dict(active_groups=913, profile_pairs_evaluated=457**2, all_profiles_retained=True,
        complete_selected_weight_screen=True, regional_logs_sha256='b'*64, log2_upper_screen=score,
        witness=deepcopy(witness), evaluation_index=0)


class PureTests(unittest.TestCase):
    def test_parent_exact_positive_scope_and_hash(self):
        with TemporaryDirectory() as directory, patch.object(mode,'HERE',Path(directory)):
            path=Path(directory)/mode.PARENT; path.parent.mkdir(parents=True)
            for changed in (None,'weight','marker','domain'):
                value=parent()
                if changed=='weight': value['active_groups']=914
                if changed=='marker': value['witness']['markers'].pop()
                if changed=='domain': value['witness']['alpha']='2'
                raw=mode.encode(value); path.write_bytes(raw)
                with patch.object(mode,'PARENT_SHA256',sha256(raw).hexdigest()):
                    if changed is None: self.assertEqual(mode.load_parent(),(value,raw))
                    else:
                        with self.assertRaises(ValueError): mode.load_parent()
            with self.assertRaises(ValueError): mode.load_parent()

    def test_disk_output_budget_relative_metadata_and_immutability(self):
        with TemporaryDirectory() as directory, patch.object(mode,'HERE',Path(directory)):
            root=Path(directory); out=root/'out'; out.mkdir(); outputs=mode.Outputs(out,normalize)
            value=dict(local=dict(path=str(root/'a.json'),sha256='b'*64), foreign=parent()['source'])
            expected=mode.encode(normalize(value,root)); usage=SimpleNamespace(free=mode.RESERVE)
            with patch.object(mode.shutil,'disk_usage',return_value=usage):
                with self.assertRaises(OSError): outputs.write('one.json',value)
                usage.free=mode.RESERVE+len(expected)
                link=outputs.write('one.json',value)
                self.assertEqual(link['path'],'out/one.json'); self.assertEqual(outputs.used,len(expected))
                self.assertEqual((out/'one.json').read_bytes(),expected)
                self.assertEqual(value['foreign'],parent()['source'])
                with self.assertRaises(FileExistsError): outputs.write('one.json',value)
                with patch.object(mode,'MAX_OUTPUT_BYTES',outputs.used),self.assertRaises(OSError):
                    outputs.write('summary.json',{})
                outputs.check(); (out/'one.json').write_bytes(b'changed')
                with self.assertRaises(RuntimeError): outputs.check()
            for name in ('../a.json','a.py','sub/a.json'):
                with self.assertRaises(ValueError): outputs.write(name,{})

    def test_linux_memory_and_hard_wall_limits_are_installed_before_yield(self):
        resource=SimpleNamespace(RLIMIT_AS=9,RLIM_INFINITY=-1,getrlimit=Mock(return_value=(-1,-1)),setrlimit=Mock())
        with patch.object(mode.sys,'platform','linux'),patch.dict('sys.modules',{'resource':resource}),\
                patch.object(mode.signal,'SIGALRM',14,create=True),patch.object(mode.signal,'ITIMER_REAL',0,create=True),\
                patch.object(mode.signal,'getsignal',return_value='old') as get,\
                patch.object(mode.signal,'signal') as signal,patch.object(mode.signal,'setitimer',create=True) as timer:
            with mode.remote_limits():
                resource.setrlimit.assert_called_once_with(9,(4<<30,4<<30))
                timer.assert_called_once_with(0,180)
                handler=signal.call_args[0][1]
                with patch.object(mode.os,'_exit') as exit: handler(14,None); exit.assert_called_once_with(124)
            self.assertEqual(timer.call_args[0],(0,0)); self.assertEqual(signal.call_args[0],(14,'old'))
        with patch.object(mode.sys,'platform','win32'),self.assertRaises(RuntimeError):
            with mode.remote_limits(): pass

    def test_full_profile_and_guard_validation(self):
        valid=row(parent()['witness']); self.assertIs(mode.full(valid),valid)
        for key,value in (('active_groups',912),('profile_pairs_evaluated',208848),
                          ('all_profiles_retained',False),('complete_selected_weight_screen',False)):
            bad=deepcopy(valid); bad[key]=value
            with self.assertRaises(ArithmeticError): mode.full(bad)
        for value in (-1,True,.5):
            with self.assertRaises(ValueError): mode.reserve_guard(mode.HERE,value)
        for path in ('../out','/tmp/out','out/sub','out.'):
            with self.assertRaises(ValueError): mode.run(path)


class MockedRuns(unittest.TestCase):
    def exercise(self, score=-60., failure=None, pending=False):
        temporary=TemporaryDirectory(); self.addCleanup(temporary.cleanup)
        root=Path(temporary.name); path=root/mode.PARENT; path.parent.mkdir(parents=True)
        source=parent(); raw=mode.encode(source); path.write_bytes(raw)
        calls=[]; adaptive_calls=[]
        def fake_search(source_path,weight,**options):
            loaded=json.loads(Path(source_path).read_bytes())
            self.assertEqual(loaded['witness'],source['witness'])
            calls.append((weight,options))
            candidate=row(loaded['witness'],score); candidate['witness']['z']='1/8'
            options['checkpoint'](candidate)
            if failure=='memory': raise MemoryError('mock bounded allocation')
            if failure=='no-best': return dict(best=None,history=[])
            return dict(best=candidate,history=[dict(evaluation_index=0)],optimizer=dict(success=False))
        def fake_adaptive(source_path,weight,**options):
            adaptive_calls.append((weight,options))
            witness=json.loads(Path(source_path).read_bytes())['witness']
            base=row(witness,score)
            result=dict(status='fit-limit',history=[base],last_complete_step=base,witness=witness,
                winning_marker_indices_at_last_complete=list(range(len(witness['markers']))),pending_step=None)
            if len(adaptive_calls)==1:
                if failure=='baseline': base['regional_logs_sha256']='c'*64
                if failure=='literals': witness['alpha']='1/4'
                if pending: result['status']='timeout'; result['pending_step']={'fit':'unscored'}
            elif failure=='replay': base['regional_logs_sha256']='d'*64
            if failure=='source': path.write_bytes(b'mutated')
            return result
        helpers=(SimpleNamespace(normalize_new_metadata=normalize,runtime_metadata=lambda:{'test':True},
            adaptive=SimpleNamespace(run=fake_adaptive)),SimpleNamespace(run_search=fake_search))
        def check(unused,original):
            if path.read_bytes()!=original: raise RuntimeError('parent changed')
        with ExitStack() as stack:
            for replacement in (patch.object(mode,'HERE',root),patch.object(mode,'PARENT_SHA256',sha256(raw).hexdigest()),
                patch.object(mode,'helpers',return_value=helpers),patch.object(mode,'source_pins',return_value={mode.NAME:'e'*64}),
                patch.object(mode,'check_pins',side_effect=check),patch.object(mode,'remote_limits',return_value=nullcontext()),
                patch.object(mode.shutil,'disk_usage',return_value=SimpleNamespace(free=mode.RESERVE+mode.MAX_OUTPUT_BYTES+4096))):
                stack.enter_context(replacement)
            if failure=='deadline': stack.enter_context(patch.object(mode,'monotonic',side_effect=[0.,180.]))
            if failure=='source':
                with self.assertRaises(RuntimeError): mode.run('out')
                result=None
            else: result=mode.run('out')
        return root,path,raw,calls,adaptive_calls,result

    def test_retune_then_optional_fit_and_fresh_subset_replay(self):
        root,path,raw,calls,adaptive,result=self.exercise()
        self.assertEqual(path.read_bytes(),raw)
        self.assertEqual(calls[0][0],913); self.assertEqual(calls[0][1]['max_evaluations'],32)
        self.assertLessEqual(calls[0][1]['max_seconds'],90)
        self.assertEqual([c[0] for c in adaptive],[913,913])
        self.assertEqual([c[1]['max_fits'] for c in adaptive],[8,0])
        self.assertTrue(all(c[1]['max_calls']==128 for c in adaptive))
        self.assertLessEqual(adaptive[0][1]['max_seconds'],45); self.assertLessEqual(adaptive[1][1]['max_seconds'],20)
        self.assertTrue(result['exact_replay_match']); self.assertEqual(result['status'],'FLOAT-budget-exhausted')
        proposal=json.loads((root/'out/proposal.json').read_bytes())
        self.assertEqual(proposal['witness']['z'],'1/8'); self.assertEqual(proposal['witness']['markers'],parent()['witness']['markers'])
        self.assertFalse(result['directed']); self.assertFalse(result['whole_code_certificate'])

    def test_target_skips_fits_but_not_two_fresh_replays(self):
        *_,calls,result=self.exercise(score=-90.)
        self.assertEqual([c[1]['max_fits'] for c in calls],[0,0])
        self.assertEqual(result['status'],'FLOAT-target-met')

    def test_pending_uses_last_complete_witness_only(self):
        *_,result=self.exercise(pending=True)
        self.assertEqual(result['adaptive_status'],'timeout'); self.assertTrue(result['exact_replay_match'])
        self.assertEqual(result['status'],'FLOAT-budget-exhausted')

    def test_fail_closed_memory_baseline_replay_and_parent_mutation(self):
        for failure in ('memory','baseline','replay','literals'):
            root,*_,result=self.exercise(failure=failure)
            self.assertEqual(result['status'],'partial-resource-or-validation-stop')
            self.assertIsNone(result['complete_score']); self.assertFalse(result['exact_replay_match'])
            self.assertFalse((root/'out/summary.json').exists())
            self.assertTrue((root/'out/evaluation_000.json').exists())
            self.assertTrue((root/'out/partial_summary.json').exists())
        root,*_=self.exercise(failure='source')
        self.assertFalse((root/'out/summary.json').exists())

    def test_no_best_does_not_invent_score_or_proposal(self):
        root,_,_,_,calls,result=self.exercise(failure='no-best')
        self.assertFalse(calls); self.assertIsNone(result['complete_score'])
        self.assertIsNone(result['exact_replay_match']); self.assertFalse((root/'out/proposal.json').exists())

    def test_exhausted_whole_job_budget_prevents_numerical_call(self):
        root,_,_,search,adaptive,result=self.exercise(failure='deadline')
        self.assertFalse(search); self.assertFalse(adaptive)
        self.assertEqual(result['error_type'],'TimeoutError')
        self.assertFalse((root/'out/summary.json').exists())


if __name__=='__main__': unittest.main()
