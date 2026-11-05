"""Pure and mocked orchestration tests; no coefficient bank or objective runs."""
from copy import deepcopy
from hashlib import sha256
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import hd_field128_gamma30_d40_r913_portable_refit_v24 as mode


def parent():
    return dict(active_groups=913, state_dimension=30, target_distance='2/5',
        construction=mode.adaptive.provider.gamma.construction(30), source_sha256={},
        witness=dict(z='1/4', alpha='1/8', grid=128, state_dimension=30,
                     markers=[['1', str(k+1)] for k in range(67)]),
        source=dict(path='C:/original/windows/ancestor.json', sha256='a'*64))


def result(source, *, timeout=False, empty=False):
    witness = deepcopy(parent()['witness'])
    complete = None if empty else dict(profile_pairs_evaluated=457**2,
        complete_selected_weight_screen=True, regional_logs_sha256='b'*64, log2_upper_screen=-80.)
    return dict(status='timeout' if timeout else 'fit-limit', witness=witness,
        source=dict(path=str(Path(source).resolve()), sha256=sha256(Path(source).read_bytes()).hexdigest()),
        last_complete_step=complete, pending_step={'stage':'fit'} if timeout else None,
        winning_marker_indices_at_last_complete=None if empty else list(range(67)))


class PureTests(unittest.TestCase):
    def test_normalize_only_new_links_inside_root(self):
        with TemporaryDirectory() as d:
            root=Path(d); value=dict(local=dict(path=str(root/'run'/'a.json'),sha256='a'*64),
                foreign=dict(path='C:/original/windows/ancestor.json',sha256='b'*64),
                bare=str(root/'not-a-link.json'), nested=[dict(path=str(root/'b.json'),sha256='c'*64)],
                already=dict(path='sub/relative.json',sha256='d'*64))
            original=deepcopy(value); changed=mode.normalize_new_metadata(value,root)
            self.assertEqual(changed['local']['path'],'run/a.json')
            self.assertEqual(changed['nested'][0]['path'],'b.json')
            self.assertEqual(changed['foreign'],value['foreign'])
            self.assertEqual(changed['bare'],value['bare'])
            self.assertEqual(changed['already'],value['already'])
            self.assertEqual(value,original)

    def test_parent_hash_scope_and_bytes_preserved(self):
        with TemporaryDirectory() as d, patch.object(mode,'HERE',Path(d)):
            path=Path(d)/mode.PARENT; path.parent.mkdir(); raw=mode.encode(parent()); path.write_bytes(raw)
            with patch.object(mode,'PARENT_SHA256',sha256(raw).hexdigest()):
                self.assertEqual(mode.load_parent(),(parent(),raw))
                self.assertEqual(path.read_bytes(),raw)
            with self.assertRaises(ValueError): mode.load_parent()
            for change in (dict(active_groups=914),dict(state_dimension=31),dict(target_distance='1/3')):
                value=parent(); value.update(change); raw=mode.encode(value); path.write_bytes(raw)
                with patch.object(mode,'PARENT_SHA256',sha256(raw).hexdigest()),self.assertRaises(ValueError):
                    mode.load_parent()

    def test_floor_and_exact_cumulative_bytes(self):
        with TemporaryDirectory() as d, patch.object(mode,'HERE',Path(d)):
            root=Path(d); out=root/'out'; out.mkdir(); usage=SimpleNamespace(free=mode.RESERVE)
            budget=mode.OutputBudget(out); value={'partial':True}; size=len(mode.encode(value))
            with patch.object(mode.shutil,'disk_usage',return_value=usage):
                with self.assertRaises(OSError): budget.write_new('a.json',value)
                usage.free=mode.RESERVE+size
                pin=budget.write_new('a.json',value)
                self.assertEqual(pin['path'],'out/a.json'); self.assertEqual(budget.used_bytes,size)
                self.assertEqual((out/'a.json').stat().st_size,size)
                with self.assertRaises(FileExistsError): budget.write_new('a.json',value)
                with patch.object(mode,'MAX_OUTPUT_BYTES',size),self.assertRaises(OSError):
                    budget.write_new('summary.json',value)
                self.assertFalse((out/'summary.json').exists()); budget.check()
                (out/'a.json').write_bytes(b'changed')
                with self.assertRaises(RuntimeError): budget.check()

    def test_guard_and_destination_rejections(self):
        for n in (-1,True,0.5):
            with self.assertRaises(ValueError): mode.reserve_guard(mode.HERE,n)
        with patch.object(mode.shutil,'disk_usage',return_value=SimpleNamespace(free=mode.RESERVE+mode.MAX_OUTPUT_BYTES-1)):
            with self.assertRaises(OSError): mode.reserve_guard(mode.HERE,mode.MAX_OUTPUT_BYTES)
        for name in ('../out','/tmp/out','C:/out','out/sub','out.'):
            with self.assertRaises(ValueError): mode.run(name)
        with TemporaryDirectory() as d:
            budget=mode.OutputBudget(d)
            for name in ('../a.json','sub/a.json','a.py'):
                with self.assertRaises(ValueError): budget.write_new(name,{})


class MockRunTests(unittest.TestCase):
    def exercise(self, first='complete', mismatch=False, source_change=False):
        temporary=TemporaryDirectory(); self.addCleanup(temporary.cleanup)
        root=Path(temporary.name); path=root/mode.PARENT; path.parent.mkdir()
        raw=mode.encode(parent()); path.write_bytes(raw); calls=[]
        def fake(source,weight,**options):
            calls.append((source,weight,options))
            row=result(source,timeout=first=='timeout' and len(calls)==1,
                       empty=first=='empty' and len(calls)==1)
            if mismatch and len(calls)==2: row['last_complete_step']['regional_logs_sha256']='c'*64
            if source_change: path.write_bytes(b'changed')
            return row
        patches=(patch.object(mode,'HERE',root),patch.object(mode,'PARENT_SHA256',sha256(raw).hexdigest()),
            patch.object(mode,'source_pins',return_value={mode.NAME:'d'*64}),
            patch.object(mode,'check_pins',side_effect=lambda pins,original:
                None if path.read_bytes()==original else (_ for _ in ()).throw(RuntimeError('parent changed'))),
            patch.object(mode.shutil,'disk_usage',return_value=SimpleNamespace(free=mode.RESERVE+mode.MAX_OUTPUT_BYTES+4096)),
            patch.object(mode.adaptive,'run',side_effect=fake))
        from contextlib import ExitStack
        with ExitStack() as stack:
            for p in patches: stack.enter_context(p)
            if mismatch or source_change:
                with self.assertRaises((ArithmeticError,RuntimeError)): mode.run('out')
                report=None
            else: report=mode.run('out')
        return root,path,raw,calls,report

    def test_exact_budgets_full_replay_and_relative_new_links(self):
        root,path,raw,calls,report=self.exercise()
        self.assertEqual(path.read_bytes(),raw); self.assertEqual(len(calls),2)
        self.assertEqual([c[1] for c in calls],[913,913])
        self.assertEqual([c[2]['max_fits'] for c in calls],[15,0])
        self.assertTrue(all(c[2]['max_seconds']==90 and c[2]['max_calls']==128 for c in calls))
        self.assertTrue(report['exact_replay_match']); self.assertFalse(report['directed'])
        self.assertFalse(report['whole_code_certificate'])
        for name in ('adaptive.json','replay.json'):
            saved=json.loads((root/'out'/name).read_bytes())
            self.assertFalse(Path(saved['source']['path']).is_absolute())
        seed=json.loads((root/'out'/'seed.json').read_bytes())
        self.assertEqual(seed['source']['path'],mode.PARENT)
        self.assertEqual(seed['witness'],parent()['witness'])

    def test_timeout_uses_only_complete_witness_and_stays_timeout(self):
        root,path,raw,calls,report=self.exercise(first='timeout')
        self.assertEqual(report['status'],'timeout'); self.assertEqual(len(calls),2)
        self.assertTrue(report['exact_replay_match']); self.assertFalse(report['directed'])

    def test_no_completed_step_does_not_invent_proposal_or_score(self):
        root,path,raw,calls,report=self.exercise(first='empty')
        self.assertEqual(len(calls),1); self.assertIsNone(report['complete_score'])
        self.assertIsNone(report['exact_replay_match']); self.assertFalse((root/'out'/'proposal.json').exists())

    def test_mismatch_and_source_mutation_fail_closed(self):
        root,*_=self.exercise(mismatch=True)
        self.assertFalse((root/'out'/'summary.json').exists())
        self.assertTrue((root/'out'/'failed_replay_summary.json').exists())
        root,*_=self.exercise(source_change=True)
        self.assertFalse((root/'out'/'summary.json').exists())


if __name__=='__main__': unittest.main()
