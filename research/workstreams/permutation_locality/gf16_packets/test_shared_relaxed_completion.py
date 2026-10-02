from concurrent.futures import FIRST_COMPLETED, Future
import copy
from fractions import Fraction as Q
import hashlib
import json
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from flint import arb, ctx
import shared_relaxed_completion as scheduler


class ControlledExecutor:
    """No threads, clocks, or proof computation: tests decide completions."""
    def __init__(self, model, choose=None, outcomes=None):
        self.model = model
        self.choose = choose or (lambda requests:min(requests, key=lambda r:r['id']))
        self.outcomes = outcomes or {}
        self.requests = {}
        self.submissions = []
        self.completions = []
        self.peak = 0

    def submit(self, function, request):
        future = Future()
        self.assert_function = function
        self.requests[future] = copy.deepcopy(request)
        self.submissions.append((request['job'][0],tuple(self.completions)))
        self.peak = max(self.peak,sum(not f.done() for f in self.requests))
        return future

    def row(self, request):
        path, cell, _, old = request['job']
        row = dict(path=path, cell=list(map(str,cell)), proposal=-60.,
                   witness={'fresh':old}, upper=[1,-60],
                   component_sha256=request['scope']['component_sha256'],
                   proof_scope=copy.deepcopy(request['scope']), request_id=request['id'])
        outcome = self.outcomes.get(path)
        if isinstance(outcome, BaseException):
            raise outcome
        if outcome is False:
            row.update(upper=None, proposal=1.)
        elif isinstance(outcome, dict):
            row.update(copy.deepcopy(outcome))
        return row

    def wait(self, futures, return_when):
        if return_when != FIRST_COMPLETED:
            raise AssertionError('a batch barrier was requested')
        request = self.choose([self.requests[f] for f in futures])
        future = next(f for f in futures if self.requests[f]['id']==request['id'])
        try:
            future.set_result(self.row(request))
        except BaseException as error:
            future.set_exception(error)
        self.completions.append(request['job'][0])
        return {future}, set(futures)-{future}


class CompletionTests(unittest.TestCase):
    def model(self):
        return SimpleNamespace(root=(Q(0),Q(1)), components=[(Q(1),Q(0),0),(Q(100),Q(1,2),1)],
            threshold=146800, q_min=33, tilt=Q(3,16), variance_shuffle=True, variance_bins=16,
            regional_count=True, inner=SimpleNamespace(__name__='birth_classes'),
            data=dict(bits=19,windows=32,updates=2))

    def run_search(self, executor=None, cover=None, **changes):
        model = self.model()
        executor = executor or ControlledExecutor(model)
        snapshots = []
        options = dict(workers=2,max_cells=16,max_depth=8,target_bits=32,precision=384,
                       checkpoint=snapshots.append,reuse_parent_witness=False,wait_for=executor.wait)
        options.update(changes)
        result = scheduler.search(model,cover or dict(leaves={},unresolved={'':{}}),executor,**options)
        return result, snapshots, executor

    def test_refills_without_waiting_for_slow_sibling(self):
        def choose(requests):
            return min(requests,key=lambda r:(r['job'][0]=='0',r['id']))
        executor=ControlledExecutor(self.model(),choose,{'':False,'1':False})
        state,snapshots,_=self.run_search(executor)
        self.assertEqual(executor.completions,['','1','10','11','0'])
        self.assertEqual(set(state['leaves']),{'0','10','11'})
        self.assertEqual(state['unresolved'],{})
        self.assertEqual(state['visited'],5)
        self.assertEqual(executor.peak,2)
        for path, completed_before in executor.submissions:
            if path in ('10','11'):
                self.assertNotIn('0',completed_before)
        self.assertTrue(any('0' in s['completion_scheduler']['in_flight']
                            and '10' in s['leaves'] for s in snapshots))
        for state in snapshots:
            scheduler.core.proof.sc.partition(self.model(),state['leaves'],state['unresolved'])
            self.assertTrue(set(state['completion_scheduler']['in_flight']) <= set(state['unresolved']))

    def test_budget_counts_submissions_and_keeps_unscheduled_children(self):
        executor=ControlledExecutor(self.model(),outcomes={'':False,'0':False,'1':False})
        state,_,_=self.run_search(executor,max_cells=2)
        self.assertEqual(len(executor.submissions),2)
        self.assertEqual(state['completion_scheduler']['submitted'],2)
        self.assertEqual(state['visited'],2)
        self.assertEqual(set(state['unresolved']),{'00','01','1'})
        self.assertEqual(state['completion_scheduler']['in_flight'],[])

    def test_prior_leaves_are_rechecked_and_prior_work_preserved(self):
        cover=dict(leaves={'0':dict(witness={'old':1},upper=[1,-99999])},
                   unresolved={'1':{}},visited=37)
        saved=copy.deepcopy(cover)
        state,snapshots,_=self.run_search(cover=cover)
        self.assertEqual(cover,saved)
        self.assertEqual(snapshots[0]['leaves'],{})
        self.assertEqual(set(snapshots[0]['unresolved']),{'0','1'})
        self.assertEqual(state['leaves']['0']['witness'],{'fresh':{'old':1}})
        self.assertEqual(state['visited'],39)

    def hint_cover(self):
        return dict(leaves={'0':dict(witness={'old':1},proposal=-60.,upper=[1,-99999])},
                    unresolved={'1':{}},visited=150,
                    component_sha256=scheduler.core.component_digest(self.model().components))

    def test_opt_in_retains_hints_without_submitting_or_using_saved_bounds(self):
        cover=self.hint_cover();original=copy.deepcopy(cover)
        state,snapshots,executor=self.run_search(cover=cover,retain_search_leaves=True)
        self.assertEqual(cover,original)
        self.assertEqual([path for path,_ in executor.submissions],['1'])
        self.assertEqual(state['visited'],151)
        self.assertTrue(state['leaves']['0']['retained_search_hint'])
        self.assertNotIn('upper',state['leaves']['0'])
        metadata=state['completion_scheduler']
        self.assertEqual(metadata['retained_search_leaf_count'],1)
        self.assertEqual(metadata['retained_search_leaf_paths'],['0'])
        self.assertEqual(metadata['freshly_checked_leaf_count'],1)
        expected=hashlib.sha256(json.dumps(cover,sort_keys=True,
            separators=(',', ':'),allow_nan=False).encode()).hexdigest()
        self.assertEqual(metadata['source_cover_sha256'],expected)
        self.assertEqual(metadata['source_visited'],150)
        self.assertIn('NOT freshly verified',metadata['status'])
        self.assertIn('EVERY leaf',metadata['status'])
        for snapshot in snapshots:
            scheduler.core.proof.sc.partition(self.model(),snapshot['leaves'],snapshot['unresolved'])

    def test_hint_retention_is_off_by_default_even_for_completed_source(self):
        cover=self.hint_cover()
        state,_,executor=self.run_search(cover=cover)
        self.assertEqual([path for path,_ in executor.submissions],['0','1'])
        self.assertFalse(state['completion_scheduler']['retain_search_leaves'])
        self.assertEqual(state['completion_scheduler']['retained_search_leaf_count'],0)
        self.assertNotIn('retained_search_hint',state['leaves']['0'])

    def test_source_leaf_missing_strict_target_is_enqueued_for_fresh_check(self):
        for proposal in (-32.,-31.,0.,10.):
            cover=self.hint_cover();cover['leaves']['0']['proposal']=proposal
            state,_,executor=self.run_search(cover=cover,retain_search_leaves=True)
            self.assertEqual([path for path,_ in executor.submissions],['0','1'])
            self.assertEqual(state['completion_scheduler']['retained_search_leaf_count'],0)
            self.assertNotIn('retained_search_hint',state['leaves']['0'])

    def test_all_retained_source_is_still_only_a_search_checkpoint(self):
        cover=self.hint_cover()
        cover['leaves']={'':cover['leaves']['0']};cover['unresolved']={}
        state,_,executor=self.run_search(cover=cover,retain_search_leaves=True)
        self.assertEqual(executor.submissions,[])
        self.assertEqual(state['visited'],150)
        self.assertEqual(state['unresolved'],{})
        self.assertEqual(state['completion_scheduler']['freshly_checked_leaf_count'],0)
        self.assertIn('NOT freshly verified',state['completion_scheduler']['status'])

    def test_bad_retained_witness_proposal_comparison_or_scope_is_rejected(self):
        for change in (dict(witness={}),dict(witness=None),dict(proposal=None),
                       dict(proposal=True),dict(proposal=float('nan')),dict(proposal=float('inf'))):
            cover=self.hint_cover();cover['leaves']['0'].update(change)
            executor=ControlledExecutor(self.model())
            with self.assertRaises(ValueError):
                self.run_search(executor,cover,retain_search_leaves=True)
            self.assertEqual(executor.submissions,[])
        for change in (dict(component_sha256='wrong'),dict(component_sha256=None),
                       dict(completion_scheduler=None),dict(completion_scheduler=[]),
                       dict(completion_scheduler=dict(proof_scope={'wrong':True}))):
            cover=self.hint_cover();cover.update(change)
            with self.assertRaises(ValueError):self.run_search(cover=cover,retain_search_leaves=True)

    def test_hint_retention_requires_boolean_and_preserves_worker_budget(self):
        with self.assertRaises(ValueError):
            self.run_search(cover=self.hint_cover(),retain_search_leaves=1)
        executor=ControlledExecutor(self.model(),outcomes={'1':False})
        state,_,_=self.run_search(executor,self.hint_cover(),retain_search_leaves=True,max_cells=1)
        self.assertEqual(len(executor.submissions),1)
        self.assertEqual(set(state['unresolved']),{'10','11'})
        self.assertEqual(set(state['leaves']),{'0'})
        self.assertEqual(state['completion_scheduler']['submitted'],1)

    def test_max_depth_stops_splitting(self):
        executor=ControlledExecutor(self.model(),outcomes={'':False,'0':False,'1':False})
        state,_,_=self.run_search(executor,max_depth=1)
        self.assertEqual(set(state['unresolved']),{'0','1'})
        self.assertEqual(len(executor.submissions),3)

    def test_worker_error_keeps_failed_and_inflight_cells(self):
        cover=dict(leaves={},unresolved={'0':{},'1':{}})
        executor=ControlledExecutor(self.model(),outcomes={'0':RuntimeError('worker failed')})
        snapshots=[]
        with self.assertRaisesRegex(RuntimeError,'worker failed'):
            self.run_search(executor,cover,checkpoint=snapshots.append)
        state=snapshots[-1]
        self.assertEqual(state['leaves'],{})
        self.assertEqual(set(state['unresolved']),{'0','1'})
        self.assertTrue(next(f for f,r in executor.requests.items() if r['job'][0]=='1').cancelled())

    def test_submit_failure_preserves_unsent_cell(self):
        executor=ControlledExecutor(self.model())
        snapshots=[]
        with patch.object(executor,'submit',side_effect=RuntimeError('submit failed')):
            with self.assertRaisesRegex(RuntimeError,'submit failed'):
                self.run_search(executor,checkpoint=snapshots.append)
        self.assertEqual(set(snapshots[-1]['unresolved']),{''})
        self.assertEqual(snapshots[-1]['completion_scheduler']['submitted'],0)

    def test_corrupted_scope_point_cell_or_endpoint_is_rejected(self):
        changes=[dict(component_sha256='wrong'),dict(request_id=99),dict(request_id=True),
                 dict(path='1'),dict(cell=['1/2','1/2']),dict(witness=None),dict(proposal=float('nan')),
                 dict(proof_scope={}),dict(upper=[1,-32]),dict(upper=[0,0]),dict(upper=[-1,-60]),
                 dict(upper=[True,-60])]
        for change in changes:
            with self.subTest(change=change):
                executor=ControlledExecutor(self.model(),outcomes={'':change})
                snapshots=[]
                with self.assertRaises(ArithmeticError):
                    self.run_search(executor,checkpoint=snapshots.append)
                self.assertEqual(snapshots[-1]['leaves'],{})
                self.assertEqual(set(snapshots[-1]['unresolved']),{''})

    def test_wrong_or_empty_completion_set_is_rejected(self):
        for invalid in (set(),{Future()}):
            with self.assertRaises(ArithmeticError):
                self.run_search(wait_for=lambda *a,**k:(invalid,set()))

    def test_interrupt_preserves_complete_partition_and_cancels_owned_work(self):
        executor=ControlledExecutor(self.model())
        snapshots=[]
        def interrupt(*args,**kwargs):
            raise KeyboardInterrupt()
        with self.assertRaises(KeyboardInterrupt):
            self.run_search(executor,checkpoint=snapshots.append,wait_for=interrupt)
        self.assertEqual(set(snapshots[-1]['unresolved']),{''})
        self.assertTrue(all(f.cancelled() for f in executor.requests))

    def test_snapshot_with_inflight_work_can_be_resumed_as_unresolved(self):
        executor=ControlledExecutor(self.model(),outcomes={'':False})
        _,snapshots,_=self.run_search(executor)
        interrupted=next(s for s in snapshots if len(s['completion_scheduler']['in_flight'])==2)
        state,_,resumed=self.run_search(cover=interrupted)
        self.assertEqual(set(state['leaves']),{'0','1'})
        self.assertEqual(len(resumed.submissions),2)
        self.assertEqual(state['visited'],3)

    def test_parent_restriction_error_preserves_parent_cell(self):
        executor=ControlledExecutor(self.model(),outcomes={'':False})
        snapshots=[]
        with patch.object(scheduler.core,'restrict_witness',side_effect=ValueError('bad restriction')):
            with self.assertRaisesRegex(ValueError,'bad restriction'):
                self.run_search(executor,reuse_parent_witness=True,checkpoint=snapshots.append)
        self.assertEqual(set(snapshots[-1]['unresolved']),{''})
        self.assertEqual(snapshots[-1]['visited'],0)

    def test_duplicate_future_is_rejected_and_second_cell_retained(self):
        executor=ControlledExecutor(self.model())
        shared=Future();snapshots=[]
        cover=dict(leaves={},unresolved={'0':{},'1':{}})
        with patch.object(executor,'submit',return_value=shared):
            with self.assertRaises(ArithmeticError):
                self.run_search(executor,cover,checkpoint=snapshots.append)
        self.assertEqual(set(snapshots[-1]['unresolved']),{'0','1'})

    def test_checkpoint_callback_cannot_mutate_internal_state(self):
        snapshots=[]
        def mutate(state):
            snapshots.append(copy.deepcopy(state))
            state['leaves'].clear()
            state['unresolved'].clear()
            state['completion_scheduler']['proof_scope']['threshold']=-1
        state,_,_=self.run_search(checkpoint=mutate)
        self.assertEqual(set(state['leaves']),{''})
        self.assertEqual(state['completion_scheduler']['proof_scope']['threshold'],146800)

    def test_parent_witness_reuse_is_opt_in(self):
        for reuse in (False,True):
            executor=ControlledExecutor(self.model(),outcomes={'':False})
            state,_,_=self.run_search(executor,reuse_parent_witness=reuse)
            children=[r for r in executor.requests.values() if r['job'][0]]
            self.assertEqual([r['job'][3] for r in children],
                             [{'fresh':None}]*2 if reuse else [None,None])

    def test_scope_rejects_different_ensemble_or_degenerate_root(self):
        for key,value in (('threshold',True),('q_min',0),('root',(Q(1,2),Q(1,2))),
                          ('variance_shuffle',False),('regional_count',False),
                          ('variance_bins',0),('data',dict(bits=19,windows=32,updates=4))):
            model=self.model();setattr(model,key,value)
            with self.assertRaises(ValueError):scheduler.scope(model,384)
        for precision in (True,192,384.):
            with self.assertRaises(ValueError):scheduler.scope(self.model(),precision)

    def test_invalid_limits_fail_before_submission(self):
        for changes in (dict(workers=True),dict(workers=0),dict(max_cells=0),dict(max_depth=0),
                        dict(target_bits=19),dict(reuse_parent_witness=1)):
            executor=ControlledExecutor(self.model())
            with self.assertRaises(ValueError):self.run_search(executor,**changes)
            self.assertEqual(executor.submissions,[])

    def test_wrapper_checks_worker_model_before_core_evaluation(self):
        model=self.model();request=dict(id=1,scope=scheduler.scope(model,384),
            job=('',model.root,32,None))
        with patch.object(scheduler.core,'_model',model,create=True),\
             patch.object(scheduler.core,'_precision',384,create=True),\
             patch.object(scheduler.core,'evaluate_worker',return_value={'value':'fresh'}) as work:
            result=scheduler.evaluate_worker(request)
            self.assertEqual(result['value'],'fresh')
            self.assertEqual(result['proof_scope'],request['scope'])
            work.assert_called_once_with(request['job'])
            request['scope']['threshold']+=1
            with self.assertRaises(ArithmeticError):scheduler.evaluate_worker(request)
            self.assertEqual(work.call_count,1)

    def test_initializer_delegates_through_the_same_named_module(self):
        arguments=(['components'],{'record':True},0,384,'logs')
        with patch.object(scheduler.core,'initialize') as initialize:
            scheduler.initialize(*arguments)
            initialize.assert_called_once_with(*arguments)
        self.assertEqual(scheduler.core.__name__,'shared_relaxed_parallel')

    def test_point_or_wrong_path_request_cannot_enter_the_worker(self):
        model=self.model()
        for path,cell in (('',(Q(1,2),Q(1,2))),('0',model.root),('x',model.root)):
            request=dict(id=1,scope=scheduler.scope(model,384),job=(path,cell,32,None))
            with patch.object(scheduler.core,'_model',model,create=True),\
                 patch.object(scheduler.core,'_precision',384,create=True),\
                 patch.object(scheduler.core,'evaluate_worker') as work:
                with self.assertRaises(ArithmeticError):scheduler.evaluate_worker(request)
                work.assert_not_called()

    def test_actual_worker_point_success_cannot_accept_interval(self):
        model=self.model();model.point_probes=[];model.split_on_midpoint=True
        model.outward=lambda cell,witness:arb(2)**(-80 if cell[0]==cell[1] else 5)
        model.proposal=lambda cell:(_ for _ in ()).throw(AssertionError('unexpected new proposal'))
        job=('',model.root,32,dict(parameters=['1/100','0','0']))
        request=dict(id=1,scope=scheduler.scope(model,384),job=job)
        previous=ctx.prec
        try:
            result=scheduler.core.evaluate(model,job,384,request['scope']['component_sha256'])
            result.update(request_id=1,proof_scope=request['scope'])
            self.assertIsNone(scheduler.checked_result(result,request)['upper'])
        finally:
            ctx.prec=previous


if __name__ == '__main__':
    unittest.main()
