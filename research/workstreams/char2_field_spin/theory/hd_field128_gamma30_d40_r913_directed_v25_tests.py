"""Pure/mocked driver tests; no coefficient bank, evaluator, or proof run."""
from contextlib import ExitStack, nullcontext
from copy import deepcopy
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

import hd_field128_gamma30_d40_r913_directed_v25 as mode
import hd_field128_gamma_portable_band_runner_v25 as frozen_controller


def proposal():
    return dict(active_groups=913, state_dimension=30, construction=mode.CONSTRUCTION,
        target_distance='2/5', marker_count=83, source_sha256={},
        witness=dict(z='1/4', alpha='1/8', grid=128, state_dimension=30,
                     markers=[['1', str(k+1)] for k in range(83)]),
        source=dict(path='C:/unchanged/windows/ancestor.json', sha256='a'*64))


class ValidationTests(unittest.TestCase):
    def test_exact_proposal_and_fixed_scope(self):
        with TemporaryDirectory() as directory, patch.object(mode,'HERE',Path(directory)):
            path=Path(directory)/'proposal.json'
            for key,value in ((None,None),('active_groups',912),('state_dimension',31),
                              ('field_order',256),('target_distance','9/20'),('marker_count',82)):
                record=proposal()
                if key: record[key]=value
                raw=json.dumps(record).encode(); path.write_bytes(raw)
                if key:
                    with self.assertRaises(ValueError): mode.load_proposal('proposal.json',sha256(raw).hexdigest())
                else:
                    got, parameters, source=mode.load_proposal('proposal.json',sha256(raw).hexdigest())
                    self.assertEqual(got,record); self.assertEqual(parameters[:2],(F(1,4),F(1,8)))
                    self.assertEqual(len(parameters[2]),83); self.assertEqual(source['path'],'proposal.json')
            record=proposal(); record['witness']['markers'][0][0]='2/2'
            raw=json.dumps(record).encode();path.write_bytes(raw)
            with self.assertRaises(ValueError): mode.load_proposal('proposal.json',sha256(raw).hexdigest())
            with self.assertRaises(ValueError): mode.load_proposal('proposal.json','0'*64)

    def test_relative_paths_cell_thresholds_and_full_output_preflight(self):
        for name in ('../a.json','/tmp/a.json','C:/a.json','a\\b.json','sub//a.json','a./b.json'):
            with self.assertRaises(ValueError): mode.relative_name(name,json_file=True)
        self.assertEqual(mode.relative_name('old/proposal.json',json_file=True),'old/proposal.json')
        for value in (False,53,57,72):
            with self.assertRaises(ValueError): mode.run('missing.json','0'*64,'out',cell_failure_bits=value)
        for nodes in (0,4097,True):
            with self.assertRaises(ValueError): mode.run('missing.json','0'*64,'out',max_nodes=nodes)
        with patch.object(mode.shutil,'disk_usage',return_value=SimpleNamespace(free=mode.RESERVE+mode.MAX_OUTPUT_BYTES-1)):
            with self.assertRaises(RuntimeError): mode.reserve_guard(mode.MAX_OUTPUT_BYTES)
        with patch.object(mode.shutil,'disk_usage',return_value=SimpleNamespace(free=mode.RESERVE+mode.MAX_OUTPUT_BYTES)):
            mode.reserve_guard(mode.MAX_OUTPUT_BYTES)

    def test_linux_limits_before_heavy_helper_imports(self):
        resource=SimpleNamespace(RLIMIT_AS=7,RLIM_INFINITY=-1,getrlimit=Mock(return_value=(-1,-1)),setrlimit=Mock())
        with patch.object(mode.sys,'platform','linux'),patch.dict('sys.modules',{'resource':resource}),\
                patch.object(mode.signal,'SIGALRM',14,create=True),patch.object(mode.signal,'ITIMER_REAL',0,create=True),\
                patch.object(mode.signal,'getsignal',return_value='old'),patch.object(mode.signal,'signal') as sig,\
                patch.object(mode.signal,'setitimer',create=True) as timer:
            with mode.linux_limits() as limit:
                self.assertEqual(limit,4<<30); resource.setrlimit.assert_called_once_with(7,(4<<30,4<<30))
                timer.assert_called_once_with(0,600)
                with patch.object(mode.os,'_exit') as exit:
                    sig.call_args[0][1](14,None); exit.assert_called_once_with(124)
            self.assertEqual(timer.call_args[0],(0,0))

    def test_source_and_witness_mutation_checks(self):
        with TemporaryDirectory() as directory, patch.object(mode,'HERE',Path(directory)):
            root=Path(directory); (root/'helper.py').write_bytes(b'pass\n'); (root/'p.json').write_bytes(b'{}')
            pins={'helper.py':sha256(b'pass\n').hexdigest()}; witness=dict(path='p.json',sha256=sha256(b'{}').hexdigest())
            mode.check_pins(pins,witness)
            (root/'p.json').write_bytes(b'changed')
            with self.assertRaises(RuntimeError): mode.check_pins(pins,witness)
            (root/'p.json').write_bytes(b'{}'); (root/'helper.py').write_bytes(b'changed')
            with self.assertRaises(RuntimeError): mode.check_pins(pins,witness)


class OrchestrationTests(unittest.TestCase):
    def exercise(self, *, cell_bits=58, fail=None, complete=True, meets=True):
        temporary=TemporaryDirectory();self.addCleanup(temporary.cleanup)
        root=Path(temporary.name); raw=json.dumps(proposal()).encode(); path=root/'proposal.json';path.write_bytes(raw)
        bank=SimpleNamespace(check_pins=Mock(),metadata=lambda:dict(witness=proposal()['witness']))
        ctor=Mock(return_value=bank); evaluator=Mock(return_value=SimpleNamespace(test_only=True))
        covers=[]
        def cover(value,**options):
            covers.append((value,options))
            options['context_checkpoint'](dict(context={'test':True},directed=False))
            options['checkpoint'](dict(sequence=0,directed=False,complete_requested_band=False))
            if fail=='cover': raise MemoryError('mock resource stop')
            if fail=='source':path.write_bytes(b'changed')
            return dict(status='complete' if complete else 'partial',complete_requested_band=complete,
                directed=complete,meets_requested_band_target=meets and complete,whole_code_certificate=False,
                requested_weight_range=[913,913],evaluated_nodes=1,log2_upper_rational='-55' if complete else None,
                stop_reason=None if complete else 'max_nodes')
        controller=SimpleNamespace(OutputBudget=frozen_controller.OutputBudget,
            cert=SimpleNamespace(AliasBank=ctor),bands=SimpleNamespace(Evaluator=evaluator,__file__=frozen_controller.bands.__file__),
            portable=SimpleNamespace(runtime_metadata=lambda:{'test':True}),cover=cover)
        source_checks=[]
        def check(pins,witness):
            source_checks.append(witness)
            if path.read_bytes()!=raw: raise RuntimeError('parent changed')
        with ExitStack() as stack:
            for item in (patch.object(mode,'HERE',root),patch.object(mode,'linux_limits',return_value=nullcontext(4<<30)),
                patch.object(mode,'helper',return_value=controller),patch.object(mode,'source_pins',return_value={mode.NAME:'f'*64}),
                patch.object(mode,'check_pins',side_effect=check),
                patch.object(mode.shutil,'disk_usage',return_value=SimpleNamespace(free=mode.RESERVE+mode.MAX_OUTPUT_BYTES+4096))):
                stack.enter_context(item)
            if fail:
                with self.assertRaises((RuntimeError,MemoryError)):
                    mode.run('proposal.json',sha256(raw).hexdigest(),'out',cell_failure_bits=cell_bits)
                result=None
            else:
                result=mode.run('proposal.json',sha256(raw).hexdigest(),'out',cell_failure_bits=cell_bits)
                # Entry, driver, bank, context, summary, exit: no checkpoint scan.
                self.assertEqual(len(source_checks),6)
        return root,path,raw,ctor,evaluator,covers,result,bank

    def test_actual_parameters_fixed913_and_aggregate54(self):
        root,path,raw,ctor,evaluator,covers,result,bank=self.exercise()
        self.assertEqual(path.read_bytes(),raw)
        args=ctor.call_args
        self.assertEqual(args.args[:2],(F(1,4),F(1,8)));self.assertEqual(len(args.args[2]),83)
        for key,value in dict(state=30,q=1<<128,distance=F(2,5),grid=128,precision=256,bits=80,
                              outer_measure='invertible-toeplitz',max_table_bytes=64<<20).items():
            self.assertEqual(args.kwargs[key],value)
        self.assertEqual(evaluator.call_args.args[1:],(913,913))
        self.assertEqual(covers[0][1]['aggregate_bits'],54);self.assertEqual(covers[0][1]['cell_failure_bits'],58)
        self.assertEqual(covers[0][1]['max_nodes'],4096)
        self.assertTrue(result['meets_requested_band_target']);self.assertFalse(result['whole_code_certificate'])
        driver=json.loads((root/'out/driver.json').read_bytes())
        self.assertFalse(driver['fresh_bank_rebuilt']);self.assertTrue(driver['fresh_bank_required'])
        self.assertTrue(json.loads((root/'out/bank.json').read_bytes())['fresh_bank_rebuilt'])
        self.assertTrue((root/'out/node_000000.json').exists());self.assertGreaterEqual(bank.check_pins.call_count,3)

    def test_explicit64_and_partial_no_success(self):
        root,_,_,_,_,covers,result,_=self.exercise(cell_bits=64,complete=False)
        self.assertEqual(covers[0][1]['cell_failure_bits'],64)
        self.assertFalse(result['complete_requested_band']);self.assertFalse(result['directed'])
        self.assertIsNone(result['log2_upper_rational'])

    def test_resource_and_source_failures_leave_no_summary(self):
        for failure in ('cover','source'):
            root,*_=self.exercise(fail=failure)
            self.assertTrue((root/'out/node_000000.json').exists())
            self.assertFalse((root/'out/summary.json').exists())

    def test_exit_requires_full_directed_aggregate_not_just_cell_completion(self):
        original=['script','--proposal','p.json','--expected-sha256','a'*64,'--output-dir','out']
        record=dict(status='complete',requested_weight_range=[913,913],evaluated_nodes=1,
            log2_upper_rational='-53',meets_requested_band_target=False,stop_reason=None,
            complete_requested_band=True,directed=True)
        with patch.object(mode.sys,'argv',original),patch.object(mode,'run',return_value=record):
            self.assertEqual(mode.main(),2)
            record['meets_requested_band_target']=True;self.assertEqual(mode.main(),0)
            record['directed']=False;self.assertEqual(mode.main(),2)


if __name__=='__main__':unittest.main()
