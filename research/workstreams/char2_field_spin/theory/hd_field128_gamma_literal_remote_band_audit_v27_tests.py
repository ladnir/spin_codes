"""Pure/mocked literal-boundary tests; no coefficients or production audit."""
from copy import deepcopy
from hashlib import sha256
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import hd_field128_gamma_literal_remote_band_audit_v27 as mode


def seed():
    return dict(schema='gamma-portable-exact-witness-v1', active_groups=960,
        construction=mode.typed.CONSTRUCTION, field_order=1<<128, state_dimension=30,
        target_distance='2/5', directed=False, selected_weight_certificate=False,
        whole_code_certificate=False,
        witness=dict(z='1/8',alpha='1/128',markers=[['2','3']],grid=128,state_dimension=30),
        source=dict(path='C:\\missing-ancestry\\never-open.json',sha256='f'*64),
        source_sha256={'deliberately_unopened_historical_helper.py':'e'*64})


def put(path, value):
    raw=(json.dumps(value,indent=2)+'\n').encode()
    path.write_bytes(raw)
    return dict(path=str(path),sha256=sha256(raw).hexdigest())


class Tests(unittest.TestCase):
    def test_exact_literal_boundary_does_not_open_history(self):
        with tempfile.TemporaryDirectory() as directory, patch.object(mode,'HERE',Path(directory)):
            path=Path(directory)/'seed.json'; pin=put(path,seed())
            context=dict(witness=seed()['witness'],witness_source=dict(path='seed.json',sha256=pin['sha256']))
            inputs=mode.typed.Inputs()
            got,unused=mode._seed(context,inputs,pin['sha256'])
            self.assertEqual(got,seed())
            self.assertEqual(len(inputs.pins),1)
            inputs.check()
            with self.assertRaisesRegex(ValueError,'explicit expected digest'):
                mode._seed(context,inputs,'0'*64)
            bad=deepcopy(context);bad['witness']['z']='1/16'
            with self.assertRaisesRegex(ValueError,'literal mismatch'):
                mode._seed(bad,mode.typed.Inputs(),pin['sha256'])
            invalid=seed();invalid['witness']['markers']=[['2/1','3']]
            badpin=put(path,invalid);context['witness_source']['sha256']=badpin['sha256']
            with self.assertRaisesRegex(ValueError,'canonical'):
                mode._seed(context,mode.typed.Inputs(),badpin['sha256'])

    def test_checkpoint_reconstruction_and_partial_flags(self):
        record=dict(checkpoints=2,nodes=[dict(id=0,status='accepted')],context_sha256='a'*64)
        flags=dict(complete_requested_band=False,covers_all_supports_in_band=False,directed=False,
            batch_source_verified_at_exit=False,whole_code_certificate=False,log2_upper_rational=None)
        rows=[dict(flags,sequence=0,context_sha256='a'*64,changed_nodes=[]),
              dict(flags,sequence=1,context_sha256='a'*64,changed_nodes=record['nodes'])]
        class Inputs:
            def load(self,path):
                return rows[int(path.stem.split('_')[1])],{}
        budget=SimpleNamespace(check=lambda:None)
        self.assertEqual(mode._checkpoints(record,Path('.'),Inputs(),budget),2)
        rows[1]['directed']=True
        with self.assertRaisesRegex(ValueError,'overclaims'):
            mode._checkpoints(record,Path('.'),Inputs(),budget)
        rows[1]['directed']=False;rows[1]['changed_nodes']=[dict(id=0,status='split')]
        with self.assertRaisesRegex(ValueError,'final tree'):
            mode._checkpoints(record,Path('.'),Inputs(),budget)

    def _fixture(self,directory):
        root=Path(directory); portable=seed(); seed_pin=put(root/'seed.json',portable)
        witness=dict(portable['witness'],precision=256,packet_mantissa_bits=80)
        source=dict(path='seed.json',sha256=seed_pin['sha256'])
        context=dict(witness=witness,witness_source=source,construction=mode.typed.CONSTRUCTION,
            field_order=1<<128,state_dimension=30,target_distance='2/5',output_length=1<<21,
            source_dimension=1<<20,outer_measure='invertible-toeplitz',bank_source_sha256={},
            all_pair_before_alpha=True,old_numerical_ceilings_reused=False,construction_unchanged=True)
        pins={mode.runner.NAME:'a'*64,mode.runner.previous.NAME:'b'*64,
              mode.runner.controller.NAME:'c'*64}
        driver=mode.runner.driver_record(pins,source,968,999,4096,960,64)
        bank=dict(context,runner_source_sha256=pins,fresh_bank_rebuilt=True,directed=False,
            whole_code_certificate=False,memory_limit=dict(mechanism='Linux RLIMIT_AS',
                requested_bytes=4<<30,applied_bytes=4<<30),runtime=dict(OPENBLAS_NUM_THREADS='1'))
        record=dict(context=context,context_sha256='d'*64,requested_weight_range=[968,999],
            max_nodes=4096,cell_failure_bits=64,runner=driver['driver'],runner_source_sha256=pins,
            controller=dict(path=mode.runner.controller.NAME,sha256=pins[mode.runner.controller.NAME]),
            controller_source_sha256={},runner_admitted=True,construction_changed=False,
            total_compute_seconds=50.,checkpoints=2,nodes=[dict(id=0,status='accepted')],
            status='complete',complete_requested_band=True,directed=True,whole_code_certificate=False,
            requested_failure_bits=54,meets_requested_band_target=True,pending_frontier=[],
            zero=False,log2_upper_rational='-60')
        summary_pin=put(root/'summary.json',record)
        put(root/'driver.json',driver);put(root/'bank.json',bank)
        put(root/'context.json',dict(context=context,context_sha256='d'*64,directed=False,
            complete_requested_band=False,whole_code_certificate=False))
        for sequence in range(2):
            put(root/f'node_{sequence:06d}.json',dict(sequence=sequence,context_sha256='d'*64,
                changed_nodes=[] if sequence==0 else record['nodes'],complete_requested_band=False,
                covers_all_supports_in_band=False,directed=False,batch_source_verified_at_exit=False,
                whole_code_certificate=False,log2_upper_rational=None))
        numeric=dict(weights=list(range(968,1000)),normalized_histograms=1,cells_checked=1,
            accepted_cells=1,positive_count_zero_cells=0,union='-60',union_enclosure='-60',zero=False)
        return record,pins,summary_pin,seed_pin,numeric

    def _audit(self,directory,*,mutation=None):
        record,pins,summary,seed_pin,numeric=self._fixture(directory)
        root=Path(directory)
        if mutation:
            mutation(root)
        with patch.object(mode,'HERE',root),patch.object(mode,'source_pins',return_value={}), \
                patch.object(mode.typed,'check_source_pins'),patch.object(mode,'_sources'), \
                patch.object(mode.typed.previous,'_metadata') as metadata, \
                patch.object(mode.runner,'source_pins',return_value=pins), \
                patch.object(mode.runner.controller,'source_pins',return_value={}), \
                patch.object(mode.typed.previous,'audit_numeric',return_value=numeric) as math:
            report=mode.audit_summary(root/'summary.json',expected_summary_sha256=summary['sha256'],
                expected_seed_sha256=seed_pin['sha256'])
            metadata.assert_called_once_with(record)
            self.assertEqual(math.call_args.args[:2],(record,4096))
            return report

    def test_mock_lifecycle_explicit_scope_and_mandatory_bank(self):
        with tempfile.TemporaryDirectory() as directory:
            report=self._audit(directory)
            self.assertEqual(report['schema'],mode.SCHEMA)
            self.assertFalse(report['historical_origin_verified'])
            self.assertFalse(report['historical_bounds_used'])
            self.assertFalse(report['calibration_history_replayed'])
            self.assertFalse(report['coefficient_replay_in_this_audit'])
            self.assertEqual(report['directed_union_precision'],384)
            self.assertEqual(report['requested_weights'],list(range(968,1000)))
            self.assertEqual(len(report['checked_file_pins']),7)
            self.assertNotIn('checked_non_json_file_pins',report)
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaises(FileNotFoundError):
                self._audit(directory,mutation=lambda root:(root/'bank.json').unlink())

    def test_bank_and_driver_tampering_fail_closed(self):
        def change(root,name,key,value):
            path=root/name;row=json.loads(path.read_bytes());row[key]=value;put(path,row)
        for name,key,value in [('bank.json','fresh_bank_rebuilt',False),
                ('bank.json','witness',{}),('driver.json','cell_failure_bits',72)]:
            with tempfile.TemporaryDirectory() as directory:
                with self.assertRaises(ValueError):
                    self._audit(directory,mutation=lambda root:change(root,name,key,value))

    def test_source_map_guard_and_output_budget(self):
        with patch.object(mode.typed,'_source_map') as scan:
            mode._sources(dict(source_sha256={'x':'y'},nested=[dict(bank_source_sha256={'a':'b'})]),{})
            self.assertEqual(scan.call_count,2)
        with patch.object(mode.shutil,'disk_usage',return_value=SimpleNamespace(free=mode.RESERVE)), \
                patch.object(mode.typed.producer.profiles,'write_new') as write:
            with self.assertRaisesRegex(ValueError,'reserve'):
                mode.write_remote(Path('audit.json'),{})
            write.assert_not_called()
        with patch.object(mode,'MAX_OUTPUT_BYTES',8):
            with self.assertRaisesRegex(ValueError,'one MiB'):
                mode.write_remote(Path('audit.json'),{'too':'large'})

    def test_frozen_source_identities(self):
        self.assertEqual(sha256(Path(mode.typed.__file__).read_bytes()).hexdigest(),mode.TYPED_SHA256)
        self.assertEqual(sha256(Path(mode.runner.__file__).read_bytes()).hexdigest(),mode.RUNNER_SHA256)


if __name__=='__main__':
    unittest.main()
