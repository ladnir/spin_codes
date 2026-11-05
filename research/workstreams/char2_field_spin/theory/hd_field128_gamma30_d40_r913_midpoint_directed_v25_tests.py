"""Pure/mocked orchestration tests; no regional tables or proof arithmetic."""
from contextlib import nullcontext
from copy import deepcopy
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

import hd_field128_gamma30_d40_r913_midpoint_directed_v25 as mode
import hd_field128_gamma30_d40_r913_directed_v25 as validation


def write_json(path, value):
    with Path(path).open('x', encoding='utf-8', newline='\n') as output:
        output.write(json.dumps(value, indent=2, allow_nan=False)+'\n')


def complete_result():
    return dict(status='complete', complete_selected_weight=True,
        covers_all_supports_at_weight=True, covers_all_even_odd_splits=True,
        directed=True, whole_code_certificate=False, active_groups=913,
        completed_even_rows=457, profile_pairs_evaluated=208849,
        pending_even_full_range=None, requested_failure_bits=54,
        meets_requested_selected_weight_target=True, log2_upper_rational='-59')


class PureTests(unittest.TestCase):
    def test_admission_requires_complete_original_weight_and_actual_ceiling(self):
        self.assertTrue(mode.admitted(complete_result()))
        for key,value in [('status','incomplete'), ('complete_selected_weight',False),
                ('covers_all_supports_at_weight',False), ('covers_all_even_odd_splits',False),
                ('directed',False), ('whole_code_certificate',True), ('active_groups',912),
                ('completed_even_rows',456), ('profile_pairs_evaluated',208848),
                ('pending_even_full_range',[456,456]), ('requested_failure_bits',40),
                ('log2_upper_rational','-53'), ('log2_upper_rational',None),
                ('meets_requested_selected_weight_target',False)]:
            result=complete_result();result[key]=value
            self.assertFalse(mode.admitted(result),key)

    def test_exact_output_budget_exclusive_containment_and_charged_failure(self):
        with TemporaryDirectory() as directory:
            root=Path(directory);guard=Mock();value={'tiny':True}
            size=len((json.dumps(value,indent=2)+'\n').encode())
            with patch.object(mode,'MAX_OUTPUT_BYTES',size):
                budget=mode.OutputBudget(root,guard,write_json)
                budget.write_new(root/'one.json',value)
                self.assertEqual(budget.used_bytes,size);guard.assert_called_once_with(size)
                with self.assertRaises(FileExistsError):budget.write_new(root/'one.json',value)
                with self.assertRaises(ValueError):budget.write_new(root/'outside'/'two.json',value)
                with self.assertRaisesRegex(RuntimeError,'2MiB'):budget.write_new(root/'two.json',value)
                self.assertFalse((root/'two.json').exists())
            broken=mode.OutputBudget(root,guard,Mock(side_effect=OSError('disk')))
            with self.assertRaises(OSError):broken.write_new(root/'broken.json',value)
            self.assertEqual(broken.used_bytes,size)
            denied=mode.OutputBudget(root,Mock(side_effect=RuntimeError('reserve')),write_json)
            with self.assertRaisesRegex(RuntimeError,'reserve'):denied.write_new(root/'denied.json',value)
            self.assertFalse((root/'denied.json').exists())

    def test_helper_identity_and_source_conflicts(self):
        with TemporaryDirectory() as directory, patch.object(mode,'HERE',Path(directory)):
            root=Path(directory);raw=b'pass\n';(root/'helper.py').write_bytes(raw)
            with patch.object(mode.importlib,'import_module',return_value='imported') as importer:
                self.assertEqual(mode.pinned_module('helper',sha256(raw).hexdigest()),'imported')
                with self.assertRaises(RuntimeError):mode.pinned_module('helper','0'*64)
                self.assertEqual(importer.call_count,1)
            (root/mode.NAME).write_bytes(raw)
            validator=SimpleNamespace(source_pins=lambda proposal:{'x.py':'1'*64})
            producer=SimpleNamespace(source_pins=lambda:{'x.py':'2'*64})
            with self.assertRaisesRegex(RuntimeError,'conflict'):
                mode.source_pins({},validator,producer)

    def test_inherited_proposal_validation_is_byte_pinned_and_strict(self):
        self.assertEqual(sha256(Path(validation.__file__).read_bytes()).hexdigest(),mode.VALIDATOR_SHA256)
        with TemporaryDirectory() as directory, patch.object(validation,'HERE',Path(directory)):
            root=Path(directory)
            proposal=dict(active_groups=913,state_dimension=30,target_distance='2/5',
                field_order=1<<128,construction=validation.CONSTRUCTION,marker_count=83,
                witness=dict(z='1/4',alpha='1/8',grid=128,state_dimension=30,
                    markers=[['1',str(k+1)] for k in range(83)]))
            raw=json.dumps(proposal).encode();(root/'proposal.json').write_bytes(raw)
            got,parameters,source=validation.load_proposal('proposal.json',sha256(raw).hexdigest())
            self.assertEqual(got,proposal);self.assertEqual(parameters[0],F(1,4))
            proposal['active_groups']=914
            raw=json.dumps(proposal).encode();(root/'proposal.json').write_bytes(raw)
            with self.assertRaises(ValueError):validation.load_proposal('proposal.json',sha256(raw).hexdigest())


class MockedRunTests(unittest.TestCase):
    def exercise(self, action='complete'):
        temporary=TemporaryDirectory();self.addCleanup(temporary.cleanup)
        root=Path(temporary.name);events=[];captured={};bank=SimpleNamespace(check_pins=Mock(),metadata=lambda:{'test':True})
        def limit():events.append('limit');return nullcontext(4<<30)
        validator=SimpleNamespace(relative_name=lambda value,**kw:value,
            load_proposal=Mock(return_value=({},(F(1,4),F(1,8),[(F(1),F(1))]*83),{'path':'p.json','sha256':'a'*64})),
            reserve_guard=Mock(),linux_limits=limit,check_pins=Mock(),MAX_SECONDS=600,
            MAX_TABLE_BYTES=64<<20,RESERVE=20<<30,Q=1<<128,CONSTRUCTION='test',
            helper=lambda:SimpleNamespace(portable=SimpleNamespace(runtime_metadata=lambda:{'mocked':True})))
        def receipt(actual,weight,**options):
            self.assertIs(actual,bank);self.assertEqual(weight,913);captured.update(options)
            options['context_checkpoint']({'context':{'test':True},'context_sha256':'c'*64})
            for first in range(2 if action=='partial' else 457):
                options['checkpoint']({'even_full':first,'context_sha256':'c'*64,'complete_selected_weight':False})
            if action=='exception':raise MemoryError('mock resource stop')
            if action=='mutation':(root/'out'/'row_0000.json').write_bytes(b'changed')
            result=complete_result()
            if action=='partial':
                result.update(status='incomplete',complete_selected_weight=False,
                    covers_all_supports_at_weight=False,covers_all_even_odd_splits=False,
                    completed_even_rows=2,profile_pairs_evaluated=914,pending_even_full_range=[2,456],
                    log2_upper_rational=None,meets_requested_selected_weight_target=False)
            return result
        producer=SimpleNamespace(provider=SimpleNamespace(AliasBank=Mock(return_value=bank)),
            weight_receipt=receipt,write_new=write_json)
        def helper(name,digest):events.append('heavy');return producer
        with patch.object(mode,'HERE',root),patch.object(mode,'validation',return_value=validator),\
                patch.object(mode,'pinned_module',side_effect=helper),\
                patch.object(mode,'source_pins',return_value={mode.NAME:'d'*64}):
            if action in ('exception','mutation'):
                with self.assertRaises((MemoryError,RuntimeError)):
                    mode.run('p.json','a'*64,'out')
                self.assertGreaterEqual(bank.check_pins.call_count,2)
                if action=='exception':self.assertFalse((root/'out'/'summary.json').exists())
                return
            result=mode.run('p.json','a'*64,'out')
        self.assertEqual(events,['limit','heavy'])
        self.assertEqual(captured['failure_bits'],54);self.assertEqual(captured['max_singletons'],64)
        self.assertTrue(0<captured['max_seconds']<=595)
        self.assertEqual(producer.provider.AliasBank.call_args.kwargs['max_table_bytes'],64<<20)
        self.assertEqual(producer.provider.AliasBank.call_args.kwargs['state'],30)
        validator.reserve_guard.assert_any_call(2<<20)
        self.assertLessEqual(validator.check_pins.call_count,8) # Not per row for wrapper closure.
        self.assertEqual(len(list((root/'out').glob('row_*.json'))),2 if action=='partial' else 457)
        self.assertFalse(json.loads((root/'out'/'driver.json').read_text())['fresh_bank_rebuilt'])
        self.assertTrue(json.loads((root/'out'/'bank.json').read_text())['fresh_bank_rebuilt'])
        self.assertFalse(result['whole_code_certificate'])
        self.assertEqual(mode.admitted(result),action=='complete')

    def test_complete_original913_full_profile_call(self):self.exercise()
    def test_partial_rows_never_admitted(self):self.exercise('partial')
    def test_resource_failure_preserves_rows_and_checks_exit(self):self.exercise('exception')
    def test_output_mutation_rejected_at_exit(self):self.exercise('mutation')


if __name__=='__main__':
    unittest.main()
