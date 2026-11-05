"""Small synthetic bookkeeping tests; no coefficients or transfer replay."""
from copy import deepcopy
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import hd_field128_gamma_selected_index_v26 as mode


def write(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    raw=(json.dumps(value,indent=2)+'\n').encode()
    path.write_bytes(raw)
    return dict(path=str(path.resolve()),sha256=sha256(raw).hexdigest())


def parent_fixture(root,mutate=None):
    pin=write(root/'event.json',dict(synthetic=True))
    entries=[dict(pin,weight_range=[r,r],log2_upper_rational='-100') for r in range(1,247)]
    entries += [dict(pin,weight_range=[247,912],log2_upper_rational='-100'),
        dict(pin,weight_range=[914,7744],log2_upper_rational='-100'),
        dict(pin,weight_range=[8192,8192],log2_upper_rational=None,zero=True)]
    parent=dict(schema='field128-gamma-partial-cover-index-v23',
        construction=mode.midpoint.previous.EXPECTED['construction'],field_order=1<<128,
        source_dimension=1<<20,output_length=1<<21,target_distance='2/5',
        whole_code_certificate=False,directed_aggregation=True,entries=entries,
        receipt_count=249,covered_weight_count=7744,pending_weight_count=448,
        covered_weight_ranges=[[1,912],[914,7744],[8192,8192]],
        pending_weight_ranges=[[913,913],[7745,8191]],
        covered_events_union_log2_upper_rational='-90',aggregation_source_sha256={},
        producer_receipt_counts={})
    if mutate:mutate(parent)
    return write(root/'parent.json',parent)


def point_fixture(root):
    # Two row files only: the frozen full457-row numeric validator is mocked
    # in these adapter-mechanics tests, never in production.
    folder=root/'point'; witness=dict(z='1/4',alpha='1/8',markers=[['1','2']])
    proposal=write(root/'proposal.json',dict(witness=witness))
    driver_path=root/'driver.py';driver_path.write_text('pass\n')
    driver_pin=dict(path='driver.py',sha256=sha256(driver_path.read_bytes()).hexdigest())
    execution={'driver.py':driver_pin['sha256']}
    rows=[dict(even_full=a,odd_full_range=[0,456],profile_pairs_evaluated=457,
        complete_even_profile_row=True,complete_selected_weight=False) for a in range(2)]
    context=dict(active_groups=913,witness=witness)
    context_sha=sha256(json.dumps(context,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    record=dict(context,rows=rows,source_sha256={},bank_source_sha256={},
        context_sha256=context_sha,witness_source=dict(path='proposal.json',sha256=proposal['sha256']),
        execution_driver=driver_pin,execution_source_sha256=execution,
        requested_failure_bits=54,meets_requested_selected_weight_target=True,
        log2_upper_rational='-59')
    summary=write(folder/'summary.json',record)
    context_pin=write(folder/'context.json',dict(context=context,context_sha256=context_sha))
    bank_pin=write(folder/'bank.json',dict(witness=witness,runner_source_sha256=execution,fresh_bank_rebuilt=True))
    driver=write(folder/'driver.json',dict(active_groups=913,marker_count=83,expected_even_rows=457,
        expected_profile_pairs=208849,requested_failure_bits=54,fresh_bank_required=True,
        fresh_bank_rebuilt=False,old_numerical_bounds_reused=False,directed=False,
        whole_code_certificate=False,driver=driver_pin,source_sha256=execution))
    files=[summary,context_pin,bank_pin,driver]
    for a,row in enumerate(rows):
        files.append(write(folder/f'row_{a:04d}.json',dict(row,active_groups=913,
            context_sha256=context_sha,directed=True,whole_code_certificate=False)))
    report=dict(schema='gamma-midpoint-selected-receipts-audit-v1',audit_result='PASS',
        construction=mode.midpoint.previous.EXPECTED['construction'],state_dimension=30,
        target_distance='2/5',selected_weights=[913],precision_bits=384,
        source_sha256=mode.midpoint.source_pins(),coefficient_replay=False,whole_code_certificate=False,
        external_checkpoint_context_reconstructed=True,external_rows_checked=457,
        context_sha256=context_sha,execution_source_sha256=execution,
        checked_file_pins=[summary,proposal],
        audit_adapter_schema='gamma-r913-portable-midpoint-audit-adapter-v1',
        adapter=dict(path=mode.portable_audit.NAME,sha256=mode.PORTABLE_AUDIT_SHA256),
        adapter_source_sha256=mode.portable_audit.source_pins(),
        portable_proposal_links_resolved_from_theory=True,
        wrapper_limits_and_exact_proposal_checked=True,external_rows_equal_complete_inline_rows=True,
        source_and_input_bytes_checked_at_exit=True,
        external_input_sha256={str(Path(p['path']).relative_to(root)):p['sha256'] for p in files},
        weight_audits=[dict(receipt=summary,proposal=proposal,active_groups=913,even_rows=457,
            profile_pairs=208849,saved_log2_upper_rational='-59',
            all_producer_midpoint_and_final_counts_match=True,saved_ceiling_dominates_row_enclosures=True)])
    return write(root/'audit.json',report),report


class AdapterTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pins=mode.source_pins()
        # Cache immutable closure observations in synthetic mechanics tests;
        # repeated recursive source discovery is not the behavior under test.
        for module in (mode.midpoint,mode.midpoint.producer,mode.portable_audit):
            value=module.source_pins()
            cached=patch.object(module,'source_pins',side_effect=lambda value=value:dict(value))
            cached.start();cls.addClassCleanup(cached.stop)

    def setUp(self):
        self.temp=TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name);self.parent=parent_fixture(self.root/'old')
        self.audit,self.report=point_fixture(self.root)
        self.patches=[patch.object(mode,'DIRECTORY',self.root),
            patch.object(mode,'PARENT_SHA256',self.parent['sha256']),
            patch.object(mode,'source_pins',side_effect=lambda:dict(self.pins)),
            patch.object(mode.selected,'complete_selected',return_value=F(-59)),
            # The synthetic driver is outside the repository; the typed loader
            # still rehashes it. Production also checks the merged repository map.
            patch.object(mode.bands.frozen.producer,'check_source_pins')]
        for p in self.patches:p.start();self.addCleanup(p.stop)

    def extend(self,bands=()):
        return mode.extend_v25(self.parent['path'],self.audit['path'],bands)

    def replace_report(self,mutate):
        mutate(self.report);self.audit=write(Path(self.audit['path']),self.report)

    def test_selected_external_context_relative_proposal_and_parent_prefix(self):
        original=Path(self.parent['path']).read_bytes();prior=json.loads(original)
        result=self.extend()
        self.assertEqual(result['entries'][:249],prior['entries'])
        self.assertEqual(Path(self.parent['path']).read_bytes(),original)
        self.assertEqual((result['receipt_count'],result['covered_weight_count'],result['pending_weight_count']),(250,7745,447))
        self.assertEqual(result['added_weight_count'],1)
        self.assertEqual(result['entries'][249]['weight_range'],[913,913])
        self.assertFalse(result['whole_code_certificate']);self.assertFalse(result['coefficient_replay'])
        self.assertEqual(result['new_audit_schemas_allowed'],
            [*mode.bands.AUDIT_SCHEMAS,'gamma-midpoint-selected-receipts-audit-v1'])
        self.assertEqual(result['new_selected_audit_adapter_schemas_allowed'],
            ['gamma-r913-portable-midpoint-audit-adapter-v1'])
        self.assertEqual(result['covered_events_union_log2_upper_rational'],
            str(mode.index.union_ceiling(mode.bands.frozen._bounds(result['entries']),384)))
        mode.selected.complete_selected.assert_called_once()
        self.assertEqual(mode.selected.complete_selected.call_args.args[1],913)

    def test_optional_typed_band_delegation_and_overlap(self):
        pin=write(self.root/'band.json',dict(synthetic=True))
        event=dict(pin,weight_range=[7745,7750],weights_evaluated=6,log2_upper_rational='-100',zero=False)
        with patch.object(mode.bands,'_band_event',return_value=event) as validator:
            result=self.extend([pin['path']]);validator.assert_called_once()
        self.assertEqual((result['added_event_count'],result['added_weight_count']),(2,7))
        for interval in ([912,913],[913,913],[7745,7750]):
            event['weight_range']=interval
            with patch.object(mode.bands,'_band_event',return_value=event):
                with self.assertRaisesRegex(ValueError,'overlaps'):
                    self.extend([pin['path'],pin['path']])

    def test_public_count_limit_and_no_false_whole_code_claim(self):
        for invalid in ('audit.json',iter([]),[None]*64):
            with self.assertRaisesRegex(ValueError,'zero to63'):self.extend(invalid)
        calls=[]
        def event(path,*args):
            lo=int(path);hi=min(8191,lo+63)
            calls.append(lo)
            return dict(path=self.parent['path'],sha256=self.parent['sha256'],
                weight_range=[lo,hi],weights_evaluated=hi-lo+1,log2_upper_rational='-100')
        with patch.object(mode.bands,'_band_event',side_effect=event):
            result=self.extend(list(range(7745,8192,64)))
        self.assertEqual(result['covered_weight_count'],8192)
        self.assertEqual(result['pending_weight_ranges'],[])
        self.assertFalse(result['whole_code_certificate'])

    def test_missing_row_context_and_tampered_bytes_fail_closed(self):
        row_key=next(k for k in self.report['external_input_sha256'] if 'row_0001' in k)
        self.replace_report(lambda d:d['external_input_sha256'].pop(row_key))
        with self.assertRaisesRegex(ValueError,'Missing external checkpoint'):self.extend()
        self.audit,self.report=point_fixture(self.root)
        (self.root/'point'/'row_0000.json').write_text('{}')
        with self.assertRaisesRegex(ValueError,'hash mismatch'):self.extend()
        self.audit,self.report=point_fixture(self.root)
        self.replace_report(lambda d:d.update(context_sha256='0'*64))
        with self.assertRaisesRegex(ValueError,'digest'):self.extend()

    def test_wrong_audit_scope_proposal_and_source_conflicts(self):
        for mutation in (lambda d:d.update(selected_weights=[914]),
                lambda d:d.update(external_rows_checked=456),
                lambda d:d.update(coefficient_replay=True),
                lambda d:d.update(audit_adapter_schema='unreviewed'),
                lambda d:d.update(adapter_source_sha256={}),
                lambda d:d['weight_audits'][0].update(all_producer_midpoint_and_final_counts_match=False),
                lambda d:d.update(execution_source_sha256={})):
            self.audit,self.report=point_fixture(self.root);self.replace_report(mutation)
            with self.assertRaises(ValueError):self.extend()

    def test_prior_understatement_and_global_union_gate(self):
        with patch.object(mode,'PARENT_SHA256','0'*64):
            with self.assertRaisesRegex(ValueError,'hash mismatch'):self.extend()
        pin=write(self.root/'band.json',dict(synthetic=True))
        event=dict(pin,weight_range=[7745,7745],weights_evaluated=1,log2_upper_rational='-54')
        with patch.object(mode.bands,'_band_event',return_value=event):
            with self.assertRaisesRegex(ValueError,'Covered-event union'):self.extend([pin['path']])
        bad=parent_fixture(self.root/'bad',lambda d:d.update(covered_events_union_log2_upper_rational='-100'))
        self.parent=bad
        with patch.object(mode,'PARENT_SHA256',bad['sha256']):
            with self.assertRaisesRegex(ValueError,'Prior union'):self.extend()


if __name__=='__main__':
    unittest.main()
