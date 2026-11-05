"""Memory contracts; no generated data."""
from copy import deepcopy
from hashlib import sha256
import json
from pathlib import Path
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import hd_field128_gamma_selected_index_v28 as mode


class RemoteIndexTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pins = mode.source_pins()
        cls.audit_sources = mode.audited.source_pins()

    def setUp(self):
        self.root = (mode.HERE/'__v28_memory_only__').resolve()
        self.files, self.read_counts = {}, {}
        self.change_on_close = None
        self.original_read, self.original_stat = Path.read_bytes, Path.stat
        self.original_is_file = Path.is_file
        for name, method in [('read_bytes', self.read), ('stat', self.stat), ('is_file', self.is_file)]:
            def dispatch(path, *args, method=method, **kwargs):
                return method(path, *args, **kwargs)
            p = patch.object(Path, name, dispatch)
            p.start(); self.addCleanup(p.stop)
        for owner, name, value in ((mode, 'source_pins', self.pins),
                                   (mode.audited, 'source_pins', self.audit_sources)):
            p = patch.object(owner, name, side_effect=lambda value=value: dict(value))
            p.start(); self.addCleanup(p.stop)
        entries = [dict(path=f'/not-local/old{r}.json', sha256='1'*64,
            weight_range=[r,r], log2_upper_rational='-100') for r in range(1,254)]
        entries += [dict(path='/not-local/rest.json', sha256='2'*64,
            weight_range=[254,7966], log2_upper_rational='-100'),
            dict(path='/not-local/terminal.json', sha256='3'*64, weight_range=[8192,8192],
                 zero=True, log2_upper_rational=None)]
        self.parent = dict(schema='field128-gamma-partial-cover-index-v27',
            construction=mode.audited.typed.CONSTRUCTION, field_order=1 << 128,
            source_dimension=1 << 20, output_length=1 << 21, target_distance='2/5',
            whole_code_certificate=False, directed_aggregation=True, entries=entries,
            receipt_count=255, covered_weight_count=7967, pending_weight_count=225,
            covered_weight_ranges=[[1,7966],[8192,8192]], pending_weight_ranges=[[7967,8191]],
            covered_events_union_log2_upper_rational='-90', aggregation_source_sha256={},
            producer_receipt_counts={}, all_receipt_bytes_match=True,
            checked_non_json_file_pins=[dict(inherited=True)], non_json_input_file_pins_verified=5)
        self.reset_parent()

    def reset_parent(self):
        self.parent_pin = self.put('parent.json', self.parent)
        self.check = dict(schema='gamma-saved-partial-index-independent-check-v27', audit_result='PASS',
            construction=mode.audited.typed.CONSTRUCTION, state_dimension=30, target_distance='2/5',
            index=self.parent_pin, full_saved_object_matches_validator_rebuild=True,
            receipt_count=255, covered_weight_count=7967, pending_weight_count=225,
            nonoverlapping_weight_ranges=True, independent_positive_union_precision_bits=384,
            covered_events_failure_bits=54, meets_covered_events_target=True, whole_code_certificate=False,
            coefficient_replay=False, normalized_regional_replay=False, transfer_recomputed=False,
            covered_events_union_log2_upper_rational=self.parent['covered_events_union_log2_upper_rational'])
        self.check_pin = self.put('check.json', self.check)

    def read(self, path):
        key = str(path)
        if key not in self.files:
            self.assertFalse(key.startswith(('/not-local','/tmp')))
            return self.original_read(path)
        self.read_counts[key] = self.read_counts.get(key,0)+1
        raw = self.files[key]
        return raw+b' ' if key == self.change_on_close and self.read_counts[key] > 1 else raw

    def stat(self, path, *args, **kwargs):
        if str(path) in self.files:
            return SimpleNamespace(st_size=len(self.files[str(path)]))
        return self.original_stat(path, *args, **kwargs)

    def is_file(self, path):
        return str(path) in self.files or self.original_is_file(path)

    def put(self, name, value):
        path = self.root/name
        raw = (json.dumps(value, sort_keys=True)+'\n').encode()
        self.files[str(path)] = raw
        return dict(path=str(path), sha256=sha256(raw).hexdigest())

    def report(self, first=7967, last=7970, *, bound='-100', zero=False):
        directory = f'/tmp/bundle/theory/band_{first}_{last}'
        pin = lambda name: dict(path=f'{directory}/{name}', sha256=sha256(name.encode()).hexdigest())
        roots = len({r % 2 for r in range(first,last+1)})
        row = dict(schema=mode.audited.SCHEMA, audit_result='PASS',
            construction=mode.audited.typed.CONSTRUCTION, state_dimension=30, field_order=1 << 128,
            target_distance='2/5', requested_weight_range=[first,last],
            requested_weights=list(range(first,last+1)), complete_requested_band=True,
            covers_all_matching_profiles=True, integer_partition_checked=True,
            every_original_weight_preparation_metadata_checked=True,
            normalized_histograms=sum(r//2-max(0,r-4096)+1 for r in range(first,last+1)),
            counts_recomputed=True, count_algorithm='unchanged frozen independent midpoint prefix recurrence',
            cells_checked=roots, accepted_cells=roots, checkpoints_checked=roots+1,
            cell_log_enclosures_checked=True, positive_count_zero_cells=0,
            directed_union_precision=384, log2_upper_rational=None if zero else bound,
            saved_band_log2_upper_rational=None if zero else bound,
            reaggregated_cell_log2_upper_rational=None if zero else bound,
            reaggregated_cell_log2_upper_enclosure=None if zero else f'[{bound} +/- 0]', zero=zero,
            context_sha256='4'*64, source_sha256=dict(self.audit_sources),
            checked_source_sha256=dict(self.audit_sources), source_and_input_bytes_checked_at_exit=True,
            frozen_typed_v2_auditor_sha256=mode.audited.TYPED_SHA256,
            frozen_numerical_auditor_sha256=mode.audited.typed.PREVIOUS_SHA256,
            exact_seed_literals_checked=True, exact_driver_record_checked=True, fresh_bank_metadata_checked=True,
            producer_json_bytes=1000, producer_json_budget_bytes=32 << 20, source_calibration_weight=first,
            historical_origin_verified=False, calibration_history_replayed=False, historical_bounds_used=False,
            source_pin_conflicts=False, coefficient_replay_in_this_audit=False,
            normalized_regional_replay=False, cell_transfer_recomputed=False,
            no_weight_interpolation=True, whole_code_certificate=False)
        for key, name in [('input_summary','summary.json'), ('context_file','context.json'),
            ('execution_driver_receipt','driver.json'), ('bank_metadata','bank.json')]:
            row[key] = pin(name)
        row['portable_seed'] = dict(path='/tmp/bundle/theory/seed.json', sha256='5'*64)
        row['checked_file_pins'] = [deepcopy(row[key]) for key in ('input_summary','context_file',
            'execution_driver_receipt','bank_metadata','portable_seed')]
        row['checked_file_pins'] += [pin(f'node_{i:06d}.json') for i in range(roots+1)]
        return row

    def spec(self, report=None, mutate_archive=None):
        report = self.report() if report is None else report
        lo, hi = report['requested_weight_range']
        audit = self.put(f'audit{lo}_{hi}.json', report)
        base = f'theory/band_{lo}_{hi}'
        receipt = dict(schema=mode.ARCHIVE_SCHEMA, host=mode.HOST,
            archive=dict(path=str(mode.ARCHIVE_ROOT/f'{lo}_{hi}.zip'), sha256='6'*64, bytes=1000),
            manifest=dict(path=str(mode.ARCHIVE_ROOT/f'{lo}_{hi}.json'), sha256='7'*64),
            audit=dict(member=base+'/independent_literal_audit.json', sha256=audit['sha256']),
            summary=dict(member=base+'/summary.json', sha256=report['input_summary']['sha256']),
            file_count=200, logical_bytes=1000000, zip_crc_checked=True,
            all_member_sha256_checked=True, exact_manifest_membership_checked=True,
            storage_scope='remote durable archive; no local tree import', restore='Verify then restore into a new directory.')
        if mutate_archive:
            mutate_archive(receipt)
        return dict(audit=audit, archive_receipt=self.put(f'archive{lo}_{hi}.json',receipt))

    def extend(self, specs):
        with patch.object(mode, 'PARENT_SHA256', self.parent_pin['sha256']), \
             patch.object(mode, 'PARENT_CHECK_SHA256', self.check_pin['sha256']):
            return mode.extend_v27(self.parent_pin['path'], self.check_pin['path'], specs)

    def test_remote_contract_parent_prefix_and_fresh_union(self):
        spec = self.spec()
        result = self.extend([spec])
        self.assertEqual(result['entries'][:255], self.parent['entries'])
        self.assertEqual((result['receipt_count'], result['covered_weight_count'],result['pending_weight_count']),
                         (256,7971,221))
        self.assertEqual(result['covered_events_union_log2_upper_rational'],
                         str(mode.index.union_ceiling(mode.frozen._bounds(result['entries']),384)))
        for key in ('new_mathematical_audits_inherited_from_remote','new_archive_verifications_inherited_from_remote'):
            self.assertTrue(result[key])
        for key in ('whole_code_certificate','coefficient_replay','normalized_regional_replay','transfer_recomputed',
                    'local_tree_replay','network_access_used','local_archive_revalidation',
                    'new_integer_partitions_rechecked','new_cell_unions_reaggregated'):
            self.assertFalse(result[key])

    def test_complete_coverage_zero_and_limits(self):
        for invalid in ([],[None]*65,'audit.json',iter([])):
            with self.assertRaisesRegex(ValueError,'one to 64'):
                self.extend(invalid)
        specs = [self.spec(self.report(r,min(r+63,8191),zero=True)) for r in range(7967,8192,64)]
        result = self.extend(specs)
        self.assertEqual(result['covered_weight_count'],8192)
        self.assertFalse(result['whole_code_certificate'])

    def test_identity_and_required_remote_flags(self):
        for field, bad in [('schema','gamma-midpoint-weight-band-audit-v1'),('audit_result','FAIL'),
            ('historical_origin_verified',True),('exact_seed_literals_checked',False),
            ('integer_partition_checked',False),('counts_recomputed',False),('whole_code_certificate',True),
            ('frozen_typed_v2_auditor_sha256','0'*64)]:
            row=self.report(); row[field]=bad
            with self.subTest(field=field), self.assertRaises(ValueError):
                self.extend([self.spec(row)])

    def test_weight_geometry_and_remote_pin_set(self):
        mutations = [lambda r:r['requested_weights'].pop(),lambda r:r.update(normalized_histograms=0),
            lambda r:r.update(accepted_cells=1),lambda r:r.update(checkpoints_checked=1),
            lambda r:r['checked_file_pins'].pop(),lambda r:r['checked_file_pins'].append(r['portable_seed']),
            lambda r:r['context_file'].update(sha256='a'*64),
            lambda r:r['checked_file_pins'][-1].update(path='/tmp/../elsewhere.json'),
            lambda r:r.update(producer_json_bytes=33 << 20)]
        for mutate in mutations:
            row=self.report(); mutate(row)
            with self.assertRaises(ValueError):
                self.extend([self.spec(row)])
        with self.assertRaises(ValueError):
            self.extend([self.spec(self.report(7967,8031))])

    def test_archive_identity_paths_and_verification(self):
        mutations=[lambda r:r.update(host='otherhost'),lambda r:r.update(zip_crc_checked=False),
            lambda r:r.update(all_member_sha256_checked=False),
            lambda r:r.update(exact_manifest_membership_checked=False),
            lambda r:r['archive'].update(path='/tmp/archive.zip'),
            lambda r:r['audit'].update(sha256='a'*64),lambda r:r['summary'].update(sha256='b'*64),
            lambda r:r['summary'].update(member='theory/another/summary.json'),
            lambda r:r['audit'].update(member='theory/../audit.json'),lambda r:r.update(file_count=1)]
        for mutate in mutations:
            with self.assertRaises(ValueError):
                self.extend([self.spec(mutate_archive=mutate)])

    def test_overlap_global_union_and_understated_remote_union(self):
        spec=self.spec()
        for specs in ([spec,spec],[self.spec(self.report(7966,7967))]):
            with self.assertRaisesRegex(ValueError,'overlaps'):
                self.extend(specs)
        with self.assertRaisesRegex(ValueError,'54-bit budget'):
            self.extend([self.spec(self.report(bound='-54'))])
        row=self.report(); row['reaggregated_cell_log2_upper_rational']='-99'
        with self.assertRaises(ValueError):
            self.extend([self.spec(row)])

    def test_local_paths_before_io(self):
        outside=mode.HERE.parent/'outside.json'
        with patch.object(Path,'resolve',side_effect=AssertionError('I/O')):
            for path in (outside,'//server/share/a.json','../outside.json'):
                with self.assertRaises(ValueError):
                    mode._load(None,path,'0'*64,'test')
        with patch.object(Path,'resolve',return_value=outside):
            with self.assertRaises(ValueError):
                mode._load(None,mode.HERE/'link.json','0'*64,'test')

    def test_hash_tamper_closing_checks_and_bad_source_paths(self):
        spec=self.spec(); path=spec['audit']['path']
        self.files[path]+=b' '
        with self.assertRaises(ValueError):
            self.extend([spec])
        spec=self.spec(); self.read_counts.clear(); self.change_on_close=spec['archive_receipt']['path']
        with self.assertRaisesRegex(RuntimeError,'immutable band artifact changed'):
            self.extend([spec])
        self.change_on_close=None
        for name in ('../other.md','/tmp/helper.py'):
            row=self.report(); row['checked_source_sha256'][name]='a'*64
            with self.assertRaises(ValueError):
                self.extend([self.spec(row)])
        row=self.report(); row['checked_source_sha256'][next(iter(self.audit_sources))]='a'*64
        with self.assertRaises(ValueError):
            self.extend([self.spec(row)])


if __name__ == '__main__':
    unittest.main()
