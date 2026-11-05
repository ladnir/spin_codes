"""Memory-only contracts; no coefficients, files, or production-index generation."""
from contextlib import ExitStack
from hashlib import sha256
import json
import unittest
from unittest.mock import patch

import hd_field128_gamma_selected_index_v31 as mode
import hd_field128_gamma_selected_index_v30_tests as controls

COMPOSITE_SOURCES = {'composite.py': 'd'*64}


def encode(value):
    return (json.dumps(value, sort_keys=True, allow_nan=False)+'\n').encode()


def composite_report(lo, hi, bound='-51', target=48):
    weights = list(range(lo, hi+1))
    roots = len({r % 2 for r in weights})
    directory, files = f'band_{lo}_{hi}', {}
    def pin(name):
        path = directory+'/'+name
        row = dict(path=path, sha256=sha256(path.encode()).hexdigest(), bytes=10)
        files[path] = row
        return row
    summary, driver, context = pin('summary.json'), pin('driver.json'), pin('context.json')
    sources = []
    for i in range(2):
        pin(f'bank_{i:02d}.json'); pin(f'component_{i:02d}_context.json')
        sources.append(pin(f'proposal_{i}.json'))
    for i in range(roots+1):pin(f'node_{i:06d}.json')
    return dict(schema=mode.composite.SCHEMA, audit_result='PASS', construction=mode.CONSTRUCTION,
        field_order=1 << 128, state_dimension=30, target_distance='2/5',
        requested_weight_range=[lo, hi], requested_weights=weights,
        complete_requested_band=True, covers_all_matching_profiles=True,
        integer_partition_checked=True, counts_recomputed=True,
        count_algorithm='frozen independent raw-midpoint rectangle prefix recurrence',
        every_original_weight_preparation_metadata_checked=True,
        normalized_histograms=sum(r//2-max(0, r-4096)+1 for r in weights),
        exact_literal_bank_context_matches=True, generic_partition_flags_preserved=True,
        component_minima_recomputed=True, complete_component_scalar_enclosures_checked=True,
        component_outer_prefactors_already_included=True, cell_log_enclosures_checked=True,
        directed_union_precision=384, meets_requested_band_target=True,
        source_sha256=COMPOSITE_SOURCES, producer_sha256=mode.composite.PRODUCER_SHA,
        numerical_auditor_sha256=mode.composite.NUMERIC_SHA,
        source_and_input_bytes_checked_at_exit=True, producer_json_budget_bytes=32 << 20,
        producer_json_bytes=1000, historical_origin_verified=False, coefficient_replay=False,
        normalized_regional_replay=False, cell_transfer_recomputed=False,
        no_weight_interpolation=True, whole_code_certificate=False,
        component_count=2, cells_checked=roots, accepted_cells=roots, checkpoints_checked=roots+1,
        accepted_witness_counts=[roots, 0], positive_count_zero_cells=0,
        requested_failure_bits=target, cell_failure_bits=54, context_sha256='f'*64,
        checked_file_pins=list(files.values()), input_summary=summary,
        execution_driver_receipt=driver, context_file=context, component_literal_sources=sources,
        zero=bound is None, log2_upper_rational=bound, saved_band_log2_upper_rational=bound,
        reaggregated_cell_log2_upper_rational=bound)


class Fixture:
    def __init__(self):
        self.data, self.specs = {}, []
        self.parent_name = 'hd_field128_gamma30_d40_partial_cover_v30_20261010.json'
        self.check_name = 'hd_field128_gamma30_d40_partial_cover_v30_20261010_independent_check.json'
        for name in (self.parent_name, self.check_name):self.data[name] = (mode.HERE/name).read_bytes()
        self.parent = json.loads(self.data[self.parent_name])
        self.links_checked = self.sources_checked = False

    def put(self, name, value):
        raw = encode(value); self.data[name] = raw
        return dict(path=name, sha256=sha256(raw).hexdigest())

    def add(self, lo, hi=None, bound='-51', *, target=48, selected=False):
        hi = lo if hi is None else hi
        report = controls.audit('single', lo, bound, target) if selected else composite_report(lo, hi, bound, target)
        audit = self.put(f'compact/audit_{lo}_{hi}.json', report)
        directory = f'r{lo}' if selected else f'band_{lo}_{hi}'
        root = mode.prior.ARCHIVE_ROOT if selected else mode.ARCHIVE_ROOT
        schema = mode.prior.ARCHIVE_SCHEMA if selected else mode.ARCHIVE_SCHEMA
        archive = dict(schema=schema, host=mode.prior.prior.HOST, zip_crc_checked=True,
            all_member_sha256_checked=True, exact_manifest_membership_checked=True,
            all_audit_input_files_included=True, audit_source_closure_included=True,
            storage_scope='remote durable archive; no local tree import',
            archive=dict(path=str(root/f'{directory}.zip'), sha256='1'*64, bytes=100000),
            manifest=dict(path=str(root/f'{directory}.json'), sha256='2'*64),
            relative_summary_path=directory+'/summary.json', relative_audit_path=directory+'/independent_audit.json',
            summary=dict(member='theory/'+directory+'/summary.json', sha256=report['input_summary']['sha256']),
            audit=dict(member='theory/'+directory+'/independent_audit.json', sha256=audit['sha256']),
            file_count=len(report['checked_file_pins'])+20, logical_bytes=1000000, restore='Restore verified members.',
            audit_source_closure_sha256=sha256(mode.canonical(report['source_sha256'])).hexdigest())
        archive_pin = self.put(f'compact/archive_{lo}_{hi}.json', archive)
        spec = dict(audit=audit, archive_receipt=archive_pin)
        self.specs.append(spec)
        return spec

    def change(self, spec, key, change):
        pin = spec[key]; value = json.loads(self.data[pin['path']]); change(value)
        spec[key] = self.put(pin['path'], value)
        if key == 'audit':
            archive = spec['archive_receipt']; value = json.loads(self.data[archive['path']])
            value['audit']['sha256'] = spec['audit']['sha256']
            spec['archive_receipt'] = self.put(archive['path'], value)

    def run(self):
        fixture = self
        class Loader:
            def __init__(self, unused):self.files = {}
        def load(loader, path, expected, role):
            name = str(path).replace('\\', '/')
            if name.startswith(mode.HERE.as_posix()+'/'):name = name[len(mode.HERE.as_posix())+1:]
            raw = fixture.data[name]; digest = sha256(raw).hexdigest()
            mode.require(digest == expected, 'Memory input pin mismatch')
            pin = dict(path=name, sha256=digest, bytes=len(raw), roles=[role]); loader.files[name] = pin
            return json.loads(raw), pin
        def links(pins):
            for pin in pins:mode.require(sha256(fixture.data[pin['path']]).hexdigest() == pin['sha256'], 'Exit pin changed')
            fixture.links_checked = True
        def sources(pins):fixture.sources_checked = True
        with ExitStack() as stack:
            for obj, name, value in ((mode, 'source_pins', lambda:{}),
                (mode, 'composite_sources', lambda:COMPOSITE_SOURCES),
                (mode.prior, 'audit_sources', lambda kind:controls.SOURCES[kind]),
                (mode.prior.single.runner, 'source_pins', lambda:controls.EXECUTION),
                (mode.index, 'Validator', Loader), (mode.prior.prior, '_load', load),
                (mode.prior.prior.prior.band, '_check_links', links),
                (mode.frozen.producer, 'check_source_pins', sources)):
                stack.enter_context(patch.object(obj, name, value))
            return mode.extend_v30(self.parent_name, self.check_name, self.specs)


class IndexTests(unittest.TestCase):
    def test_full_domain_stays_partial_certificate_and_preserves_parent(self):
        f = Fixture(); f.add(521, 547); f.add(549, 575)
        for r in (917, 918, 919, 1146, 1147, 1148, 1149):f.add(r, selected=True)
        result = f.run()
        self.assertEqual((result['receipt_count'], result['covered_weight_count'], result['pending_weight_count']),
                         (274, 8192, 0))
        self.assertEqual(result['entries'][:265], f.parent['entries'])
        self.assertEqual(len(result['latest_audited_multiweight_events']), 2)
        self.assertEqual(len(result['latest_explicit_selected_events']), 7)
        self.assertEqual(result['added_weight_count'], 61)
        self.assertEqual(result['covered_events_failure_bits'], 40)
        self.assertFalse(result['whole_code_certificate'] or result['completion_audit_performed'])
        self.assertTrue(f.links_checked and f.sources_checked)

    def test_zero_event_and_arbitrary_partial_extension(self):
        f = Fixture(); f.add(525, 530, None)
        result = f.run()
        self.assertEqual(result['covered_weight_count'], 8137)
        self.assertTrue(result['added_events_zero'])
        self.assertIsNone(result['added_events_union_log2_upper_rational'])
        self.assertIsNone(result['entries'][-1]['log2_upper_rational'])
        self.assertFalse(result['whole_code_certificate'])

    def test_pin_overlap_unknown_schema_and_source_guard(self):
        f = Fixture(); f.add(521, 547); f.data[f.parent_name] += b' '
        with self.assertRaisesRegex(ValueError, 'pin mismatch'):f.run()
        f = Fixture(); f.add(521, 547); f.add(530, 535)
        with self.assertRaisesRegex(ValueError, 'overlap'):f.run()
        f = Fixture(); f.add(548)
        with self.assertRaisesRegex(ValueError, 'overlap'):f.run()
        f = Fixture(); spec = f.add(521, 547)
        f.change(spec, 'audit', lambda d:d.update(schema='invented'))
        with self.assertRaisesRegex(ValueError, 'schema'):f.run()
        with patch.object(mode, 'PRIOR_SHA', '0'*64):
            with self.assertRaisesRegex(ValueError, 'Reviewed'):mode.source_pins()

    def test_actual_union_not_individual_claim(self):
        f = Fixture(); f.add(521, 547, '-40', target=40)
        with self.assertRaisesRegex(ValueError, 'union misses'):f.run()

    def test_composite_metadata_count_partition_and_ceiling_rejections(self):
        changes = [lambda d:d.update(complete_requested_band=False),
            lambda d:d.update(whole_code_certificate=True), lambda d:d.update(requested_weights=[521]),
            lambda d:d.update(normalized_histograms=1), lambda d:d.update(cells_checked=1),
            lambda d:d.update(accepted_witness_counts=[0, 0]),
            lambda d:d.update(component_outer_prefactors_already_included=False),
            lambda d:d.update(component_minima_recomputed=False), lambda d:d.update(source_sha256={}),
            lambda d:d.update(requested_failure_bits=39), lambda d:d.update(cell_failure_bits=47),
            lambda d:d.update(reaggregated_cell_log2_upper_rational='-49'),
            lambda d:d.update(zero=True), lambda d:d['checked_file_pins'].pop(),
            lambda d:d['component_literal_sources'].reverse() or d['component_literal_sources'].append(
                d['component_literal_sources'][0])]
        for change in changes:
            f = Fixture(); spec = f.add(521, 547); f.change(spec, 'audit', change)
            with self.subTest(change=change):
                with self.assertRaises((ValueError, KeyError)):f.run()

    def test_archive_bindings_and_single_validator_retained(self):
        changes = [lambda d:d['archive'].update(path='/tmp/archive.zip'),
            lambda d:d.update(schema=mode.prior.ARCHIVE_SCHEMA),
            lambda d:d.update(relative_summary_path='../summary.json'),
            lambda d:d.update(all_audit_input_files_included=False),
            lambda d:d.update(audit_source_closure_sha256='0'*64),
            lambda d:d['summary'].update(sha256='0'*64),
            lambda d:d['audit'].update(member='theory/other/audit.json')]
        for change in changes:
            f = Fixture(); spec = f.add(549, 575); f.change(spec, 'archive_receipt', change)
            with self.assertRaises(ValueError):f.run()
        f = Fixture(); spec = f.add(917, selected=True)
        f.change(spec, 'audit', lambda d:d.update(context_projection_overwritten_keys=['scope', 'witness']))
        with self.assertRaises(ValueError):f.run()


if __name__ == '__main__':unittest.main()
