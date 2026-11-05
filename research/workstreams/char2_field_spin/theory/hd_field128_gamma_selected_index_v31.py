"""Return-only partial index: frozen selected audits plus explicit composite audits.

Remote count/tree/arithmetic and archive checks are inherited, not replayed.
Even full weight coverage remains whole_code_certificate=False here.
"""
from copy import deepcopy
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path, PurePosixPath

import hd_field128_gamma_selected_index_v30 as prior
import hd_field128_gamma_composite_weight_band_audit_v30 as composite

HERE, NAME = Path(__file__).resolve().parent, Path(__file__).name
index, frozen = prior.index, prior.frozen
require, same, canonical = prior.require, prior.same, prior.canonical
PARENT_SHA = '875c572f4ae84b5bc139a4c21be9b0b4e18aed51ab1f288c45c125a260aab88f'
CHECK_SHA = '22357f681f345f94f13d72394df478984f8756d5b68fbaf0db917d0b6e9b604a'
PRIOR_SHA = '6e914264b8b42c162c93f17d174a979b1b450b99cfa107b4161907e981edc78b'
COMPOSITE_SHA = '8e2449ba8758431a74d7d18b883468cdf7dff114100c4a37f73386be38c70af5'
ARCHIVE_SCHEMA = 'field128-durable-remote-composite-archive-v31'
ARCHIVE_ROOT = PurePosixPath('/home/prindal/spin-proof-archives/field128-v31-20261010')
PARENT_EVENTS, PARENT_WEIGHTS, UNION_BITS = 265, 8131, 40
CONSTRUCTION = prior.CONSTRUCTION


def composite_sources():
    return composite.source_pins(composite.producer.helpers())


def source_pins():
    for module, expected in ((prior, PRIOR_SHA), (composite, COMPOSITE_SHA)):
        require(sha256(Path(module.__file__).read_bytes()).hexdigest() == expected,
                'Reviewed index/auditor changed')
    pins = index.source_closure(NAME)
    for closure in (prior.source_pins(), composite_sources()):
        frozen._merge(pins, closure)
    return dict(sorted(pins.items()))


def composite_event(report, pins):
    expected = composite_sources()
    same(report, dict(schema=composite.SCHEMA, audit_result='PASS', construction=CONSTRUCTION,
        field_order=1 << 128, state_dimension=30, target_distance='2/5',
        complete_requested_band=True, covers_all_matching_profiles=True,
        integer_partition_checked=True, counts_recomputed=True,
        count_algorithm='frozen independent raw-midpoint rectangle prefix recurrence',
        every_original_weight_preparation_metadata_checked=True,
        exact_literal_bank_context_matches=True, generic_partition_flags_preserved=True,
        component_minima_recomputed=True, complete_component_scalar_enclosures_checked=True,
        component_outer_prefactors_already_included=True, cell_log_enclosures_checked=True,
        directed_union_precision=384, meets_requested_band_target=True,
        source_sha256=expected, producer_sha256=composite.PRODUCER_SHA,
        numerical_auditor_sha256=composite.NUMERIC_SHA,
        source_and_input_bytes_checked_at_exit=True, producer_json_budget_bytes=32 << 20,
        historical_origin_verified=False, coefficient_replay=False,
        normalized_regional_replay=False, cell_transfer_recomputed=False,
        no_weight_interpolation=True, whole_code_certificate=False), 'composite audit')
    lo, hi = index.interval(report['requested_weight_range'])
    require(hi <= 8191 and hi-lo+1 <= 64, 'Composite band outside 1..8191 or wider than 64')
    weights = list(range(lo, hi+1))
    require(report['requested_weights'] == weights, 'Missing or reordered original weights')
    histograms = sum(r//2-max(0, r-4096)+1 for r in weights)
    require(type(report['normalized_histograms']) is int and report['normalized_histograms'] == histograms,
            'Original-weight histogram domain differs')
    integer = prior.prior._integer
    count = integer(report['component_count'], 2, 4, 'component count')
    cells = integer(report['cells_checked'], 1, 4096, 'evaluated cells')
    leaves = integer(report['accepted_cells'], 1, cells, 'accepted cells')
    require(cells == 2*leaves-len({r % 2 for r in weights}), 'Incomplete binary partition')
    require(type(report['checkpoints_checked']) is int and report['checkpoints_checked'] == cells+1,
            'Wrong checkpoint count')
    winners = report['accepted_witness_counts']
    require(isinstance(winners, list) and len(winners) == count
        and all(type(n) is int and n >= 0 for n in winners) and sum(winners) == leaves,
        'Component winners omit accepted cells')
    integer(report['positive_count_zero_cells'], 0, cells, 'zero cells')
    target = integer(report['requested_failure_bits'], 40, 4096, 'band target')
    integer(report['cell_failure_bits'], target, 4096, 'cell target')
    integer(report['producer_json_bytes'], 1, 32 << 20, 'producer JSON bytes')
    prior.prior._digest(report['context_sha256'])
    checked = prior.pin_map(report, absolute=False)
    for key in ('input_summary', 'execution_driver_receipt', 'context_file'):
        prior.linked(report[key], checked)
    directory = PurePosixPath(report['input_summary']['path']).parent
    for key, name in (('input_summary', 'summary.json'), ('execution_driver_receipt', 'driver.json'),
                      ('context_file', 'context.json')):
        require(report[key]['path'] == str(directory/name), 'Wrong composite metadata layout')
    required = {str(directory/name) for name in ('summary.json', 'driver.json', 'context.json')}
    for i in range(count):
        required.update(str(directory/name) for name in (f'bank_{i:02d}.json', f'component_{i:02d}_context.json'))
    required.update(str(directory/f'node_{i:06d}.json') for i in range(cells+1))
    proposals = report['component_literal_sources']
    require(isinstance(proposals, list) and len(proposals) == count, 'Missing literal sources')
    source_names = set()
    for pin in proposals:
        prior.linked(pin, checked)
        require(pin['path'] not in source_names and pin['path'] not in required,
                'Duplicate or overlapping literal source')
        source_names.add(pin['path'])
    require(set(checked) == required | source_names, 'Incomplete or extra composite file pins')
    zero = report['zero']
    require(type(zero) is bool, 'Missing explicit zero flag')
    keys = ('log2_upper_rational', 'saved_band_log2_upper_rational', 'reaggregated_cell_log2_upper_rational')
    if zero:
        require(all(report[key] is None for key in keys), 'Zero event has a logarithmic ceiling')
        bound = None
    else:
        bound, saved, fresh = map(prior.fraction, (report[key] for key in keys))
        require(fresh <= saved == bound <= -target, 'Composite ceiling understated or target missed')
    frozen._merge(pins, expected)
    return lo, hi, bound


def composite_archive(receipt, report, audit_pin):
    same(receipt, dict(schema=ARCHIVE_SCHEMA, host=prior.prior.HOST, zip_crc_checked=True,
        all_member_sha256_checked=True, exact_manifest_membership_checked=True,
        all_audit_input_files_included=True, audit_source_closure_included=True,
        storage_scope='remote durable archive; no local tree import'), 'composite archive')
    for key, suffix in (('archive', '.zip'), ('manifest', '.json')):
        path = prior.prior._posix(receipt[key]['path'], absolute=True)
        require(ARCHIVE_ROOT in path.parents and path.suffix == suffix, 'Wrong composite archive root/type')
        prior.prior._digest(receipt[key]['sha256'])
    prior.prior._integer(receipt['archive']['bytes'], 1, 1 << 40, 'archive bytes')
    prior.prior._integer(receipt['file_count'], len(report['checked_file_pins'])+1, 100000, 'members')
    prior.prior._integer(receipt['logical_bytes'], audit_pin['bytes']+report['producer_json_bytes'],
                        1 << 40, 'logical bytes')
    require(isinstance(receipt['restore'], str) and receipt['restore'].strip(), 'Missing restore instructions')
    require(receipt['audit_source_closure_sha256'] == sha256(canonical(report['source_sha256'])).hexdigest(),
            'Archived source closure differs')
    summary = prior.prior._posix(receipt['relative_summary_path'], absolute=False)
    audit = prior.prior._posix(receipt['relative_audit_path'], absolute=False)
    require(summary == PurePosixPath(report['input_summary']['path']) and summary.name == 'summary.json'
        and audit != summary and audit.parent == summary.parent and audit.suffix == '.json',
        'Wrong archive summary/audit layout')
    for key, relative, digest in (('summary', summary, report['input_summary']['sha256']),
                                  ('audit', audit, audit_pin['sha256'])):
        require(receipt[key]['member'] == 'theory/'+str(relative)
            and receipt[key]['sha256'] == digest, 'Archived bytes/member differ')


def extend_v30(parent_path, parent_check_path, event_specs):
    require(isinstance(event_specs, (list, tuple)) and 1 <= len(event_specs) <= 64,
            'Expected 1..64 explicit audit/archive pairs')
    loader, pins = index.Validator(F(2, 5)), source_pins()
    parent, parent_pin = prior.prior._load(loader, parent_path, PARENT_SHA, 'byte-pinned-v30-index')
    check, check_pin = prior.prior._load(loader, parent_check_path, CHECK_SHA, 'byte-pinned-v30-check')
    same(parent, dict(schema='field128-gamma-partial-cover-index-v30-selected', construction=CONSTRUCTION,
        field_order=1 << 128, source_dimension=1 << 20, output_length=1 << 21, target_distance='2/5',
        receipt_count=PARENT_EVENTS, covered_weight_count=PARENT_WEIGHTS, pending_weight_count=61,
        covered_events_failure_bits=40, meets_covered_events_target=True,
        directed_aggregation=True, whole_code_certificate=False), 'frozen parent')
    same(check, dict(schema='gamma-saved-partial-index-independent-check-v30', audit_result='PASS',
        construction=CONSTRUCTION, state_dimension=30, target_distance='2/5',
        full_saved_object_matches_validator_rebuild=True, parent_entry_count=260,
        parent_canonical_entry_bytes_preserved=True, parent_entries_and_order_preserved=True,
        receipt_count=PARENT_EVENTS, covered_weight_count=PARENT_WEIGHTS, pending_weight_count=61,
        nonoverlapping_weight_ranges=True, independent_positive_union_precision_bits=384,
        covered_events_failure_bits=40, meets_covered_events_target=True, whole_code_certificate=False,
        covered_events_union_log2_upper_rational=parent['covered_events_union_log2_upper_rational']), 'parent check')
    require(check['index']['sha256'] == PARENT_SHA
        and Path(check['index']['path']).name == Path(parent_pin['path']).name, 'Parent check binding differs')
    frozen._merge(pins, parent['aggregation_source_sha256'])
    entries, covered = deepcopy(parent['entries']), set()
    require(len(entries) == PARENT_EVENTS, 'Wrong parent entry count')
    for entry in entries:
        lo, hi = index.interval(entry['weight_range'])
        weights = set(range(lo, hi+1))
        require(not covered & weights, 'Overlapping parent entries')
        covered |= weights
    same(parent, dict(covered_weight_ranges=index.ranges(covered), covered_weight_count=len(covered),
        pending_weight_ranges=index.ranges(set(range(1, 8193))-covered)), 'parent domain')
    require(check['pending_weight_ranges'] == parent['pending_weight_ranges'], 'Parent pending check differs')
    require(index.union_ceiling(frozen._bounds(entries), 384) <= prior.fraction(parent['covered_events_union_log2_upper_rational']),
            'Parent union understated')
    added = []
    for spec in event_specs:
        require(isinstance(spec, dict) and set(spec) == {'audit', 'archive_receipt'}, 'Expected audit/archive pair')
        for pin in spec.values():
            require(isinstance(pin, dict) and set(pin) == {'path', 'sha256'}, 'Expected explicit path/SHA256')
        report, audit_pin = prior.prior._load(loader, spec['audit']['path'], spec['audit']['sha256'], 'remote-audit')
        archive, archive_pin = prior.prior._load(loader, spec['archive_receipt']['path'],
                                              spec['archive_receipt']['sha256'], 'durable-archive')
        is_composite = report.get('schema') == composite.SCHEMA
        if is_composite:
            lo, hi, bound = composite_event(report, pins)
            composite_archive(archive, report, audit_pin)
            kind = 'remote-literal-composite-band'
        else:
            selected, lo, bound = prior.selected_event(report, pins)
            prior.archive_receipt(archive, report, audit_pin, selected)
            hi, kind = lo, 'remote-literal-selected-'+selected
        weights = set(range(lo, hi+1))
        require(not covered & weights, 'New event overlaps prior/new coverage')
        covered |= weights
        entry = dict(audit_pin, audit=audit_pin, receipt_kind=kind, weight_range=[lo, hi],
            weights_evaluated=len(weights), zero=bound is None,
            log2_upper_rational=None if bound is None else str(bound),
            requested_failure_bits=report['requested_failure_bits'], actual_requested_target_met=True,
            remote_summary=report['input_summary'], remote_evidence_archive_receipt=archive_pin,
            archive=archive['archive'], manifest=archive['manifest'],
            audit_archive_member=archive['audit'], summary_archive_member=archive['summary'],
            context_sha256=report['context_sha256'], historical_origin_verified=False,
            mathematical_audit_inherited_from_remote=True, archive_verification_inherited_from_remote=True,
            local_tree_rechecked=False, local_archive_revalidated=False)
        entries.append(entry); added.append(entry)
    require(canonical(entries[:PARENT_EVENTS]) == canonical(parent['entries']), 'Parent entry prefix changed')
    union = index.union_ceiling(frozen._bounds(entries), 384)
    require(union <= -UNION_BITS, 'Actual covered-event union misses 40 bits')
    bounds = frozen._bounds(added)
    result = deepcopy(parent)
    result.update(schema='field128-gamma-partial-cover-index-v31-composite-selected', parent_index=parent_pin,
        parent_independent_check=check_pin, entries=entries, receipt_count=len(entries),
        covered_weight_ranges=index.ranges(covered), covered_weight_count=len(covered),
        pending_weight_ranges=index.ranges(set(range(1, 8193))-covered), pending_weight_count=8192-len(covered),
        covered_events_union_log2_upper_rational=str(union), aggregation_precision_bits=384,
        covered_events_failure_bits=40, parent_covered_events_failure_bits=40,
        meets_covered_events_target=True, whole_code_certificate=False,
        added_event_count=len(added), added_weight_count=sum(row['weights_evaluated'] for row in added),
        added_events_zero=not bounds,
        added_events_union_log2_upper_rational=str(index.union_ceiling(bounds, 384)) if bounds else None,
        latest_added_events=added,
        latest_explicit_selected_events=[row for row in added if row['receipt_kind'] != 'remote-literal-composite-band'],
        latest_audited_multiweight_events=[row for row in added if row['receipt_kind'] == 'remote-literal-composite-band'],
        prior_entry_records_and_order_preserved=True, completion_audit_performed=False,
        input_file_pins_verified=len(loader.files), checked_local_input_pins=list(loader.files.values()),
        aggregation_source_sha256=dict(sorted(pins.items())), helper_source_files_verified=len(pins),
        new_audit_schemas_allowed=[prior.single.SCHEMA, prior.multiple.SCHEMA, composite.SCHEMA],
        new_selected_audit_adapter_schemas_allowed=[prior.single.SCHEMA, prior.multiple.SCHEMA],
        provenance_scope='Inherited remote audits and durable archives; compact bytes/source pins checked locally',
        aggregation_method='Directed 384-bit union of actual ceilings on nonoverlapping weight ranges; target 40 bits',
        scope=f'{len(covered)} covered; {8192-len(covered)} pending; separate completion audit required')
    result['producer_receipt_counts'] = dict(parent['producer_receipt_counts'])
    for row in added:
        kind = row['receipt_kind']
        result['producer_receipt_counts'][kind] = result['producer_receipt_counts'].get(kind, 0)+1
    prior.prior.prior.band._check_links(tuple(loader.files.values()))
    frozen.producer.check_source_pins(pins)
    return result
