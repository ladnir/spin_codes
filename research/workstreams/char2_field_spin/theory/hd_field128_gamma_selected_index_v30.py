"""Return-only selected-event extension; inherited remote audits, fresh40-bit union."""
from copy import deepcopy
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path, PurePosixPath

import hd_field128_gamma_selected_index_v28 as prior
import hd_field128_gamma_literal_midpoint_audit_v30 as single
import hd_field128_gamma_profile_min_audit_v29 as multiple

HERE, NAME = Path(__file__).resolve().parent, Path(__file__).name
index, frozen = prior.index, prior.frozen
require, same = prior.require, prior.same
PARENT_SHA = 'ea1900e2f869a258f16404bd34c3882675a78e960e2703a3754820cc4e63737d'
CHECK_SHA = '3a68ec34469cf68f381b4b40cdedf8e1c16220940f14507953845dc23c805056'
PRIOR_SHA = 'ff1d62fbeee9b232319cea920a26a7485bacb735ce0aa9ae59b006e2dc1246d2'
SINGLE_SHA = '75f57e3ac55f247edf9a9e6e1fdee884cbb4a4a770bc46d77e21cd7c79df8e21'
MULTIPLE_SHA = '86ec57966c0676297de6067b982b9595f5e1b065332cc13a00c529be99927d19'
ARCHIVE_SCHEMA = 'field128-durable-remote-selected-archive-v30'
ARCHIVE_ROOT = PurePosixPath('/home/prindal/spin-proof-archives/field128-v30-20261010')
CONSTRUCTION = prior.audited.typed.CONSTRUCTION
PARENT_EVENTS, PARENT_WEIGHTS, UNION_BITS = 260, 8126, 40


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()


def audit_sources(kind):
    return single.source_pins() if kind == 'single' else multiple.source_pins()[1]


def source_pins():
    for module, expected in ((prior, PRIOR_SHA), (single, SINGLE_SHA), (multiple, MULTIPLE_SHA)):
        require(sha256(Path(module.__file__).read_bytes()).hexdigest() == expected, 'Reviewed helper changed')
    pins = index.source_closure(NAME)
    for closure in (prior.source_pins(), audit_sources('single'), audit_sources('multiple')):
        frozen._merge(pins, closure)
    return dict(sorted(pins.items()))


def fraction(value):
    require(isinstance(value, str) and str(F(value)) == value, 'Expected canonical rational ceiling')
    return F(value)


def pin_map(report, absolute):
    result = {}
    require(isinstance(report['checked_file_pins'], list), 'Missing checked file pins')
    for item in report['checked_file_pins']:
        require(isinstance(item, dict), 'Invalid checked file pin')
        path = str(prior._posix(item['path'], absolute=absolute))
        digest = prior._digest(item['sha256'])
        require(path not in result, 'Duplicate remote checked file')
        result[path] = digest
    return result


def linked(pin, checked):
    require(isinstance(pin, dict) and checked.get(pin.get('path')) == prior._digest(pin.get('sha256')),
            'Required remote artifact absent from audit pins')


def selected_event(report, pins):
    kind = ('single' if report.get('schema') == single.SCHEMA else
            'multiple' if report.get('schema') == multiple.SCHEMA else None)
    require(kind is not None, 'Unreviewed selected audit schema')
    expected_sources = audit_sources(kind)
    same(report, dict(audit_result='PASS', construction=CONSTRUCTION, field_order=1 << 128,
        state_dimension=30, target_distance='2/5', complete_selected_weight=True,
        covers_all_matching_profiles=True, historical_origin_verified=False,
        normalized_regional_replay=False, whole_code_certificate=False,
        meets_requested_selected_weight_target=True, source_sha256=expected_sources), 'selected audit')
    target = prior._integer(report['requested_failure_bits'], 40, 4096, 'selected target')
    if kind == 'single':
        require(isinstance(report['selected_weights'], list) and len(report['selected_weights']) == 1,
                'Expected one selected weight')
        r = prior._integer(report['selected_weights'][0], 1, 8191, 'original weight')
    else:
        r = prior._integer(report['active_groups'], 2, 4095, 'profile-min original weight')
    first, last = max(0, r-4096), r//2
    rows, profiles = last-first+1, (last-first+1)**2
    checked = pin_map(report, absolute=kind == 'single')
    linked(report['input_summary'], checked)
    directory = PurePosixPath(report['input_summary']['path']).parent
    require(PurePosixPath(report['input_summary']['path']).name == 'summary.json', 'Wrong summary layout')
    if kind == 'single':
        same(report, dict(counts_recomputed=True, exact_seed_literals_checked=True,
            exact_driver_record_checked=True, fresh_bank_metadata_checked=True,
            bank_context_scope_checked=True, selected_receipt_scope_checked=True,
            context_projection_overwritten_keys=['scope'],
            external_checkpoint_context_reconstructed=True, external_rows_checked=rows,
            external_rows_equal_complete_inline_rows=True, precision_bits=384,
            numerical_auditor_sha256=single.FROZEN_SHA, execution_source_sha256=single.runner.source_pins(),
            source_and_input_bytes_checked_at_exit=True, producer_json_budget_bytes=2 << 20,
            calibration_history_replayed=False, historical_bounds_used=False,
            coefficient_replay=False, transfer_recomputed=False, coverage_index_modified=False), 'single audit')
        same(report['execution_driver'], dict(path=single.runner.NAME, sha256=single.RUNNER_SHA), 'execution driver')
        prior._integer(report['producer_json_bytes'], 1, 2 << 20, 'producer bytes')
        prior._integer(report['source_calibration_weight'], 1, 8192, 'calibration weight')
        require(len(report['weight_audits']) == 1, 'Missing selected count audit')
        row = report['weight_audits'][0]
        same(row, dict(active_groups=r, even_rows=rows, profile_pairs=profiles,
            all_producer_midpoint_and_final_counts_match=True, saved_ceiling_dominates_row_enclosures=True,
            requested_failure_bits=target, meets_requested_target=True), 'single count audit')
        prior._digest(row['exact_count_table_sha256'])
        bound = fraction(report['saved_log2_upper_rational'])
        require(fraction(row['saved_log2_upper_rational']) == bound
            and fraction(row['saved_row_reaggregation_log2_upper_rational']) <= bound,
            'Saved ceiling understates rows')
        require(index.union_ceiling([bound], 384) <= fraction(report['selected_events_union_log2_upper_rational']),
                'Union understated')
        for key, name in (('execution_driver_receipt', 'driver.json'), ('bank_metadata', 'bank.json'),
                          ('context_file', 'context.json')):
            linked(report[key], checked)
            require(report[key]['path'] == str(directory/name), 'Wrong single artifact layout')
        linked(report['portable_seed'], checked)
        require(len(checked) == rows+5, 'Missing single inputs')
        for full in range(first, last+1):
            require(str(directory/f'row_{full:04d}.json') in checked, 'Missing original-r row checkpoint')
    else:
        same(report, dict(producer_sha256=multiple.PRODUCER_SHA, even_rows=rows, profile_pairs=profiles,
            exact_literal_bank_context_matches=True, original_r_prefixes_recomputed=True,
            coefficient_replay_in_this_audit=False, transfer_or_winner_replay=False,
            row_partition_and_count_policy_checked=True, directed_union_precision=384,
            second_outer_prefix_applied=False), 'profile-min audit')
        wins = report['witness_winning_profile_counts']
        require(isinstance(wins, list) and 2 <= len(wins) <= 4
            and all(type(n) is int and n >= 0 for n in wins), 'Invalid profile witness counters')
        positives = prior._integer(report['positive_count_profiles'], 0, profiles, 'positive profiles')
        require(sum(wins) == positives, 'Profile witness counters omit positive profiles')
        fresh = report['fresh_saved_row_union_384']
        require(type(fresh['zero']) is bool, 'Missing profile-min zero flag')
        if fresh['zero']:
            require(report['log2_upper_rational'] is None and fresh['log2_upper_rational'] is None,
                    'Zero profile-min event has a logarithmic ceiling')
            bound = None
        else:
            bound = fraction(report['log2_upper_rational'])
            require(fraction(fresh['log2_upper_rational']) <= bound, 'Profile-min row union understated')
        require(len(checked) == rows+2*len(wins)+3, 'Missing profile-min inputs')
        for name in ['driver.json', 'context.json']+[f'bank_{i:02d}.json' for i in range(len(wins))]:
            require(str(directory/name) in checked, 'Missing profile-min metadata')
        for i in range(rows):
            require(str(directory/f'row_{i:06d}.json') in checked, 'Missing profile-min row checkpoint')
    prior._digest(report['context_sha256'])
    require(bound is None or bound <= -target, 'Selected target missed')
    frozen._merge(pins, expected_sources)
    return kind, r, bound


def archive_receipt(receipt, report, audit_pin, kind):
    same(receipt, dict(schema=ARCHIVE_SCHEMA, host=prior.HOST, zip_crc_checked=True,
        all_member_sha256_checked=True, exact_manifest_membership_checked=True,
        all_audit_input_files_included=True, audit_source_closure_included=True,
        storage_scope='remote durable archive; no local tree import'), 'selected archive receipt')
    for key, suffix in (('archive', '.zip'), ('manifest', '.json')):
        path = prior._posix(receipt[key]['path'], absolute=True)
        require(ARCHIVE_ROOT in path.parents and path.suffix == suffix, 'Unapproved durable archive root/type')
        prior._digest(receipt[key]['sha256'])
    prior._integer(receipt['archive']['bytes'], 1, 1 << 40, 'archive bytes')
    prior._integer(receipt['file_count'], len(report['checked_file_pins'])+1, 100000, 'archive members')
    prior._integer(receipt['logical_bytes'], audit_pin['bytes']+report.get('producer_json_bytes', 0),
                   1 << 40, 'archive logical bytes')
    require(isinstance(receipt['restore'], str) and receipt['restore'].strip(), 'Missing restore instructions')
    require(receipt['audit_source_closure_sha256'] == sha256(canonical(report['source_sha256'])).hexdigest(),
            'Closure differs')
    summary_rel = prior._posix(receipt['relative_summary_path'], absolute=False)
    audit_rel = prior._posix(receipt['relative_audit_path'], absolute=False)
    require(summary_rel.name == 'summary.json' and audit_rel.suffix == '.json'
        and audit_rel != summary_rel and audit_rel.parent == summary_rel.parent, 'Wrong archive sibling layout')
    for key, relative in (('summary', summary_rel), ('audit', audit_rel)):
        require(receipt[key]['member'] == 'theory/'+str(relative), 'Wrong theory archive member')
        prior._digest(receipt[key]['sha256'])
    require(receipt['summary']['sha256'] == report['input_summary']['sha256']
        and receipt['audit']['sha256'] == audit_pin['sha256'], 'Archive receipt has different audit/summary bytes')
    source = prior._posix(report['input_summary']['path'], absolute=kind == 'single')
    if kind == 'single':
        member = PurePosixPath(receipt['summary']['member'])
        require(source.parts[-len(member.parts):] == member.parts, 'Summary/member mismatch')
    else:
        require(source == summary_rel, 'Relative summary/member mismatch')


def extend_v29(parent_path, parent_check_path, event_specs):
    require(isinstance(event_specs, (list, tuple)) and 1 <= len(event_specs) <= 64,
            'Expected1..64 audit/archive pairs')
    loader, pins = index.Validator(F(2, 5)), source_pins()
    parent, parent_pin = prior._load(loader, parent_path, PARENT_SHA, 'byte-pinned-v29-index')
    check, check_pin = prior._load(loader, parent_check_path, CHECK_SHA, 'byte-pinned-v29-independent-check')
    same(parent, dict(schema='field128-gamma-partial-cover-index-v28-remote-evidence',
        construction=CONSTRUCTION, field_order=1 << 128, source_dimension=1 << 20,
        output_length=1 << 21, target_distance='2/5', receipt_count=PARENT_EVENTS,
        covered_weight_count=PARENT_WEIGHTS, pending_weight_count=66,
        covered_events_failure_bits=54, meets_covered_events_target=True,
        directed_aggregation=True, whole_code_certificate=False), 'frozen v29 index')
    same(check, dict(schema='gamma-saved-partial-index-independent-check-v29', audit_result='PASS',
        full_saved_object_matches_validator_rebuild=True, receipt_count=PARENT_EVENTS,
        covered_weight_count=PARENT_WEIGHTS, pending_weight_count=66,
        independent_positive_union_precision_bits=384, covered_events_failure_bits=54,
        meets_covered_events_target=True, whole_code_certificate=False,
        covered_events_union_log2_upper_rational=parent['covered_events_union_log2_upper_rational']), 'v29 audit')
    require(check['index']['sha256'] == PARENT_SHA
        and Path(check['index']['path']).name == Path(parent_pin['path']).name, 'Parent audit binding differs')
    frozen._merge(pins, parent['aggregation_source_sha256'])
    entries, covered = deepcopy(parent['entries']), set()
    require(len(entries) == PARENT_EVENTS, 'Wrong parent entry count')
    for entry in entries:
        lo, hi = index.interval(entry['weight_range']); weights = set(range(lo, hi+1))
        require(not covered & weights, 'Parent weight ranges overlap'); covered |= weights
    same(parent, dict(covered_weight_ranges=index.ranges(covered), covered_weight_count=len(covered),
        pending_weight_ranges=index.ranges(set(range(1,8193))-covered)), 'parent coverage')
    require(index.union_ceiling(frozen._bounds(entries), 384) <= F(parent['covered_events_union_log2_upper_rational']),
            'Parent union understated')
    added = []
    for spec in event_specs:
        require(isinstance(spec, dict) and set(spec) == {'audit', 'archive_receipt'}, 'Expected audit/archive pins')
        for pin in spec.values():
            require(isinstance(pin, dict) and set(pin) == {'path', 'sha256'}, 'Expected explicit path/SHA256')
        report, audit_pin = prior._load(loader, spec['audit']['path'], spec['audit']['sha256'], 'remote-selected-audit')
        archive, archive_pin = prior._load(loader, spec['archive_receipt']['path'], spec['archive_receipt']['sha256'],
                                         'durable-selected-archive')
        kind, r, bound = selected_event(report, pins)
        archive_receipt(archive, report, audit_pin, kind)
        require(r not in covered, 'Selected weight overlap'); covered.add(r)
        entry = dict(audit_pin, receipt_kind='remote-literal-selected-'+kind, audit=audit_pin,
            weight_range=[r,r], weights_evaluated=1, zero=bound is None,
            log2_upper_rational=None if bound is None else str(bound),
            requested_failure_bits=report['requested_failure_bits'], actual_requested_target_met=True,
            remote_summary=report['input_summary'], remote_evidence_archive_receipt=archive_pin,
            archive=archive['archive'], manifest=archive['manifest'],
            audit_archive_member=archive['audit'], summary_archive_member=archive['summary'],
            context_sha256=report['context_sha256'], historical_origin_verified=False,
            mathematical_audit_inherited_from_remote=True, archive_verification_inherited_from_remote=True,
            local_tree_rechecked=False, local_archive_revalidated=False)
        entries.append(entry); added.append(entry)
    require(canonical(entries[:PARENT_EVENTS]) == canonical(parent['entries']), 'Parent entries changed')
    union = index.union_ceiling(frozen._bounds(entries), 384)
    require(union <= -UNION_BITS, 'Covered union misses40 bits')
    new_bounds = frozen._bounds(added)
    result = deepcopy(parent)
    result.update(schema='field128-gamma-partial-cover-index-v30-selected', parent_index=parent_pin,
        parent_independent_check=check_pin, entries=entries, receipt_count=len(entries),
        covered_weight_ranges=index.ranges(covered), covered_weight_count=len(covered),
        pending_weight_ranges=index.ranges(set(range(1,8193))-covered), pending_weight_count=8192-len(covered),
        covered_events_union_log2_upper_rational=str(union),
        covered_events_failure_bits=UNION_BITS, meets_covered_events_target=True, whole_code_certificate=False,
        parent_covered_events_failure_bits=54, added_event_count=len(added), added_weight_count=len(added),
        added_events_zero=not new_bounds,
        added_events_union_log2_upper_rational=str(index.union_ceiling(new_bounds,384)) if new_bounds else None,
        latest_added_events=added, latest_explicit_selected_events=added, latest_audited_multiweight_events=[],
        input_file_pins_verified=len(loader.files), checked_local_input_pins=list(loader.files.values()),
        aggregation_source_sha256=dict(sorted(pins.items())), helper_source_files_verified=len(pins),
        new_audit_schemas_allowed=[single.SCHEMA,multiple.SCHEMA],
        new_selected_audit_adapter_schemas_allowed=[single.SCHEMA,multiple.SCHEMA],
        provenance_scope='Inherited parent/remote audits and archives; no tree replay',
        aggregation_method='Directed384 union of actual ceilings; global40 budget',
        scope=f'{len(covered)} covered; {8192-len(covered)} pending; not whole-code')
    result['producer_receipt_counts'] = dict(parent['producer_receipt_counts'])
    for kind in ('single','multiple'):
        name = 'remote-literal-selected-'+kind
        result['producer_receipt_counts'][name] = result['producer_receipt_counts'].get(name,0)+sum(
            row['receipt_kind'] == name for row in added)
    prior.prior.band._check_links(tuple(loader.files.values()))
    frozen.producer.check_source_pins(pins)
    return result
