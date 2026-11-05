"""Return-only index. Remote mathematical/archive checks are inherited, not replayed.
Checks compact JSON/source bytes locally; no network, trees, coefficients or whole-code claim.
"""
from copy import deepcopy
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path, PurePosixPath

import hd_field128_gamma_selected_index_v27 as prior
import hd_field128_gamma_literal_remote_band_audit_v27 as audited

index, frozen = prior.index, prior.frozen
require, same = prior.require, prior.same
HERE = Path(__file__).resolve().parent
PARENT_SHA256 = 'bb5873cc96c861d3c9f08806bdbae63c5778fd9d3d57999a2631f529dd5e3a5f'
PARENT_CHECK_SHA256 = '2bd521821963ef03eeeec5ad55ff0945b0a6ccb90ca616e0842e4dda8e80dead'
PRIOR_SHA256 = 'db4d2621547fde7f33ae2aa88a3b360304e12c14778df46751eea59b5ded2c47'
AUDITOR_SHA256 = '7ea96c6f607b198b98aa56883fa46346ed46e6b10cdeccf85ab839e5fc3c511f'
PARENT_EVENTS, PARENT_WEIGHTS = 255, 7967
ARCHIVE_SCHEMA = 'field128-durable-remote-band-archive-v1'
ARCHIVE_ROOT = PurePosixPath('/home/prindal/spin-proof-archives/field128-v28-20261010')
HOST = 'peach48.devcore4.com:9022'
MAX_LOCAL_FILE_BYTES, MAX_LOCAL_INPUT_BYTES = 1 << 20, 16 << 20


def source_pins():
    for module, digest in ((prior, PRIOR_SHA256), (audited, AUDITOR_SHA256)):
        require(sha256(Path(module.__file__).read_bytes()).hexdigest() == digest,
                'Frozen index/auditor identity changed')
    pins = index.source_closure(Path(__file__).name)
    for closure in (prior.source_pins(), audited.source_pins()):
        frozen._merge(pins, closure)
    return dict(sorted(pins.items()))


def _digest(value):
    require(isinstance(value, str) and len(value) == 64
            and all(c in '0123456789abcdef' for c in value), 'Expected lowercase SHA256')
    return value


def _integer(value, lo, hi, label):
    require(type(value) is int and lo <= value <= hi, 'Invalid '+label)
    return value


def _posix(value, *, absolute):
    require(isinstance(value, str) and value and '\\' not in value and ':' not in value,
            'Expected a canonical POSIX evidence path')
    path = PurePosixPath(value)
    require(path.is_absolute() is absolute and str(path) == value
            and '..' not in path.parts and '.' not in path.parts, 'Unsafe evidence path')
    return path


def _remote_pin(value):
    require(isinstance(value, dict) and {'path', 'sha256'} <= value.keys(), 'Missing remote file pin')
    _posix(value['path'], absolute=True)
    _digest(value['sha256'])
    return value['path'], value['sha256']


def _load(loader, path, digest, role):
    path = Path(path)
    path = path if path.is_absolute() else HERE/path
    require(path.is_relative_to(HERE) and '..' not in path.parts, 'JSON input outside theory')
    path = path.resolve()
    require(path.is_relative_to(HERE), 'JSON input resolves outside theory')
    require(path.suffix.lower() == '.json' and path.stat().st_size <= MAX_LOCAL_FILE_BYTES,
            'Expected compact local JSON (at most one MiB)')
    value, pin = loader.load(path, expected=_digest(digest), role=role)
    require(sum(row['bytes'] for row in loader.files.values()) <= MAX_LOCAL_INPUT_BYTES,
            'Compact local inputs exceed 16 MiB')
    return value, pin


def _literal_audit(report, pins):
    same(report, dict(schema=audited.SCHEMA, audit_result='PASS', construction=audited.typed.CONSTRUCTION,
        state_dimension=30, field_order=1 << 128, target_distance='2/5',
        complete_requested_band=True, covers_all_matching_profiles=True, integer_partition_checked=True,
        every_original_weight_preparation_metadata_checked=True, counts_recomputed=True,
        count_algorithm='unchanged frozen independent midpoint prefix recurrence',
        cell_log_enclosures_checked=True, directed_union_precision=384,
        source_sha256=audited.source_pins(), frozen_typed_v2_auditor_sha256=audited.TYPED_SHA256,
        frozen_numerical_auditor_sha256=audited.typed.PREVIOUS_SHA256,
        exact_seed_literals_checked=True, exact_driver_record_checked=True, fresh_bank_metadata_checked=True,
        source_and_input_bytes_checked_at_exit=True, source_pin_conflicts=False,
        historical_origin_verified=False, calibration_history_replayed=False, historical_bounds_used=False,
        coefficient_replay_in_this_audit=False, normalized_regional_replay=False,
        cell_transfer_recomputed=False, no_weight_interpolation=True, whole_code_certificate=False),
        'explicit literal-only remote mathematical audit')
    lo, hi = index.interval(report['requested_weight_range'])
    require(hi <= 8191 and hi-lo+1 <= 64, 'Remote band is outside 1..8191 or wider than 64')
    weights = list(range(lo, hi+1))
    require(report['requested_weights'] == weights, 'Remote audit omits or reorders integer weights')
    histograms = sum(r//2-max(0, r-4096)+1 for r in weights)
    require(type(report['normalized_histograms']) is int and report['normalized_histograms'] == histograms,
            'Original-weight normalized histogram count differs')
    cells = _integer(report['cells_checked'], 1, 8192, 'remote cell count')
    leaves = _integer(report['accepted_cells'], 1, cells, 'accepted cell count')
    require(cells == 2*leaves-len({r % 2 for r in weights}), 'Complete binary partition size differs')
    require(type(report['checkpoints_checked']) is int and report['checkpoints_checked'] == cells+1,
            'Remote checkpoint count differs')
    _integer(report['positive_count_zero_cells'], 0, cells, 'positive-count zero cells')
    _integer(report['source_calibration_weight'], 1, 8192, 'seed calibration weight')
    _digest(report['context_sha256'])
    require(report['producer_json_budget_bytes'] == 32 << 20, 'Wrong producer output budget')
    _integer(report['producer_json_bytes'], 1, 32 << 20, 'producer JSON bytes')
    checked = {}
    for row in report['checked_file_pins']:
        path, digest = _remote_pin(row)
        require(path not in checked, 'Duplicate remote checked input')
        checked[path] = digest
    require(len(checked) == report['checkpoints_checked']+5, 'Remote checked-input set is incomplete')
    for key in ('input_summary', 'context_file', 'execution_driver_receipt', 'bank_metadata', 'portable_seed'):
        path, digest = _remote_pin(report[key])
        require(checked.get(path) == digest, 'Required remote file absent from audit pins: '+key)
    directory = _posix(report['input_summary']['path'], absolute=True).parent
    for key, name in (('input_summary','summary.json'), ('context_file','context.json'),
                      ('execution_driver_receipt','driver.json'), ('bank_metadata','bank.json')):
        require(report[key]['path'] == str(directory/name), 'Remote artifact layout differs: '+key)
    for i in range(report['checkpoints_checked']):
        require(str(directory/f'node_{i:06d}.json') in checked, 'Missing remote checkpoint pin')
    zero = report['zero']
    require(type(zero) is bool, 'Invalid remote zero flag')
    fields = ('log2_upper_rational', 'saved_band_log2_upper_rational', 'reaggregated_cell_log2_upper_rational')
    if zero:
        require(all(report[key] is None for key in fields)
                and report['reaggregated_cell_log2_upper_enclosure'] is None, 'Zero audit has a finite logarithm')
        bound = None
    else:
        values = []
        for key in fields:
            value = report[key]
            require(isinstance(value, str) and str(F(value)) == value, 'Noncanonical audit bound')
            values.append(F(value))
        bound, saved, reaggregated = values
        require(bound == saved and reaggregated <= saved <= -54, 'Remote band misses 54 bits or understates its union')
        require(isinstance(report['reaggregated_cell_log2_upper_enclosure'], str)
                and report['reaggregated_cell_log2_upper_enclosure'], 'Missing remote union enclosure')
    frozen._merge(pins, report['source_sha256'])
    require(isinstance(report['checked_source_sha256'], dict), 'Missing remote checked source closure')
    for name, digest in report['source_sha256'].items():
        require(report['checked_source_sha256'].get(name) == digest, 'Incomplete checked source closure')
    for name, digest in report['checked_source_sha256'].items():
        _digest(digest)
        # Only reviewed closure names may include ../ proof notes.
        if name not in report['source_sha256']:
            path = _posix(name, absolute=False)
            require(len(path.parts) == 1 and path.suffix == '.py', 'Invalid extra checked source path')
    frozen._merge(pins, report['checked_source_sha256'])
    return lo, hi, bound


def _archive(receipt, report, audit_pin):
    same(receipt, dict(schema=ARCHIVE_SCHEMA, host=HOST, zip_crc_checked=True,
        all_member_sha256_checked=True, exact_manifest_membership_checked=True,
        storage_scope='remote durable archive; no local tree import'), 'durable remote archive receipt')
    for key, suffix in (('archive','.zip'), ('manifest','.json')):
        path = _posix(receipt[key]['path'], absolute=True)
        require(ARCHIVE_ROOT in path.parents and path.suffix == suffix, 'Archive is outside the reviewed durable directory')
        _digest(receipt[key]['sha256'])
    _integer(receipt['archive']['bytes'], 1, 1 << 40, 'archive bytes')
    _integer(receipt['file_count'], len(report['checked_file_pins'])+1, 100000, 'archive file count')
    _integer(receipt['logical_bytes'], report['producer_json_bytes']+audit_pin['bytes'], 1 << 40, 'archive logical bytes')
    require(isinstance(receipt['restore'], str) and receipt['restore'].strip(), 'Missing restore instructions')
    members = {}
    for key in ('audit', 'summary'):
        member = _posix(receipt[key]['member'], absolute=False)
        require(member.parts[0] == 'theory' and member.suffix == '.json', 'Expected a theory-relative JSON member')
        _digest(receipt[key]['sha256'])
        members[key] = member
    require(members['audit'] != members['summary'] and members['audit'].parent == members['summary'].parent,
            'Audit and summary archive members are not sibling artifacts')
    require(receipt['audit']['sha256'] == audit_pin['sha256'], 'Archive audit member differs from explicit local audit')
    remote_summary = _posix(report['input_summary']['path'], absolute=True)
    require(members['summary'].name == 'summary.json'
            and remote_summary.parts[-len(members['summary'].parts):] == members['summary'].parts
            and receipt['summary']['sha256'] == report['input_summary']['sha256'],
            'Archive summary member differs from the audited summary')


def extend_v27(parent_path, parent_check_path, event_specs):
    """Return an extension for 1..64 pinned remote audits plus archive receipts."""
    require(isinstance(event_specs, (list, tuple)) and 1 <= len(event_specs) <= 64,
            'Expected one to 64 explicit remote audit/archive pairs')
    loader, pins = index.Validator(F(2, 5)), source_pins()
    parent, parent_pin = _load(loader, parent_path, PARENT_SHA256, 'byte-pinned-v27-index')
    check, check_pin = _load(loader, parent_check_path, PARENT_CHECK_SHA256, 'byte-pinned-v27-independent-check')
    same(parent, dict(schema='field128-gamma-partial-cover-index-v27', construction=audited.typed.CONSTRUCTION,
        field_order=1 << 128, source_dimension=1 << 20, output_length=1 << 21, target_distance='2/5',
        receipt_count=PARENT_EVENTS, covered_weight_count=PARENT_WEIGHTS, pending_weight_count=8192-PARENT_WEIGHTS,
        whole_code_certificate=False, directed_aggregation=True), 'frozen v27 index')
    same(check, dict(schema='gamma-saved-partial-index-independent-check-v27', audit_result='PASS',
        construction=audited.typed.CONSTRUCTION, state_dimension=30, target_distance='2/5',
        full_saved_object_matches_validator_rebuild=True, receipt_count=PARENT_EVENTS,
        covered_weight_count=PARENT_WEIGHTS, pending_weight_count=8192-PARENT_WEIGHTS,
        nonoverlapping_weight_ranges=True, independent_positive_union_precision_bits=384,
        covered_events_failure_bits=54, meets_covered_events_target=True, whole_code_certificate=False,
        coefficient_replay=False, normalized_regional_replay=False, transfer_recomputed=False,
        covered_events_union_log2_upper_rational=parent['covered_events_union_log2_upper_rational']), 'v27 independent check')
    require(check['index']['sha256'] == PARENT_SHA256
            and Path(check['index']['path']).name == Path(parent_pin['path']).name, 'Parent audit points to another index')
    frozen._merge(pins, parent['aggregation_source_sha256'])
    entries, covered = deepcopy(parent['entries']), set()
    require(len(entries) == PARENT_EVENTS, 'Wrong parent entry count')
    for entry in entries:
        lo, hi = index.interval(entry['weight_range'])
        weights = set(range(lo, hi+1))
        require(not covered & weights, 'Parent ranges overlap')
        covered |= weights
    same(parent, dict(covered_weight_ranges=index.ranges(covered), covered_weight_count=len(covered),
        pending_weight_ranges=index.ranges(set(range(1,8193))-covered)), 'parent coverage')
    require(index.union_ceiling(frozen._bounds(entries), 384) <= F(parent['covered_events_union_log2_upper_rational']),
            'Parent union understated')
    added = []
    for spec in event_specs:
        require(isinstance(spec, dict) and set(spec) == {'audit', 'archive_receipt'}, 'Expected explicit audit/archive pins')
        for key in spec:
            require(isinstance(spec[key], dict) and set(spec[key]) == {'path','sha256'}, 'Expected local path and SHA256')
        report, audit_pin = _load(loader, spec['audit']['path'], spec['audit']['sha256'], 'explicit-remote-literal-audit')
        receipt, archive_pin = _load(loader, spec['archive_receipt']['path'], spec['archive_receipt']['sha256'], 'explicit-durable-archive-receipt')
        lo, hi, bound = _literal_audit(report, pins)
        _archive(receipt, report, audit_pin)
        weights = set(range(lo, hi+1))
        require(not covered & weights, 'New remote band overlaps prior or newly added coverage')
        covered |= weights
        entry = dict(audit_pin, receipt_kind='remote-literal-audited-midpoint-band', audit=audit_pin,
            weight_range=[lo,hi], weights_evaluated=hi-lo+1, zero=bound is None,
            log2_upper_rational=None if bound is None else str(bound), remote_summary=report['input_summary'],
            remote_evidence_archive_receipt=archive_pin, archive=receipt['archive'], manifest=receipt['manifest'],
            audit_archive_member=receipt['audit'], summary_archive_member=receipt['summary'],
            context_sha256=report['context_sha256'], historical_origin_verified=False,
            mathematical_audit_inherited_from_remote=True, archive_verification_inherited_from_remote=True,
            local_tree_rechecked=False, local_archive_revalidated=False)
        entries.append(entry); added.append(entry)
    require(entries[:PARENT_EVENTS] == parent['entries'], 'Parent entry prefix changed')
    union = index.union_ceiling(frozen._bounds(entries), 384)
    require(union <= -54, 'Covered-event union fails 54-bit budget')
    added_bounds = frozen._bounds(added)
    result = deepcopy(parent)
    result.pop('all_receipt_bytes_match', None)
    result.pop('covered_events_union_log2_upper_enclosure', None)
    result.update(schema='field128-gamma-partial-cover-index-v28-remote-evidence',
        parent_index=parent_pin, parent_independent_check=check_pin, entries=entries, receipt_count=len(entries),
        covered_weight_ranges=index.ranges(covered), covered_weight_count=len(covered),
        pending_weight_ranges=index.ranges(set(range(1,8193))-covered), pending_weight_count=8192-len(covered),
        covered_events_union_log2_upper_rational=str(union), aggregation_precision_bits=384,
        covered_events_failure_bits=54, meets_covered_events_target=True, whole_code_certificate=False,
        added_event_count=len(added), added_weight_count=sum(row['weights_evaluated'] for row in added),
        added_events_zero=not added_bounds,
        added_events_union_log2_upper_rational=str(index.union_ceiling(added_bounds,384)) if added_bounds else None,
        latest_added_events=added, latest_audited_multiweight_events=added, latest_explicit_selected_events=[],
        prior_entry_records_and_order_preserved=True, prior_semantic_audits_inherited=True,
        prior_entry_bytes_rechecked_locally=False, prior_non_json_source_pins_rechecked=False,
        new_count_audits_inherited=True, new_mathematical_audits_inherited_from_remote=True,
        new_archive_verifications_inherited_from_remote=True, new_integer_partitions_rechecked=False,
        new_cell_unions_reaggregated=False, new_selected_rows_and_prefix_reaggregated=False,
        local_archive_revalidation=False, local_tree_replay=False, network_access_used=False,
        source_pin_authentication_claim=False, coefficient_replay=False, normalized_regional_replay=False,
        transfer_recomputed=False, input_file_pins_verified=len(loader.files),
        non_json_input_file_pins_verified=0, prior_non_json_input_records_inherited=True,
        checked_local_input_pins=list(loader.files.values()), all_local_input_bytes_match=True,
        aggregation_source_sha256=dict(sorted(pins.items())), helper_source_files_verified=len(pins),
        new_audit_schemas_allowed=[audited.SCHEMA], new_selected_audit_adapter_schemas_allowed=[],
        provenance_scope='Parent audit inherited; new remote literal audits and durable-archive receipts explicitly byte-pinned locally',
        aggregation_method=f'Directed384-bit union of {len(entries)} events on nonoverlapping weight ranges; new remote band audits inherited',
        scope=f'Only {len(covered)} listed weights; {8192-len(covered)} pending; no whole-code certificate')
    name = 'hd_field128_gamma_midpoint_weight_bands'
    result['producer_receipt_counts'][name] = result['producer_receipt_counts'].get(name,0)+len(added)
    prior.band._check_links(tuple(loader.files.values()))
    frozen.producer.check_source_pins(pins)
    return result
