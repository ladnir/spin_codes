"""Explicit typed-provenance extension of byte-pinned v21; no coefficient replay.

Legacy audits retain their frozen validator. Typed audits have a distinct
identity and recheck Python/Markdown byte links without parsing them as JSON.
The count audit is inherited; integer partitions and saved unions are rechecked.
"""

from copy import deepcopy
from fractions import Fraction as F
from hashlib import sha256
import json
import os
from pathlib import Path

from flint import arb, ctx

import hd_field128_gamma_selected_index_v22 as previous
import hd_field128_gamma_midpoint_weight_band_typed_audit as audited
import hd_field128_gamma_midpoint_weight_band_typed_v2_audit as audited_v2

frozen = previous.frozen
index, band, midpoint = previous.index, previous.band, previous.midpoint
producer = frozen.producer
require, same, fraction = previous.require, previous.same, previous.fraction
DIRECTORY = Path(__file__).resolve().parent
PARENT_SHA256 = previous.PARENT_SHA256
PARENT_EVENT_COUNT = previous.PARENT_EVENT_COUNT
FAILURE_BITS = previous.FAILURE_BITS


def source_pins():
    return dict(sorted({**index.source_closure(Path(__file__).name),
                        **previous.source_pins(), **audited.source_pins(),
                        **audited_v2.source_pins()}.items()))


def _typed_band_event(path, loader, pins, source_inputs, audit_provider, schema):
    report, audit_pin = loader.load(path, role='new-independent-weight-band-audit')
    same(report, dict(schema=schema,
        audit_result='PASS', complete_requested_band=True, covers_all_matching_profiles=True,
        integer_partition_checked=True, counts_recomputed=True, cell_log_enclosures_checked=True,
        directed_union_precision=384, coefficient_replay_in_this_audit=False,
        normalized_regional_replay=False, source_sha256=audit_provider.source_pins(),
        whole_code_certificate=False, typed_provenance_checked=True,
        source_and_input_bytes_checked_at_exit=True, cell_transfer_recomputed=False,
        frozen_numerical_auditor_sha256=audit_provider.PREVIOUS_SHA256), 'typed independent band audit')
    if audit_provider is audited_v2:
        same(report, dict(nested_source_maps_checked=True,
            frozen_typed_v1_auditor_sha256=audited_v2.TYPED_V1_SHA256), 'typed v2 audit identity')
    summary, summary_pin = loader.dependency(report['input_summary'], audit_pin, 'audited-weight-band')
    same(summary, dict(schema='gamma-midpoint-weight-band-cover-v1', status='complete',
        complete_requested_band=True, covers_all_supports_in_band=True, covers_all_matching_profiles=True,
        every_integer_weight_prepared=True, directed=True, numerical_backend='actual-all-pair-bank',
        no_weight_interpolation=True, whole_code_certificate=False, pending_frontier=[],
        stop_reason=None, batch_source_verified_at_entry=True, batch_source_verified_at_exit=True), 'complete band')
    first, last = index.interval(summary['requested_weight_range'])
    require(last < 8192 and last-first+1 <= 64, 'Band outside producer scope')
    weights, context = list(range(first, last+1)), summary['context']
    same(summary, dict(requested_weights=weights), 'every requested integer weight')
    same(report, dict(requested_weight_range=[first, last], requested_weights=weights), 'audited weight range')
    expected = dict(midpoint.previous.EXPECTED, directed=False,
                    schema='gamma-midpoint-weight-band-context-v1', pairs=4096,
                    repetitions=256, paired_region_power=128, requested_weight_range=[first, last],
                    requested_weights=weights, source_sha256=producer.source_pins(),
                    bank_source_sha256=producer.provider.source_pins(),
                    normalized_every_original_weight=True, count_pointwise_minimum_used=False,
                    count_only_complement_above_pairs=True, no_profile_pair_table=True,
                    no_weight_interpolation=True, complete_requested_band=False)
    same(context, expected, 'gamma30/d40 band context')
    parity = {str(eps): [r for r in weights if r % 2 == eps] for eps in (0, 1)
              if any(r % 2 == eps for r in weights)}
    same(context, dict(parity_weights=parity), 'parity weight partition')
    digest = sha256(json.dumps(context, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()
    same(summary, dict(context_sha256=digest), 'canonical context')
    same(report, dict(context_sha256=digest), 'audit context binding')
    witness = context['witness']
    require(0 < fraction(witness['z']) < 1 and 0 < fraction(witness['alpha']) <= 1, 'Invalid exact marks')
    require(isinstance(witness['markers'], list) and 1 <= len(witness['markers']) <= 128
            and all(isinstance(row, list) and len(row) == 2 and all(fraction(x) > 0 for x in row)
                    for row in witness['markers']), 'Invalid exact positive markers')
    same(witness, dict(state_dimension=30, inner_width=64, gamma_extension_degree=16,
                       gamma_base_field_order=256, joint_rank_return_cap=True,
                       all_pair_before_alpha=True), 'witness construction')
    require(type(witness['precision']) is int and 128 <= witness['precision'] <= 4096
            and type(witness['grid']) is int and witness['grid'] in (64, 128, 256)
            and type(witness['packet_mantissa_bits']) is int and 32 <= witness['packet_mantissa_bits'] <= 512,
            'Invalid directed arithmetic parameters')
    roots = frozen._roots(first, last)
    descriptors = [cell.descriptor() for cell in roots]
    same(context, dict(root_cells=descriptors), 'context roots')
    same(summary, dict(root_cells=descriptors), 'summary roots')
    bits = summary['cell_failure_bits']
    require(type(bits) is int and 0 <= bits <= 4096, 'Invalid cell probability target')
    accepted, pending = producer.audit_tree(summary['nodes'], roots, bits, digest)
    require(not pending, 'Incomplete integer partition')
    same(summary, dict(accepted_leaves=accepted, tree_nodes=len(summary['nodes']),
                       evaluated_nodes=len(summary['nodes'])), 'complete tree coverage')
    nodes = {row['id']: row for row in summary['nodes']}
    for node in nodes.values():
        receipt = node['receipt']
        require(receipt['cell'] == node['cell'] and receipt['context_sha256'] == digest,
                'Cell/context receipt mismatch')
        require(type(receipt['count_upper']) is int and receipt['count_upper'] >= 0,
                'Invalid integer count cap')
        if receipt.get('zero') is True:
            same(receipt, dict(log2_upper_rational=None, log2_upper_enclosure=None), 'zero cell')
        else:
            same(receipt, dict(zero=False), 'positive cell')
            index.scalar_ceiling(receipt, 384)
    leaves = [nodes[i]['receipt'] for i in accepted]
    bounds = frozen._bounds(leaves)
    zero = not bounds
    same(summary, dict(zero=zero), 'band zero classification')
    with ctx.workprec(384):
        if zero:
            same(summary, dict(log2_upper_rational=None, log2_upper_enclosure=None), 'zero band')
            value = fresh = saved = None
        else:
            value = producer.tail.log2arb(sum((arb(2)**producer.base.rational(x) for x in bounds), arb(0)))
            fresh = producer.base.rational_ceiling(value)
            saved = index.scalar_ceiling(summary, 384)
            require(producer.base.rational(saved) >= value.upper(), 'Band ceiling misses fresh384 leaf union')
            require(saved <= -54, 'Band fails54-bit admission target')
    same(report, dict(log2_upper_rational=None if zero else str(saved),
                       saved_band_log2_upper_rational=None if zero else str(saved),
                       reaggregated_cell_log2_upper_rational=None if zero else str(fresh)),
         'audit numerical binding')
    target = summary['requested_failure_bits']
    require(type(target) is int and 0 <= target <= 4096, 'Invalid requested band target')
    same(summary, dict(meets_requested_band_target=zero or saved <= -target), 'band target flag')
    # Rehash every file the independent audit used. This is provenance, not
    # authentication or a substitute for coefficient replay.
    require(isinstance(report['checked_file_pins'], list), 'Missing audit input pins')
    checked = set()
    for link in report['checked_file_pins']:
        unused, pin = loader.dependency(link, audit_pin, 'independent-audit-input')
        checked.add((pin['path'], pin['sha256']))
    require((summary_pin['path'], summary_pin['sha256']) in checked, 'Audit did not pin its summary')
    require(isinstance(report['checked_non_json_file_pins'], list), 'Missing typed helper pins')
    for link in report['checked_non_json_file_pins']:
        require(isinstance(link, dict) and isinstance(link.get('path'), str),
                'Invalid typed helper link')
        path = Path(link['path'])
        require(path.is_absolute(), 'Typed helper path must be absolute')
        pin = source_inputs.load_source(path, link['sha256'])
        same(link, pin, 'typed helper classification')
        name = Path(os.path.relpath(path.resolve(), DIRECTORY)).as_posix()
        require(report['checked_source_sha256'].get(name) == pin['sha256'],
                'Typed helper missing from checked source closure')
    execution = report.get('execution_driver_receipt')
    if execution is not None:
        unused, pin = loader.dependency(execution, audit_pin, 'typed-execution-driver')
        require((pin['path'], pin['sha256']) in checked, 'Audit did not pin execution metadata')
    source = context.get('witness_source')
    if source is not None:
        require(isinstance(source, dict), 'Invalid witness provenance')
        path = Path(source['path'])
        unused, pin = loader.load(path if path.is_absolute() else DIRECTORY/path,
                                  expected=source['sha256'], role='band-exact-witness-source')
        require((pin['path'], pin['sha256']) in checked, 'Audit did not pin witness provenance')
    frozen._merge(pins, report['source_sha256'])
    frozen._merge(pins, report['checked_source_sha256'])
    frozen._merge(pins, context['source_sha256'])
    frozen._merge(pins, context['bank_source_sha256'])
    return dict(summary_pin, weight_range=[first, last], weights_evaluated=len(weights),
                zero=zero, log2_upper_rational=None if zero else str(saved),
                receipt_kind='typed-audited-midpoint-multiweight-cell-band', audit=audit_pin,
                context_sha256=digest)



def _band_event(path, loader, pins, source_inputs):
    report, unused = loader.load(path, role='explicit-band-audit-dispatch')
    schema = report.get('schema')
    if schema == 'gamma-midpoint-weight-band-independent-audit-v1':
        return frozen._band_event(path, loader, pins)
    if schema == 'gamma-midpoint-weight-band-typed-independent-audit-v1':
        return _typed_band_event(path, loader, pins, source_inputs, audited, schema)
    require(schema == 'gamma-midpoint-weight-band-typed-independent-audit-v2',
            'Unsupported independent audit identity')
    return _typed_band_event(path, loader, pins, source_inputs, audited_v2, schema)


def extend_v21(parent_path, audit_paths):
    """Retain all 227 v21 events and append 1..64 explicitly audited bands."""
    require(isinstance(audit_paths, (list, tuple)) and 1 <= len(audit_paths) <= 64,
            'Expected one to 64 explicit audited bands')
    loader, pins = index.Validator(F(2, 5)), source_pins()
    source_inputs = audited.Inputs()
    parent, parent_pin = loader.load(parent_path, expected=PARENT_SHA256, role='frozen-v21-index')
    same(parent, dict(schema='field128-gamma-partial-cover-index-v21',
        construction=midpoint.previous.EXPECTED['construction'], field_order=1 << 128,
        source_dimension=1 << 20, output_length=1 << 21, target_distance='2/5',
        whole_code_certificate=False, directed_aggregation=True,
        receipt_count=PARENT_EVENT_COUNT), 'frozen v21')
    require(len(parent['entries']) == PARENT_EVENT_COUNT, 'Unexpected parent event count')
    frozen._merge(pins, parent['aggregation_source_sha256'])
    entries, covered = deepcopy(parent['entries']), set()
    for entry in entries:
        path = Path(entry['path'])
        loader.load(path if path.is_absolute() else DIRECTORY/path,
                    expected=entry['sha256'], role='prior-event-bytes')
        if 'audit' in entry:
            loader.dependency(entry['audit'], parent_pin, 'inherited-prior-audit')
        lo, hi = index.interval(entry['weight_range'])
        weights = set(range(lo, hi+1))
        require(not covered & weights, 'Prior event overlap')
        covered |= weights
    same(parent, dict(covered_weight_ranges=index.ranges(covered), covered_weight_count=len(covered),
        pending_weight_ranges=index.ranges(set(range(1, 8193))-covered),
        pending_weight_count=8192-len(covered), receipt_count=len(entries)), 'prior coverage')
    require(index.union_ceiling(frozen._bounds(entries), 384)
            <= fraction(parent['covered_events_union_log2_upper_rational']), 'Prior union understated')
    added = []
    for path in audit_paths:
        entry = _band_event(path, loader, pins, source_inputs)
        lo, hi = entry['weight_range']
        weights = set(range(lo, hi+1))
        require(not covered & weights, 'New band overlaps existing or newly added coverage')
        covered |= weights
        entries.append(entry)
        added.append(entry)
    entries.sort(key=lambda row: row['weight_range'])
    union = index.union_ceiling(frozen._bounds(entries), 384)
    require(union <= -FAILURE_BITS, 'Covered-event union fails 54-bit budget')
    new_bounds = frozen._bounds(added)
    result = deepcopy(parent)
    result.update(schema='field128-gamma-partial-cover-index-v22-typed', parent_index=parent_pin, entries=entries,
        receipt_count=len(entries), covered_weight_ranges=index.ranges(covered), covered_weight_count=len(covered),
        pending_weight_ranges=index.ranges(set(range(1, 8193))-covered), pending_weight_count=8192-len(covered),
        aggregation_precision_bits=384, covered_events_union_log2_upper_rational=str(union),
        covered_events_failure_bits=FAILURE_BITS, meets_covered_events_target=True,
        aggregation_source_sha256=dict(sorted(pins.items())), added_event_count=len(added),
        added_weight_count=sum(row['weights_evaluated'] for row in added),
        added_events_zero=not new_bounds,
        added_events_union_log2_upper_rational=str(index.union_ceiling(new_bounds, 384)) if new_bounds else None,
        latest_audited_multiweight_events=added,
        prior_semantic_audits_inherited=True, new_count_audits_inherited=True,
        new_integer_partitions_rechecked=True, new_cell_unions_reaggregated=True,
        coefficient_replay=False, normalized_regional_replay=False,
        input_file_pins_verified=len(loader.files),
        non_json_input_file_pins_verified=len(source_inputs.non_json_pins),
        checked_non_json_file_pins=list(source_inputs.non_json_pins.values()),
        typed_and_legacy_audit_identities_explicit=True, helper_source_files_verified=len(pins),
        all_receipt_bytes_match=True, all_helper_source_closures_match=True, source_pin_conflicts=False,
        no_weight_interpolation=True,
        zero_event_semantics='zero=true and null log ceiling; omitted from positive union only',
        aggregation_method=f'Directed384-bit union over {len(entries)} nonoverlapping weight ranges; each complete new band counted once',
        scope=f'Only {len(covered)} listed weights; {8192-len(covered)} remain pending; no whole-code certificate')
    counts = result['producer_receipt_counts']
    name = 'hd_field128_gamma_midpoint_weight_bands'
    counts[name] = counts.get(name, 0)+len(added)
    band._check_links(tuple(loader.files.values()))
    source_inputs.check()
    frozen.producer.check_source_pins(pins)
    return result

