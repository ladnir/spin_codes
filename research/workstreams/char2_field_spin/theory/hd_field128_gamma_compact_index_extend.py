"""Append audited compact bands to a gamma30/d40 partial index, without replay.

The frozen schema resolver validates each new event. This small version-
independent accumulator inherits prior semantic audits while rechecking
all prior event bytes, source pins, coverage, and the directed384-bit union.
"""

from copy import deepcopy
from fractions import Fraction as F
from pathlib import Path
import re

from flint import arb, ctx

import hd_field128_gamma_near_face_index_update as frozen


index, band = frozen.index, frozen.band
require, same, fraction = frozen.require, frozen.same, frozen.fraction
DIRECTORY = Path(__file__).resolve().parent


def source_pins():
    return index.source_closure(Path(__file__).name)


def _merge(pins, extra):
    for name, digest in extra.items():
        require(name not in pins or pins[name] == digest, 'Conflicting source pin: '+name)
        pins[name] = digest


def _inline_context(context, proposal_path, loader):
    frozen.audited._verify_witness(context)
    proposal, pin = loader.load(proposal_path, role='immutable-inline-proposal')
    require(context['witness_source'] == proposal, 'Inline proposal differs from immutable proposal bytes')
    same(proposal, dict(schema='field128-gamma-exact-proposal-v1', construction=frozen.audited.EXPECTED['construction'],
                        state_dimension=30, target_distance='2/5', directed=False, whole_code_certificate=False), 'inline proposal')
    require(type(proposal['calibration_active_groups']) is int and 1 <= proposal['calibration_active_groups'] <= 8192,
            'Invalid proposal calibration weight')
    same(context['witness'], proposal['witness'], 'proposal exact bank literals')
    require(proposal['marker_count'] == len(context['witness']['markers']), 'Proposal marker count differs')
    source = proposal['source']
    path = Path(source['path'])
    loader.load(path if path.is_absolute() else DIRECTORY/path, expected=source['sha256'], role='proposal-source-literals')
    return pin


def audit_inline_band(path, proposal_path):
    """Audit original rows; resolve an inline proposal against its saved bytes."""
    own_pins, loader = source_pins(), index.Validator(F(2, 5))
    summary, summary_pin = loader.load(path)
    same(summary, dict(schema='gamma-near-face-band-summary-v1', construction=frozen.audited.EXPECTED['construction'],
                       complete_requested_range=True, covers_all_supports_in_range=True,
                       every_integer_weight_evaluated=True, covers_all_even_odd_splits=True,
                       pending_weight_ranges=[], partial_weight_receipt=None, failure=None, directed=True,
                       whole_code_certificate=False, regional_caches_built=1, common_actual_bank=True,
                       no_weight_interpolation=True, source_pins_unchanged_at_exit=True), 'complete inline band')
    context, context_pin = loader.dependency(summary['context'], summary_pin, 'shared-band-context')
    proposal_pin = _inline_context(context, proposal_path, loader)
    lo, hi = index.interval(summary['requested_weight_range'])
    same(context, dict(requested_weight_range=[lo, hi]), 'context range')
    same(summary, dict(weights_evaluated=hi-lo+1, completed_weight_range=[lo, hi]), 'complete range')
    links, reports, bounds = summary['weight_receipts'], [], []
    require(isinstance(links, list) and len(links) == hi-lo+1, 'Incomplete weight links')
    with ctx.workprec(384):
        for weight, link in zip(range(lo, hi+1), links):
            row, row_pin = loader.dependency(link, summary_pin, 'complete-weight-rows')
            require(row['context'] == summary['context'], 'Weight context mismatch')
            same(row, dict(active_groups=weight, status='complete', complete_selected_weight=True,
                           covers_all_supports_at_weight=True, covers_all_even_odd_splits=True,
                           common_potential_used=False, directed=True, whole_code_certificate=False, failure=None,
                           requires_linked_context=True, source_pins_checked_before_publication=True,
                           count_provenance=context['count_provenance'],
                           requested_failure_bits=context['requested_failure_bits']), 'weight scope')
            same(link, {key: row[key] for key in ('active_groups', 'status', 'complete_selected_weight',
                        'log2_upper_rational', 'completed_even_rows', 'pending_even_full_range')}, 'weight link')
            report = frozen.audited.audit_rows(row, context, precision=384)
            bound = fraction(row['log2_upper_rational'])
            same(row, dict(meets_requested_selected_weight_target=bound <= -context['requested_failure_bits']), 'weight target')
            reports.append(dict(report, receipt=row_pin))
            bounds.append(bound)
        union, saved = index.union_ceiling(bounds, 384), fraction(summary['log2_upper_rational'])
        require(saved >= union, 'Saved band ceiling misses directed384 union')
        same(summary, dict(completed_union_log2_upper_rational=str(saved), requested_failure_bits=context['requested_failure_bits'],
                           requested_range_meets_target=saved <= -context['requested_failure_bits'],
                           weights_above_target=[lo+i for i, value in enumerate(bounds) if value > -context['requested_failure_bits']]),
             'band union/target')
    require(source_pins() == own_pins, 'Audit source changed')
    band.selected.check_source_pins(context['source_sha256'])
    band._check_links(tuple(loader.files.values()))
    return dict(schema='gamma-inline-proposal-compact-band-audit-v1', audit_result='PASS',
                construction=frozen.audited.EXPECTED['construction'], state_dimension=30, target_distance='2/5',
                weight_range=[lo, hi], summary=summary_pin, context=context_pin, resolved_inline_proposal=proposal_pin,
                weight_audits=reports, count_cutoff=context['count_provenance']['maximum_refined_singletons'],
                exact_marker_count=len(context['witness']['markers']), precision_bits=384,
                saved_band_log2_upper_rational=str(saved), directed_union_log2_upper_rational=str(union),
                complete_range_event=True, coefficient_replay=False, whole_code_certificate=False,
                checked_file_pins=list(loader.files.values()), source_sha256=own_pins,
                checked_producer_source_sha256=context['source_sha256'],
                audit_limit='Frozen independent exact-count/coverage/row-prefix audit; inline proposal equals saved file; '
                            'no coefficient/FFT replay or authentication of original execution')


def _event(summary_path, audit_path, loader):
    report, report_pin = loader.load(audit_path, role='compact-band-audit')
    if report.get('schema') != 'gamma-inline-proposal-compact-band-audit-v1':
        return frozen._event('compact-band', summary_path, audit_path, loader)
    same(report, dict(audit_result='PASS', construction=frozen.audited.EXPECTED['construction'], state_dimension=30,
                      target_distance='2/5', complete_range_event=True, coefficient_replay=False,
                      precision_bits=384, source_sha256=source_pins()), 'inline audit')
    band._check_links(report['checked_file_pins'])
    summary, pin = loader.load(summary_path, expected=report['summary']['sha256'])
    same(summary, dict(schema='gamma-near-face-band-summary-v1', complete_requested_range=True,
                       covers_all_supports_in_range=True, every_integer_weight_evaluated=True,
                       requested_weight_range=report['weight_range'],
                       log2_upper_rational=report['saved_band_log2_upper_rational']), 'inline audit binding')
    same(report['context'], summary['context'], 'context path/hash identity')
    context, unused = loader.dependency(summary['context'], pin, 'audited-band-context')
    proposal_pin = _inline_context(context, report['resolved_inline_proposal']['path'], loader)
    require(proposal_pin['sha256'] == report['resolved_inline_proposal']['sha256'], 'Resolved proposal hash mismatch')
    bound = fraction(summary['log2_upper_rational'])
    require(bound <= -54, 'New event fails54-bit admission target')
    return dict(path=pin['path'], sha256=pin['sha256'], weight_range=report['weight_range'],
                log2_upper_rational=str(bound), receipt_kind='compact-band', audit=report_pin,
                context=summary['context']), report['source_sha256']


def extend_compact_index(parent_path, events):
    """events is a bounded sequence of (complete summary path, audit path)."""
    require(isinstance(events, (list, tuple)) and 1 <= len(events) <= 32,
            'Expected one to32 explicitly audited compact bands')
    loader = index.Validator(F(2, 5))
    parent, parent_pin = loader.load(parent_path)
    match = re.fullmatch(r'field128-gamma-partial-cover-index-v([1-9][0-9]*)', parent.get('schema', ''))
    require(match is not None and 12 <= int(match[1]) <= 999, 'Expected a gamma index version in12..999')
    version = int(match[1])+1
    same(parent, dict(construction=frozen.audited.EXPECTED['construction'], field_order=1 << 128,
                      source_dimension=1 << 20, output_length=1 << 21, target_distance='2/5',
                      whole_code_certificate=False, directed_aggregation=True), 'prior gamma index')
    pins, entries, covered = dict(parent['aggregation_source_sha256']), deepcopy(parent['entries']), set()
    _merge(pins, source_pins())
    for entry in entries:
        path = Path(entry['path'])
        loader.load(path if path.is_absolute() else DIRECTORY/path, expected=entry['sha256'], role='prior-event-bytes')
        lo, hi = index.interval(entry['weight_range'])
        weights = set(range(lo, hi+1))
        require(not covered & weights, 'Prior events overlap')
        covered |= weights
    require(index.ranges(covered) == parent['covered_weight_ranges'] and len(covered) == parent['covered_weight_count'],
            'Prior covered ranges/count disagree')
    pending = set(range(1, 8193))-covered
    require(index.ranges(pending) == parent['pending_weight_ranges'] and len(pending) == parent['pending_weight_count'],
            'Prior pending ranges/count disagree')
    require(len(entries) == parent['receipt_count'], 'Prior event count disagrees')
    old_union = index.union_ceiling([fraction(row['log2_upper_rational']) for row in entries], 384)
    require(old_union <= fraction(parent['covered_events_union_log2_upper_rational']), 'Prior union understated')
    added, new_weights = [], 0
    for summary_path, audit_path in events:
        entry, closure = _event(summary_path, audit_path, loader)
        _merge(pins, closure)
        lo, hi = index.interval(entry['weight_range'])
        weights = set(range(lo, hi+1))
        require(not covered & weights, 'New compact band overlaps covered weights')
        covered |= weights
        entries.append(entry)
        added.append(entry)
        new_weights += len(weights)
    entries.sort(key=lambda row: row['weight_range'])
    with ctx.workprec(384):
        value = frozen.tail.log2arb(sum((arb(2)**frozen.base.rational(fraction(row['log2_upper_rational']))
                                        for row in entries), arb(0)))
        union, enclosure = str(frozen.base.rational_ceiling(value)), str(value)
    result = deepcopy(parent)
    result.update(schema=f'field128-gamma-partial-cover-index-v{version}', parent_index=parent_pin,
                  entries=entries, receipt_count=len(entries), covered_weight_ranges=index.ranges(covered),
                  pending_weight_ranges=index.ranges(set(range(1, 8193))-covered), covered_weight_count=len(covered),
                  pending_weight_count=8192-len(covered), aggregation_precision_bits=384,
                  aggregation_source_sha256=dict(sorted(pins.items())),
                  covered_events_union_log2_upper_rational=union, covered_events_union_log2_upper_enclosure=enclosure,
                  added_event_count=len(added), added_weight_count=new_weights,
                  added_events_union_log2_upper_rational=str(index.union_ceiling([fraction(row['log2_upper_rational']) for row in added], 384)),
                  latest_explicit_near_face_events=added, prior_semantic_audits_inherited=True,
                  all_receipt_bytes_match=True, all_helper_source_closures_match=True, source_pin_conflicts=False,
                  input_file_pins_verified=len(loader.files), helper_source_files_verified=len(pins),
                  compact_context_resolution_explicit=True, no_weight_interpolation=True,
                  aggregation_method=f'Directed384-bit union of{len(entries)} disjoint events; terminal once; no interpolation',
                  scope=f'Only{len(covered)} listed weights;{8192-len(covered)} remain pending; no whole-code certificate')
    counts = result['producer_receipt_counts']
    counts['hd_field128_gamma_near_face_band'] = counts.get('hd_field128_gamma_near_face_band', 0)+len(added)
    band._check_links(tuple(loader.files.values()))
    band.selected.check_source_pins(pins)
    return result
