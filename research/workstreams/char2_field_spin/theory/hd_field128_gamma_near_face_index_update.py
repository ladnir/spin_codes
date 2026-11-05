"""Explicit saved-row slices and disjoint gamma30/d40 index updates.

No coefficient arithmetic is performed. Complete compact-band audits are
resolved by their actual schema; a slice retains its original context and
weight files. The selected768 event requires its independent saved-row audit.
"""

from copy import deepcopy
from fractions import Fraction as F
from pathlib import Path

from flint import arb, ctx

import hd_field128_gamma_near_face_band_audit as audited


band, index, base, tail = audited.band, audited.index, audited.base, audited.tail
require, same, fraction = audited.require, audited.same, audited.fraction
DIRECTORY = Path(__file__).resolve().parent


def source_pins():
    return index.source_closure(Path(__file__).name)


def checked_band_audit(path, loader):
    audit, audit_pin = loader.load(path, role='explicit-compact-band-audit')
    same(audit, dict(schema='gamma-near-face-compact-band-audit-v1', audit_result='PASS',
                    construction=audited.EXPECTED['construction'], state_dimension=30,
                    target_distance='2/5', complete_range_event=True, coefficient_replay=False,
                    source_sha256=audited.source_pins()), 'compact audit')
    require(audit['precision_bits'] >= 384, 'Insufficient audit precision')
    band._check_links(audit['checked_file_pins'])
    summary, summary_pin = loader.dependency(audit['summary'], audit_pin, 'audited-band-summary')
    same(summary, dict(schema='gamma-near-face-band-summary-v1', complete_requested_range=True,
                       covers_all_supports_in_range=True, every_integer_weight_evaluated=True,
                       requested_weight_range=audit['weight_range'],
                       log2_upper_rational=audit['saved_band_log2_upper_rational']), 'audit binding')
    same(audit['context'], summary['context'], 'context path/hash identity')
    context, unused = loader.dependency(summary['context'], summary_pin, 'audited-band-context')
    audited._verify_witness(context)
    return audit, audit_pin, summary, summary_pin


def _slice(parent, parent_pin, audit, audit_pin, first, last):
    lo, hi = index.interval(parent['requested_weight_range'])
    require(type(first) is int and type(last) is int and lo <= first <= last <= hi, 'Invalid saved-row slice')
    links = parent['weight_receipts'][first-lo:last-lo+1]
    require([row['active_groups'] for row in links] == list(range(first, last+1)), 'Slice weights incomplete')
    require(all(row['complete_selected_weight'] and row['status'] == 'complete' for row in links), 'Partial slice weight')
    reports = audit['weight_audits'][first-lo:last-lo+1]
    require(len(reports) == len(links), 'Missing selected-weight audits')
    for report, link in zip(reports, links):
        same(report, dict(active_groups=link['active_groups'], saved_log2_upper_rational=link['log2_upper_rational']),
             'slice row audit')
        require(report['receipt']['sha256'] == link['sha256'], 'Slice audited row hash mismatch')
    bound = index.union_ceiling([fraction(row['log2_upper_rational']) for row in links], 384)
    target = parent['requested_failure_bits']
    return dict(schema='gamma-near-face-saved-weight-slice-v1', construction=audited.EXPECTED['construction'],
                state_dimension=30, target_distance='2/5', parent_band=parent_pin, parent_band_audit=audit_pin,
                context=parent['context'], weight_range=[first, last], weight_receipts=deepcopy(links),
                weights_evaluated=len(links), complete_requested_range=True, covers_all_supports_in_range=True,
                every_integer_weight_evaluated=True, all_rows_inherited_from_verified_complete_weights=True,
                log2_upper_rational=str(bound), aggregation_precision_bits=384, requested_failure_bits=target,
                requested_range_meets_target=bound <= -target,
                weights_above_target=[row['active_groups'] for row in links if fraction(row['log2_upper_rational']) > -target],
                excluded_parent_weights=index.ranges(set(range(lo, hi+1))-set(range(first, last+1))),
                parent_event_included=False, no_weight_interpolation=True, directed=True,
                coefficient_replay=False, whole_code_certificate=False, source_sha256=source_pins(),
                scope='Only the explicit complete saved weights in this slice; parent event and all other weights excluded')


def slice_receipt(parent_audit_path, first, last):
    loader = index.Validator(F(2, 5))
    audit, audit_pin, parent, parent_pin = checked_band_audit(parent_audit_path, loader)
    result = _slice(parent, parent_pin, audit, audit_pin, first, last)
    band._check_links(tuple(loader.files.values()))
    band.selected.check_source_pins(result['source_sha256'])
    return result


def _event(kind, path, audit_path, loader):
    if kind == 'compact-band':
        report, report_pin, record, pin = checked_band_audit(audit_path, loader)
        provided, unused = loader.load(path, expected=pin['sha256'])
        require(provided == record, 'Wrong compact summary')
        bounds = record['requested_weight_range']
        pins = report['source_sha256']
    elif kind == 'saved-slice':
        record, pin = loader.load(path)
        same(record, dict(schema='gamma-near-face-saved-weight-slice-v1', source_sha256=source_pins()), 'slice producer')
        report, report_pin, parent, parent_pin = checked_band_audit(record['parent_band_audit']['path'], loader)
        require(report_pin['sha256'] == record['parent_band_audit']['sha256'], 'Wrong parent audit hash')
        require(record == _slice(parent, parent_pin, report, report_pin, *record['weight_range']), 'Derived slice differs from verified rows')
        bounds, pins = record['weight_range'], record['source_sha256']
    elif kind == 'selected768':
        record, pin = loader.load(path)
        report, report_pin = loader.load(audit_path)
        same(report, dict(schema='field128-gamma-r768-directed-receipt-audit-v1', audit_pass=True,
                          construction=audited.EXPECTED['construction'], state_dimension=30,
                          target_distance='2/5', active_groups=768, compact_rows_checked=385,
                          profile_pairs_checked=148225, metadata_lineage_and_source_pins_checked=True), 'selected768 audit')
        same(record, dict(audited.EXPECTED, active_groups=768, complete_selected_weight=True,
                          covers_all_supports_at_weight=True, covers_all_even_odd_splits=True), 'selected768 receipt')
        require(report['input_sha256']['summary.json'] == pin['sha256'], 'Selected768 audit/summary hash mismatch')
        require(report['numeric']['stored_final_ceiling'] == record['log2_upper_rational'], 'Selected768 ceiling mismatch')
        require(report['numeric']['audit_precision'] >= 384 and report['numeric']['stored_ceiling_dominates_reaggregated_enclosure'],
                'Incomplete selected768 numerical audit')
        for name, digest in report['input_sha256'].items():
            target = Path(name)
            loader.load(target if target.is_absolute() else Path(report['receipt_directory'])/target,
                        expected=digest, role='selected768-audited-input')
        same(record, dict(source_sha256=band.selected.source_pins(), bank_source_sha256=band.provider.source_pins()),
             'selected768 producer pins')
        bounds, pins = [768, 768], record['source_sha256']
    else:
        raise ValueError('Unsupported explicit event schema')
    bound = fraction(record['log2_upper_rational'])
    require(bound <= -54, 'New event fails54-bit admission target')
    return dict(path=pin['path'], sha256=pin['sha256'], weight_range=bounds,
                log2_upper_rational=str(bound), receipt_kind=kind, audit=report_pin,
                **({'context': record['context']} if 'context' in record else {})), pins


def extend_index(parent_path, events, *, version):
    """Preserve prior semantic audits; recheck their bytes, pins, ranges, and union."""
    require(version in (10, 11, 12), 'This adapter only creates v10, v11, or v12')
    loader = index.Validator(F(2, 5))
    parent, parent_pin = loader.load(parent_path)
    same(parent, dict(schema=f'field128-gamma-partial-cover-index-v{version-1}',
                      construction=audited.EXPECTED['construction'], field_order=1 << 128,
                      source_dimension=1 << 20, output_length=1 << 21, target_distance='2/5',
                      whole_code_certificate=False), 'prior gamma index')
    require(events and len(events) <= 2, 'Expected one or two explicitly audited events')
    entries, covered = deepcopy(parent['entries']), set()
    pins = dict(parent['aggregation_source_sha256'])
    for name, digest in source_pins().items():
        require(name not in pins or pins[name] == digest, 'Conflicting helper source: '+name)
        pins[name] = digest
    for entry in entries:
        target = Path(entry['path'])
        loader.load(target if target.is_absolute() else DIRECTORY/target, expected=entry['sha256'], role='prior-event-bytes')
        lo, hi = index.interval(entry['weight_range'])
        weights = set(range(lo, hi+1))
        require(not covered & weights, 'Prior overlap')
        covered |= weights
    require(index.ranges(covered) == parent['covered_weight_ranges'] and len(covered) == parent['covered_weight_count'],
            'Prior covered ranges/count differ')
    require(index.ranges(set(range(1, 8193))-covered) == parent['pending_weight_ranges'], 'Prior pending ranges differ')
    require(len(entries) == parent['receipt_count'], 'Prior event count differs')
    old_union = index.union_ceiling([fraction(row['log2_upper_rational']) for row in entries], 384)
    require(old_union <= fraction(parent['covered_events_union_log2_upper_rational']), 'Prior union understated')
    added, new_weight_count = [], 0
    for kind, path, audit_path in events:
        entry, closure = _event(kind, path, audit_path, loader)
        for name, digest in closure.items():
            require(name not in pins or pins[name] == digest, 'Conflicting event source: '+name)
            pins[name] = digest
        lo, hi = index.interval(entry['weight_range'])
        weights = set(range(lo, hi+1))
        require(not covered & weights, 'New event overlaps covered weights')
        covered |= weights
        entries.append(entry)
        added.append(entry)
        new_weight_count += len(weights)
    entries.sort(key=lambda row: row['weight_range'])
    with ctx.workprec(384):
        value = tail.log2arb(sum((arb(2)**base.rational(fraction(row['log2_upper_rational'])) for row in entries), arb(0)))
        union, enclosure = str(base.rational_ceiling(value)), str(value)
    result = deepcopy(parent)
    result.update(schema=f'field128-gamma-partial-cover-index-v{version}', parent_index=parent_pin,
                  entries=entries, receipt_count=len(entries), covered_weight_ranges=index.ranges(covered),
                  pending_weight_ranges=index.ranges(set(range(1, 8193))-covered), covered_weight_count=len(covered),
                  pending_weight_count=8192-len(covered), aggregation_precision_bits=384,
                  aggregation_source_sha256=dict(sorted(pins.items())),
                  covered_events_union_log2_upper_rational=union, covered_events_union_log2_upper_enclosure=enclosure,
                  added_event_count=len(added), added_weight_count=new_weight_count,
                  added_events_union_log2_upper_rational=str(index.union_ceiling([fraction(row['log2_upper_rational']) for row in added], 384)),
                  latest_explicit_near_face_events=added, prior_semantic_audits_inherited=True,
                  all_receipt_bytes_match=True, all_helper_source_closures_match=True, source_pin_conflicts=False,
                  input_file_pins_verified=len(loader.files), helper_source_files_verified=len(pins),
                  failed_752_767_parent_band_included=False, compact_context_resolution_explicit=True,
                  aggregation_method=f'Directed384-bit union of{len(entries)} disjoint events; no interpolation; terminal once',
                  scope=f'Only{len(covered)} listed weights;{8192-len(covered)} weights remain pending; no whole-code certificate')
    for entry in added:
        name = {'selected768': 'hd_field128_gamma_near_face_certificate', 'compact-band': 'hd_field128_gamma_near_face_band',
                'saved-slice': 'hd_field128_gamma_near_face_index_update:saved-slice'}[entry['receipt_kind']]
        result['producer_receipt_counts'][name] = result['producer_receipt_counts'].get(name, 0)+1
    band._check_links(tuple(loader.files.values()))
    band.selected.check_source_pins(pins)
    return result
