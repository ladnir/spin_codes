"""Explicit v14 addition of five reviewed selected-weight events, not bands."""

from copy import deepcopy
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

from flint import arb, ctx

import hd_field128_gamma_midpoint_receipt_audit as midpoint
import hd_field128_gamma_compact_index_extend as prior


index, band, base, tail = prior.index, prior.band, midpoint.base, midpoint.tail
require, same, fraction = prior.require, prior.same, prior.fraction
DIRECTORY = Path(__file__).resolve().parent


def source_pins():
    return dict(sorted({**index.source_closure(Path(__file__).name), **midpoint.source_pins()}.items()))


def complete_selected(record, weight, pins):
    n = weight//2+1
    same(record, dict(midpoint.previous.EXPECTED, active_groups=weight, status='complete', failure=None,
                      complete_selected_weight=True, covers_all_supports_at_weight=True,
                      covers_all_even_odd_splits=True, repetitions=256, paired_region_power=128,
                      common_potential_used=False, source_pins_unchanged_at_exit=True,
                      completed_even_rows=n, expected_even_rows=n, completed_even_full_range=[0, n-1],
                      pending_even_full_range=None, profile_pairs_evaluated=n*n, expected_profile_pairs=n*n,
                      source_sha256=pins, coverage_driver_sha256=pins,
                      bank_source_sha256=band.provider.source_pins()), 'complete selected event')
    require(len(record['rows']) == n, 'Incomplete inline rows')
    with ctx.workprec(384):
        bound = index.scalar_ceiling(record, 384)
        require(bound <= -54, 'Selected event misses54-bit admission target')
        prefix = midpoint.previous.outer_prefix(record, weight)
        total = arb(0)
        for a, row in enumerate(record['rows']):
            same(row, dict(even_full=a, odd_full_range=[0, n-1], profile_pairs_evaluated=n,
                           complete_even_profile_row=True, covers_all_odd_profiles_for_this_even_profile=True,
                           complete_selected_weight=False, covers_all_supports_at_weight=False,
                           common_potential_used=False), 'row coverage')
            value = midpoint.previous._log_record(row['row_support'], 384)
            if value is not None:
                total += arb(2)**value.upper()
        require(not total.is_zero() and base.rational(bound) >= (prefix+tail.log2arb(total)).upper(),
                'Selected ceiling misses fresh prefix and saved rows')
    return bound


def near_event(audit_path, expected_weight, loader):
    report, audit_pin = loader.load(audit_path, role='inherited-near-face-audit')
    r, n = expected_weight, expected_weight//2+1
    same(report, dict(schema=f'field128-gamma-r{r}-directed-receipt-audit-v1', audit_pass=True,
                      construction=midpoint.previous.EXPECTED['construction'], state_dimension=30,
                      field_order=str(1 << 128), target_distance='2/5', active_groups=r,
                      compact_rows_checked=n, profile_pairs_checked=n*n,
                      metadata_lineage_and_source_pins_checked=True, coefficient_replay=False), 'near-face audit')
    folder = Path(report['receipt_directory'])
    objects = {}
    for name, digest in report['input_sha256'].items():
        path = Path(name)
        objects[name], unused = loader.load(path if path.is_absolute() else folder/path, expected=digest, role='audited-near-face-input')
    record, pin = loader.load(folder/'summary.json', expected=report['input_sha256']['summary.json'])
    bound = complete_selected(record, r, band.selected.source_pins())
    require(record['count_provenance'] == band.selected.count_provenance(64), 'Wrong near-face count policy')
    same(report['numeric'], dict(audit_precision=384, stored_final_ceiling=str(bound),
                                 stored_ceiling_dominates_reaggregated_enclosure=True), 'near-face arithmetic audit')
    context_record = objects['context.json']
    context = context_record['context']
    context_hash = sha256(json.dumps(context, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()
    require(context_hash == context_record['context_sha256'] == record['context_sha256'] == report['context_sha256'],
            'Context digest mismatch')
    same(record, {k: v for k, v in context.items() if k != 'scope'}, 'context metadata')
    same(record['witness'], objects['proposal.json']['witness'], 'proposal literals')
    require(record['witness'] == objects['bank.json']['witness'], 'Bank witness differs')
    for a, row in enumerate(record['rows']):
        key = str(Path('rows')/f'even_{a:04d}.json')
        external = objects[key]
        same(external, row, 'external row content')
        same(external, dict(active_groups=r, context_sha256=context_hash), 'external row context')
    return dict(pin, weight_range=[r, r], log2_upper_rational=str(bound),
                receipt_kind='audited-near-face-selected', audit=audit_pin), report['source_sha256']


def midpoint_events(audit_path, loader):
    report, audit_pin = loader.load(audit_path, role='inherited-midpoint-audit')
    same(report, dict(schema='gamma-midpoint-selected-receipts-audit-v1', audit_result='PASS',
                      construction=midpoint.previous.EXPECTED['construction'], state_dimension=30,
                      target_distance='2/5', precision_bits=384, source_sha256=midpoint.source_pins(),
                      coefficient_replay=False, whole_code_certificate=False), 'midpoint audit')
    band._check_links(report['checked_file_pins'])
    entries = []
    for weight in (511, 758):
        matches = [row for row in report['weight_audits'] if row['active_groups'] == weight]
        require(len(matches) == 1, 'Missing or duplicated midpoint audit weight')
        verified = matches[0]
        record, pin = loader.dependency(verified['receipt'], audit_pin, 'audited-midpoint-receipt')
        bound = complete_selected(record, weight, midpoint.producer.source_pins())
        same(verified, dict(saved_log2_upper_rational=str(bound), all_producer_midpoint_and_final_counts_match=True,
                            saved_ceiling_dominates_row_enclosures=True), 'midpoint arithmetic/count audit')
        proposal, unused = loader.dependency(verified['proposal'], audit_pin, 'midpoint-proposal')
        same(record['witness'], proposal['witness'], 'midpoint exact proposal')
        entries.append(dict(pin, weight_range=[weight, weight], log2_upper_rational=str(bound),
                            receipt_kind='audited-midpoint-selected', audit=audit_pin))
    return entries, report['source_sha256']


def extend_v13(parent_path, near_audits, midpoint_audit):
    require(set(near_audits) == {128, 384, 1536}, 'Exactly the three requested near-face audits are required')
    loader = index.Validator(F(2, 5))
    parent, parent_pin = loader.load(parent_path)
    same(parent, dict(schema='field128-gamma-partial-cover-index-v13', construction=midpoint.previous.EXPECTED['construction'],
                      field_order=1 << 128, source_dimension=1 << 20, output_length=1 << 21,
                      target_distance='2/5', whole_code_certificate=False, directed_aggregation=True), 'frozen v13')
    entries, covered, pins = deepcopy(parent['entries']), set(), dict(parent['aggregation_source_sha256'])
    prior._merge(pins, source_pins())
    for entry in entries:
        path = Path(entry['path'])
        loader.load(path if path.is_absolute() else DIRECTORY/path, expected=entry['sha256'], role='prior-event-bytes')
        lo, hi = index.interval(entry['weight_range'])
        weights = set(range(lo, hi+1))
        require(not covered & weights, 'Prior overlap')
        covered |= weights
    require(index.ranges(covered) == parent['covered_weight_ranges'] and len(covered) == parent['covered_weight_count'], 'Prior coverage differs')
    require(index.ranges(set(range(1, 8193))-covered) == parent['pending_weight_ranges'], 'Prior pending ranges differ')
    require(len(entries) == parent['receipt_count'], 'Prior event count differs')
    require(index.union_ceiling([fraction(x['log2_upper_rational']) for x in entries], 384)
            <= fraction(parent['covered_events_union_log2_upper_rational']), 'Prior union understated')
    added = []
    for weight, path in sorted(near_audits.items()):
        entry, closure = near_event(path, weight, loader)
        prior._merge(pins, closure)
        added.append(entry)
    middle, closure = midpoint_events(midpoint_audit, loader)
    prior._merge(pins, closure)
    added.extend(middle)
    for entry in added:
        weight = entry['weight_range'][0]
        require(entry['weight_range'] == [weight, weight] and weight not in covered, 'New selected event overlaps')
        covered.add(weight)
        entries.append(entry)
    entries.sort(key=lambda row: row['weight_range'])
    with ctx.workprec(384):
        value = tail.log2arb(sum((arb(2)**base.rational(fraction(row['log2_upper_rational'])) for row in entries), arb(0)))
        union, enclosure = str(base.rational_ceiling(value)), str(value)
    result = deepcopy(parent)
    result.update(schema='field128-gamma-partial-cover-index-v14', parent_index=parent_pin, entries=entries,
                  receipt_count=len(entries), covered_weight_ranges=index.ranges(covered),
                  pending_weight_ranges=index.ranges(set(range(1, 8193))-covered), covered_weight_count=len(covered),
                  pending_weight_count=8192-len(covered), aggregation_precision_bits=384,
                  covered_events_union_log2_upper_rational=union, covered_events_union_log2_upper_enclosure=enclosure,
                  aggregation_source_sha256=dict(sorted(pins.items())), added_event_count=5, added_weight_count=5,
                  added_events_union_log2_upper_rational=str(index.union_ceiling([fraction(x['log2_upper_rational']) for x in added], 384)),
                  latest_explicit_selected_events=added, prior_semantic_audits_inherited=True,
                  input_file_pins_verified=len(loader.files), helper_source_files_verified=len(pins),
                  all_receipt_bytes_match=True, all_helper_source_closures_match=True, source_pin_conflicts=False,
                  excluded_midpoint_weights=[756, 757, 768], no_weight_interpolation=True,
                  aggregation_method=f'Directed384-bit union of{len(entries)} disjoint events; five selected weights, no fake bands',
                  scope=f'Only{len(covered)} listed weights;{8192-len(covered)} remain pending; no whole-code certificate')
    counts = result['producer_receipt_counts']
    counts['hd_field128_gamma_near_face_certificate'] = counts.get('hd_field128_gamma_near_face_certificate', 0)+3
    counts['hd_field128_gamma_midpoint_certificate'] = counts.get('hd_field128_gamma_midpoint_certificate', 0)+2
    band._check_links(tuple(loader.files.values()))
    band.selected.check_source_pins(pins)
    return result
