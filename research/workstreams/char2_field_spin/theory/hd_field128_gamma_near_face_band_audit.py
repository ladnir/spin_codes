"""Resolve compact gamma30/d40 bands and audit saved rows without coefficient replay.

This consumer explicitly resolves context and weight byte hashes. It checks
exact profile counts, complete coverage, fresh outer prefixes, and directed
384-bit reaggregation. Source/byte pins are provenance, not authentication
of the original coefficient execution. Partial bands cannot become events.
"""

from copy import deepcopy
from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb
from pathlib import Path

from flint import arb, ctx

import hd_field128_gamma_near_face_band as band
import hd_field128_lift_partial_index as index


DIRECTORY = Path(__file__).resolve().parent
base, tail = band.base, band.tail
require, same, fraction = index.require, index.same, index.fraction
EXPECTED = dict(construction=band.provider.local.construction(30), state_dimension=30,
                field_order=1 << 128, field_bits=128, source_dimension=1 << 20,
                output_length=1 << 21, outer_groups=8192, target_distance='2/5',
                forbidden_weight=838860, outer_measure='invertible-toeplitz',
                projective_divisor=str((1 << 128)-1), inner_width=64, gamma_extension_degree=16,
                gamma_base_field_order=256, initial_state='zero', final_state='unrestricted',
                pair_swaps=True, state_persists_across_regions=True, full_field_random_coefficients=True,
                all_pair_before_alpha=True, counts_outside_alpha=True, no_coefficient_truncation=True,
                directed=True, whole_code_certificate=False)


def source_pins():
    return index.source_closure(Path(__file__).name)


def _runs(length, weight):
    """Independent binomial formula; no frozen run recurrence is called."""
    result = [0]*(weight+1)
    if weight in (0, length):
        result[0] = 1
        return result
    left, right = 1, 1
    for k in range(1, min(weight, length-weight)+1):
        numerator = length*left*right
        result[k], remainder = divmod(numerator, k)
        require(not remainder, 'Nonintegral cyclic run formula')
        left = left*(weight-k)//k
        right = right*(length-weight-k)//k
    return result


class IndependentCounts:
    def __init__(self, pairs, weight, cutoff):
        self.r, self.full = weight, range(max(0, weight-pairs), weight//2+1)
        self.runs = _runs(2*pairs, weight)
        self.contracted = _runs(pairs, weight//2) if weight % 2 == 0 else None
        self.marginals, self.near = {}, {}
        for f in self.full:
            singles = weight-2*f
            self.marginals[f] = comb(pairs, f)*comb(pairs-f, singles)*(1 << singles)
            if singles <= cutoff:
                prefix = [0]
                for value in _runs(pairs-singles, f):
                    prefix.append(prefix[-1]+value)
                self.near[f] = singles, comb(pairs, singles)*(1 << singles), prefix

    def count(self, a, b):
        k = self.r-a-b
        if 2*a == self.r or 2*b == self.r:
            result = self.contracted[k] if k < len(self.contracted) else 0
        else:
            result = min(self.runs[k], self.marginals[a], self.marginals[b])
        for f in (a, b):
            if f in self.near:
                s, factor, prefix = self.near[f]
                lo, hi = min(max(0, k-s), len(prefix)-1), min(k+1, len(prefix)-1)
                result = min(result, factor*(prefix[hi]-prefix[lo]))
        return result


def outer_prefix(context, weight):
    witness = context['witness']
    alpha, z = map(base.rational, (fraction(witness['alpha']), fraction(witness['z'])))
    q = context['field_order']
    density = tail.log2arb(base.rational(base.outer_density(q, context['outer_measure'])))
    return alpha*(weight*density-tail.log2arb(arb(q-1))
                  -context['forbidden_weight']*tail.log2arb(z))


def _log_record(record, precision):
    require(isinstance(record, dict), 'Missing support record')
    if record.get('zero') is True:
        same(record, dict(zero=True, log2_upper_rational=None, log2_upper_enclosure=None), 'zero support')
        return None
    same(record, dict(zero=False), 'positive support')
    index.scalar_ceiling(record, precision)
    return arb(record['log2_upper_enclosure'])


def audit_rows(row, context, *, pairs=4096, precision=384):
    """Generic small-pair audit for tests; public resolved bands use4096."""
    weight = row['active_groups']
    cutoff = context['count_provenance']['maximum_refined_singletons']
    counts = IndependentCounts(pairs, weight, cutoff)
    producer_counts = band.selected.near.NearFaceCounts(pairs, weight, max_singletons=cutoff)
    full, rows = counts.full, row['rows']
    same(row, dict(completed_even_rows=len(full), completed_even_full_range=[full.start, full.stop-1],
                   pending_even_full_range=None, expected_even_rows=len(full), expected_profile_pairs=len(full)**2,
                   profile_pairs_evaluated=len(full)**2), 'selected-weight coverage')
    require(len(rows) == len(full), 'Missing or duplicated even rows')
    total, rounded, positives, count_hash = arb(0), arb(0), 0, sha256()
    for first, saved in zip(full, rows):
        nonzero = 0
        for second in full:
            value = counts.count(first, second)
            require(value == producer_counts.count(first, second), 'Independent and producer exact profile counts differ')
            nonzero += bool(value)
            encoded = value.to_bytes((value.bit_length()+7)//8, 'big')
            count_hash.update(len(encoded).to_bytes(4, 'big'))
            count_hash.update(encoded)
        same(saved, dict(even_full=first, odd_full_range=[full.start, full.stop-1],
                         profile_pairs_evaluated=len(full), profiles_with_positive_count=nonzero,
                         complete_even_profile_row=True, covers_all_odd_profiles_for_this_even_profile=True,
                         complete_selected_weight=False, covers_all_supports_at_weight=False,
                         common_potential_used=False), 'even-row coverage/counts')
        logarithm = _log_record(saved['row_support'], precision)
        if logarithm is not None:
            total += arb(2)**logarithm.upper()
            rounded += arb(2)**base.rational(fraction(saved['row_support']['log2_upper_rational']))
        else:
            require(nonzero == 0, 'Zero row with positive support profiles')
        cumulative = _log_record(saved['completed_support'], precision)
        if total.is_zero():
            require(cumulative is None, 'Nonzero cumulative support for zero rows')
        else:
            require(cumulative is not None, 'Missing cumulative support')
            require(base.rational(fraction(saved['completed_support']['log2_upper_rational']))
                    >= tail.log2arb(total).upper(), 'Cumulative ceiling misses saved row enclosures')
        positives += nonzero
    same(row, dict(profiles_with_positive_count=positives), 'positive profile total')
    require(not total.is_zero(), 'This auditor requires a positive selected-weight support sum')
    support = _log_record(row['support_sum'], precision)
    require(support is not None, 'Missing completed support sum')
    require(base.rational(fraction(row['support_sum']['log2_upper_rational'])) >= tail.log2arb(total).upper(),
            'Completed support ceiling misses its saved rows')
    prefix = outer_prefix(context, weight)
    require(arb(row['log2_outer_prefix_enclosure']).overlaps(prefix), 'Saved outer prefix differs from fresh prefix')
    bound = index.scalar_ceiling(row, precision)
    value = prefix+tail.log2arb(total)
    require(base.rational(bound) >= value.upper(), 'Selected-weight ceiling misses fresh prefix plus saved rows')
    require(arb(row['log2_upper_enclosure']).overlaps(prefix+support), 'Summary and support/prefix enclosures disagree')
    return dict(active_groups=weight, even_rows=len(rows), profile_pairs=len(full)**2,
                profiles_with_positive_count=positives, exact_count_table_sha256=count_hash.hexdigest(),
                count_hash_encoding='row-major unsigned big-endian counts with4-byte big-endian byte lengths',
                saved_log2_upper_rational=str(bound),
                row_enclosure_reaggregation_log2_upper_rational=str(base.rational_ceiling(value)),
                rounded_row_reaggregation_log2_upper_rational=str(base.rational_ceiling(prefix+tail.log2arb(rounded))),
                fresh_outer_prefix=True, saved_ceiling_dominates_row_enclosures=True)


def _verify_witness(context):
    same(context, EXPECTED, 'gamma30/d40 context')
    same(context, dict(schema='gamma-near-face-band-context-v1', repetitions=256, paired_region_power=128,
                       common_potential_used=False, no_weight_interpolation=True,
                       complete_requested_range=False, covers_all_supports_in_range=False), 'context scope')
    witness = context['witness']
    require(0 < fraction(witness['z']) < 1 and 0 < fraction(witness['alpha']) <= 1, 'Invalid exact marks')
    markers = witness['markers']
    require(isinstance(markers, list) and 1 <= len(markers) <= 128
            and all(isinstance(pair, list) and len(pair) == 2 and all(fraction(x) > 0 for x in pair)
                    for pair in markers), 'Invalid exact positive marker bank')
    same(witness, dict(state_dimension=30, inner_width=64, gamma_extension_degree=16,
                       gamma_base_field_order=256, joint_rank_return_cap=True, all_pair_before_alpha=True), 'witness scope')
    require(type(witness['precision']) is int and 128 <= witness['precision'] <= 4096, 'Invalid bank precision')
    require(type(witness['packet_mantissa_bits']) is int and 32 <= witness['packet_mantissa_bits'] <= 512,
            'Invalid packet mantissa bits')
    require(type(witness['grid']) is int and witness['grid'] in (64, 128, 256), 'Invalid grid')
    cutoff = context['count_provenance']['maximum_refined_singletons']
    require(context['count_provenance'] == band.selected.count_provenance(cutoff), 'Unknown count provenance')
    same(context, dict(source_sha256=band.source_pins(), coverage_driver_sha256=band.source_pins(),
                       bank_source_sha256=band.provider.source_pins()), 'producer closures')


def audit_band(path, *, precision=384):
    require(type(precision) is int and 384 <= precision <= 4096, 'Audit precision must be384..4096')
    own_pins = source_pins()
    loader = index.Validator(F(2, 5))
    summary, summary_pin = loader.load(path)
    same(summary, dict(schema='gamma-near-face-band-summary-v1', construction=EXPECTED['construction'],
                       complete_requested_range=True, covers_all_supports_in_range=True,
                       every_integer_weight_evaluated=True, covers_all_even_odd_splits=True,
                       pending_weight_ranges=[], partial_weight_receipt=None, failure=None,
                       directed=True, whole_code_certificate=False, regional_caches_built=1,
                       common_actual_bank=True, no_weight_interpolation=True,
                       source_pins_unchanged_at_exit=True), 'complete compact band')
    context, context_pin = loader.dependency(summary['context'], summary_pin, 'shared-band-context')
    _verify_witness(context)
    lo, hi = index.interval(summary['requested_weight_range'])
    same(context, dict(requested_weight_range=[lo, hi]), 'context band range')
    same(summary, dict(weights_evaluated=hi-lo+1, completed_weight_range=[lo, hi]), 'band count')
    links = summary['weight_receipts']
    require(isinstance(links, list) and len(links) == hi-lo+1, 'Incomplete weight links')
    witness_source = context.get('witness_source')
    if witness_source:
        # The producer defines relative witness provenance against the theory directory.
        witness_path = Path(witness_source['path'])
        loader.load(witness_path if witness_path.is_absolute() else DIRECTORY/witness_path,
                    expected=witness_source['sha256'], role='literal-witness-provenance')
    reports, bounds = [], []
    with ctx.workprec(precision):
        for weight, link in zip(range(lo, hi+1), links):
            row, row_pin = loader.dependency(link, summary_pin, 'complete-weight-rows')
            require(row['context'] == summary['context'], 'Weight context differs from shared context')
            same(row, dict(active_groups=weight, status='complete', complete_selected_weight=True,
                           covers_all_supports_at_weight=True, covers_all_even_odd_splits=True,
                           common_potential_used=False, directed=True, whole_code_certificate=False,
                           failure=None, requires_linked_context=True, source_pins_checked_before_publication=True,
                           count_provenance=context['count_provenance'],
                           requested_failure_bits=context['requested_failure_bits']), 'weight scope')
            same(link, {key: row[key] for key in ('active_groups', 'status', 'complete_selected_weight',
                       'log2_upper_rational', 'completed_even_rows', 'pending_even_full_range')}, 'weight link')
            report = audit_rows(row, context, precision=precision)
            bound = fraction(row['log2_upper_rational'])
            same(row, dict(meets_requested_selected_weight_target=bound <= -context['requested_failure_bits']),
                 'weight target flag')
            reports.append(dict(report, receipt=row_pin))
            bounds.append(bound)
        union = index.union_ceiling(bounds, precision)
        saved = fraction(summary['log2_upper_rational'])
        require(saved >= union, 'Saved band ceiling misses directed384 union')
        same(summary, dict(completed_union_log2_upper_rational=str(saved),
                           requested_failure_bits=context['requested_failure_bits'],
                           requested_range_meets_target=saved <= -context['requested_failure_bits'],
                           weights_above_target=[lo+i for i, bound in enumerate(bounds)
                                                 if bound > -context['requested_failure_bits']]), 'band union/target')
    require(source_pins() == own_pins, 'Audit source closure changed')
    band.selected.check_source_pins(context['source_sha256'])
    band._check_links(tuple(loader.files.values()))
    return dict(schema='gamma-near-face-compact-band-audit-v1', audit_result='PASS',
                construction=EXPECTED['construction'], state_dimension=30, target_distance='2/5',
                weight_range=[lo, hi], summary=summary_pin, context=context_pin,
                weight_audits=reports, count_cutoff=context['count_provenance']['maximum_refined_singletons'],
                exact_marker_count=len(context['witness']['markers']), precision_bits=precision,
                saved_band_log2_upper_rational=str(saved), directed_union_log2_upper_rational=str(union),
                complete_range_event=True, coefficient_replay=False, whole_code_certificate=False,
                checked_file_pins=list(loader.files.values()), source_sha256=own_pins,
                checked_producer_source_sha256=context['source_sha256'],
                audit_limit='Exact count/coverage, linked bytes/sources, saved-row enclosures and fresh prefix only; '
                            'no coefficient/FFT replay or authentication of original execution')


def extend_v9(parent_path, audit, audit_pin, *, precision=384):
    """Add this one audited disjoint event; old receipt semantics remain inherited."""
    loader = index.Validator(F(2, 5))
    parent, parent_pin = loader.load(parent_path)
    same(parent, dict(schema='field128-gamma-partial-cover-index-v9', construction=EXPECTED['construction'],
                       field_order=EXPECTED['field_order'], source_dimension=EXPECTED['source_dimension'],
                       output_length=EXPECTED['output_length'], target_distance='2/5', whole_code_certificate=False),
         'parent index')
    same(audit, dict(schema='gamma-near-face-compact-band-audit-v1', audit_result='PASS',
                      construction=EXPECTED['construction'], complete_range_event=True), 'new audit')
    loaded_audit, unused = loader.load(audit_pin['path'], expected=audit_pin['sha256'], role='new-band-audit')
    require(loaded_audit == audit, 'Audit differs from its immutable file')
    band._check_links(audit['checked_file_pins'])
    resolved, unused = loader.load(audit['summary']['path'], expected=audit['summary']['sha256'], role='new-band-summary')
    same(resolved, dict(complete_requested_range=True, covers_all_supports_in_range=True,
                        requested_weight_range=audit['weight_range'], context=audit['context'],
                        log2_upper_rational=audit['saved_band_log2_upper_rational']), 'audit/summary binding')
    pins = dict(parent['aggregation_source_sha256'])
    for name, digest in audit['source_sha256'].items():
        require(name not in pins or pins[name] == digest, 'Conflicting source pins: '+name)
        pins[name] = digest
    band.selected.check_source_pins(pins)
    entries, covered = deepcopy(parent['entries']), set()
    for entry in entries:
        path = Path(entry['path'])
        loader.load(path if path.is_absolute() else DIRECTORY/path, expected=entry['sha256'], role='prior-event-bytes')
        lo, hi = index.interval(entry['weight_range'])
        added = set(range(lo, hi+1))
        require(not covered & added, 'Prior events overlap')
        covered |= added
    require(index.ranges(covered) == parent['covered_weight_ranges'], 'Prior ranges disagree with events')
    require(index.ranges(set(range(1, 8193))-covered) == parent['pending_weight_ranges'], 'Prior pending ranges disagree')
    require(len(entries) == parent['receipt_count'], 'Prior receipt count mismatch')
    old_union = index.union_ceiling([fraction(entry['log2_upper_rational']) for entry in entries], precision)
    require(old_union <= fraction(parent['covered_events_union_log2_upper_rational']), 'Prior union is understated')
    lo, hi = audit['weight_range']
    new_weights = set(range(lo, hi+1))
    require(not covered & new_weights, 'New event overlaps a covered weight')
    new_bound = fraction(audit['saved_band_log2_upper_rational'])
    require(new_bound <= -54, 'New event does not meet the requested54-bit admission target')
    summary = audit['summary']
    entries.append(dict(path=summary['path'], sha256=summary['sha256'], weight_range=[lo, hi],
                        log2_upper_rational=str(new_bound), receipt_kind='complete-near-face-compact-band',
                        weights_evaluated=hi-lo+1, context=audit['context'], audit=audit_pin))
    entries.sort(key=lambda entry: entry['weight_range'])
    covered |= new_weights
    with ctx.workprec(precision):
        value = tail.log2arb(sum((arb(2)**base.rational(fraction(entry['log2_upper_rational']))
                                 for entry in entries), arb(0)))
        union, enclosure = str(base.rational_ceiling(value)), str(value)
    result = deepcopy(parent)
    result.update(schema='field128-gamma-partial-cover-index-v10', entries=entries,
                  parent_index=dict(parent_pin, schema=parent['schema'],
                                    role='Immutable v9; prior event bytes and source closure rechecked; prior semantic audits inherited'),
                  covered_weight_ranges=index.ranges(covered), pending_weight_ranges=index.ranges(set(range(1, 8193))-covered),
                  covered_weight_count=len(covered), pending_weight_count=8192-len(covered), receipt_count=len(entries),
                  aggregation_precision_bits=precision, aggregation_source_sha256=dict(sorted(pins.items())),
                  covered_events_union_log2_upper_rational=union, covered_events_union_log2_upper_enclosure=enclosure,
                  added_event_count=1, added_weight_count=hi-lo+1, added_events_union_log2_upper_rational=str(new_bound),
                  latest_compact_band_audit=audit_pin, compact_context_resolution_explicit=True,
                  all_receipt_bytes_match=True, all_helper_source_closures_match=True, source_pin_conflicts=False,
                  input_file_pins_verified=len(loader.files)+len(audit['checked_file_pins']),
                  helper_source_files_verified=len(pins),
                  aggregation_method=f'Directed Arb union over{len(entries)} disjoint events; terminal counted once; '
                                     'one new complete context-resolved integer band, no interpolation',
                  scope=f'Only the listed{len(covered)} weights are covered; all{8192-len(covered)} other weights remain pending. '
                        'This index is not a whole-code certificate.')
    result['producer_receipt_counts']['hd_field128_gamma_near_face_band'] = 1
    band._check_links(tuple(loader.files.values()))
    band.selected.check_source_pins(pins)
    return result


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--summary', type=Path, required=True)
    parser.add_argument('--audit-output', type=Path, required=True)
    parser.add_argument('--parent-v9', type=Path)
    parser.add_argument('--index-output', type=Path)
    args = parser.parse_args()
    require(bool(args.parent_v9) == bool(args.index_output), 'Supply both parent-v9 and index-output')
    require(not args.audit_output.exists() and (not args.index_output or not args.index_output.exists()),
            'Outputs must be unused')
    audit = audit_band(args.summary)
    pin = band.write_new(args.audit_output, audit)
    print(json.dumps(dict(audit=pin, result=audit['audit_result'], weight_range=audit['weight_range'],
                         union=audit['directed_union_log2_upper_rational'])), flush=True)
    if args.parent_v9:
        result = extend_v9(args.parent_v9, audit, pin)
        output = band.write_new(args.index_output, result)
        print(json.dumps(dict(index=output, covered=result['covered_weight_count'], pending=result['pending_weight_count'],
                              events=result['receipt_count'], union=result['covered_events_union_log2_upper_rational'])), flush=True)


if __name__ == '__main__':
    main()
