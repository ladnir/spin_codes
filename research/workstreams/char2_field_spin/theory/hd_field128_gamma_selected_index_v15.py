"""Add only the externally audited r3072 midpoint event to frozen v14."""

from copy import deepcopy
from fractions import Fraction as F
from pathlib import Path

import hd_field128_gamma_selected_index_v14 as frozen

index, band, midpoint = frozen.index, frozen.band, frozen.midpoint
require, same, fraction = frozen.require, frozen.same, frozen.fraction


def source_pins():
    return dict(sorted({**index.source_closure(Path(__file__).name), **midpoint.source_pins()}.items()))


def extend_v14(parent_path, audit_path):
    loader, pins = index.Validator(F(2, 5)), source_pins()
    parent, parent_pin = loader.load(parent_path)
    same(parent, dict(schema='field128-gamma-partial-cover-index-v14', construction=midpoint.previous.EXPECTED['construction'],
                      field_order=1 << 128, source_dimension=1 << 20, output_length=1 << 21,
                      target_distance='2/5', whole_code_certificate=False, directed_aggregation=True), 'frozen v14')
    frozen.prior._merge(pins, parent['aggregation_source_sha256'])
    entries, covered = deepcopy(parent['entries']), set()
    for entry in entries:
        loader.load(entry['path'], expected=entry['sha256'], role='prior-event-bytes')
        if 'audit' in entry:
            loader.dependency(entry['audit'], parent_pin, 'inherited-prior-audit')
        lo, hi = index.interval(entry['weight_range'])
        weights = set(range(lo, hi+1))
        require(not covered & weights, 'Prior event overlap')
        covered |= weights
    same(parent, dict(covered_weight_ranges=index.ranges(covered), covered_weight_count=len(covered),
                      pending_weight_ranges=index.ranges(set(range(1, 8193))-covered), pending_weight_count=8192-len(covered),
                      receipt_count=len(entries)), 'prior coverage')
    require(index.union_ceiling([fraction(x['log2_upper_rational']) for x in entries], 384)
            <= fraction(parent['covered_events_union_log2_upper_rational']), 'Prior union understated')
    report, audit_pin = loader.load(audit_path, role='inherited-r3072-midpoint-audit')
    same(report, dict(schema='gamma-midpoint-selected-receipts-audit-v1', audit_result='PASS', selected_weights=[3072],
                      construction=midpoint.previous.EXPECTED['construction'], state_dimension=30, target_distance='2/5',
                      precision_bits=384, source_sha256=midpoint.source_pins(), coefficient_replay=False,
                      whole_code_certificate=False, external_checkpoint_context_reconstructed=True, external_rows_checked=1537,
                      exact_positive_marker_order_preserved=True, source_positive_proposals_match_float_artifact=True), 'r3072 audit')
    require(len(report['weight_audits']) == 1, 'Expected exactly one audited selected weight')
    verified = report['weight_audits'][0]
    record, pin = loader.dependency(verified['receipt'], audit_pin, 'audited-r3072-receipt')
    bound = frozen.complete_selected(record, 3072, midpoint.producer.source_pins())
    same(verified, dict(active_groups=3072, even_rows=1537, profile_pairs=2362369, saved_log2_upper_rational=str(bound),
                        all_producer_midpoint_and_final_counts_match=True, saved_ceiling_dominates_row_enclosures=True), 'count/arithmetic audit')
    proposal, unused = loader.dependency(verified['proposal'], audit_pin, 'r3072-exact-proposal')
    same(record['witness'], proposal['witness'], 'exact proposal bank literals')
    same(report, dict(context_sha256=record['context_sha256'], exact_positive_marker_count=len(record['witness']['markers'])), 'context binding')
    external = report['external_input_sha256']
    directory = Path(report['receipt_directory'])
    required = [directory/name for name in ('summary.json', 'proposal.json', 'context.json', 'bank.json')]
    required += [directory/'rows'/f'even_{a:04d}.json' for a in range(1537)]
    require(all(str(path) in external for path in required), 'Incomplete external context/checkpoint pins')
    for path, digest in external.items():
        loader.load(path, expected=digest, role='audited-external-checkpoint')
    require(external[str(Path(pin['path']))] == pin['sha256'], 'Audit external summary differs')
    band._check_links(report['checked_file_pins'])
    frozen.prior._merge(pins, report['source_sha256'])
    frozen.prior._merge(pins, report['additional_source_sha256'])
    require(3072 not in covered, 'r3072 overlaps prior coverage')
    added = dict(pin, weight_range=[3072, 3072], log2_upper_rational=str(bound),
                 receipt_kind='audited-midpoint-selected', audit=audit_pin)
    entries.append(added)
    entries.sort(key=lambda row: row['weight_range'])
    covered.add(3072)
    union = index.union_ceiling([fraction(row['log2_upper_rational']) for row in entries], 384)
    result = deepcopy(parent)
    result.pop('covered_events_union_log2_upper_enclosure', None)
    result.update(schema='field128-gamma-partial-cover-index-v15', parent_index=parent_pin, entries=entries,
                  receipt_count=len(entries), covered_weight_ranges=index.ranges(covered), covered_weight_count=len(covered),
                  pending_weight_ranges=index.ranges(set(range(1, 8193))-covered), pending_weight_count=8192-len(covered),
                  aggregation_precision_bits=384, covered_events_union_log2_upper_rational=str(union),
                  aggregation_source_sha256=dict(sorted(pins.items())), added_event_count=1, added_weight_count=1,
                  added_events_union_log2_upper_rational=str(bound), latest_explicit_selected_events=[added],
                  prior_semantic_audits_inherited=True, new_semantic_audit_inherited=True, coefficient_replay=False,
                  input_file_pins_verified=len(loader.files), helper_source_files_verified=len(pins),
                  all_receipt_bytes_match=True, all_helper_source_closures_match=True, source_pin_conflicts=False,
                  no_weight_interpolation=True, excluded_overlapping_new_weights=[768],
                  aggregation_method=f'Directed384-bit union of {len(entries)} disjoint events; add only selected r3072',
                  scope=f'Only {len(covered)} listed weights; {8192-len(covered)} remain pending; no whole-code certificate')
    counts = result['producer_receipt_counts']
    counts['hd_field128_gamma_midpoint_certificate'] = counts.get('hd_field128_gamma_midpoint_certificate', 0)+1
    band._check_links(tuple(loader.files.values()))
    band.selected.check_source_pins(pins)
    return result
