"""Append typed-audited bands to byte-pinned v22; never publish automatically.

Prior event records retain their contents and order. New counts are inherited
from explicit typed audits; saved partitions and unions use frozen validators.
Even complete integer coverage remains a partial index, not a whole-code proof.
"""

from copy import deepcopy
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path

import hd_field128_gamma_selected_index_v22_typed as previous


frozen, index, band, midpoint = previous.frozen, previous.index, previous.band, previous.midpoint
require, same, fraction = previous.require, previous.same, previous.fraction
DIRECTORY = Path(__file__).resolve().parent
PARENT_SHA256 = 'e160fb6865347307ed3e4ccc41b95fc86551ddd7e39f5ac470ea5bef92d54e4d'
PARENT_EVENT_COUNT, PARENT_WEIGHT_COUNT = 242, 7599
PREVIOUS_SHA256 = 'cab5044bb2735f2a7b03d399a0a79c693674f7c2cf136265c18500851fa873f2'
FAILURE_BITS = 54
AUDIT_SCHEMAS = ('gamma-midpoint-weight-band-typed-independent-audit-v1',
                 'gamma-midpoint-weight-band-typed-independent-audit-v2')


def source_pins():
    require(sha256(Path(previous.__file__).read_bytes()).hexdigest() == PREVIOUS_SHA256,
            'Frozen v22 typed validator changed')
    return dict(sorted({**index.source_closure(Path(__file__).name),
                        **previous.source_pins()}.items()))


def _band_event(path, loader, pins, source_inputs):
    report, unused = loader.load(path, role='explicit-v23-typed-band-audit')
    require(report.get('schema') in AUDIT_SCHEMAS, 'Expected an explicit typed v1/v2 band audit')
    return previous._band_event(path, loader, pins, source_inputs)


def extend_v22(parent_path, audit_paths):
    """Return a checked v23 object for 1..64 new bands; do not write any file."""
    require(isinstance(audit_paths, (list, tuple)) and 1 <= len(audit_paths) <= 64,
            'Expected one to 64 explicit audited bands')
    loader, pins = index.Validator(F(2, 5)), source_pins()
    source_inputs = previous.audited.Inputs()
    parent, parent_pin = loader.load(parent_path, expected=PARENT_SHA256, role='frozen-v22-typed-index')
    same(parent, dict(schema='field128-gamma-partial-cover-index-v22-typed',
        construction=midpoint.previous.EXPECTED['construction'], field_order=1 << 128,
        source_dimension=1 << 20, output_length=1 << 21, target_distance='2/5',
        whole_code_certificate=False, directed_aggregation=True,
        receipt_count=PARENT_EVENT_COUNT, covered_weight_count=PARENT_WEIGHT_COUNT,
        pending_weight_count=8192-PARENT_WEIGHT_COUNT), 'frozen v22 typed index')
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
    # Do not sort: the exact prior records remain the unchanged entry prefix.
    require(entries[:PARENT_EVENT_COUNT] == parent['entries'], 'Prior event records changed')
    union = index.union_ceiling(frozen._bounds(entries), 384)
    require(union <= -FAILURE_BITS, 'Covered-event union fails 54-bit budget')
    new_bounds = frozen._bounds(added)
    result = deepcopy(parent)
    result.update(schema='field128-gamma-partial-cover-index-v23', parent_index=parent_pin, entries=entries,
        whole_code_certificate=False, receipt_count=len(entries),
        covered_weight_ranges=index.ranges(covered), covered_weight_count=len(covered),
        pending_weight_ranges=index.ranges(set(range(1, 8193))-covered), pending_weight_count=8192-len(covered),
        aggregation_precision_bits=384, covered_events_union_log2_upper_rational=str(union),
        covered_events_failure_bits=FAILURE_BITS, meets_covered_events_target=True,
        aggregation_source_sha256=dict(sorted(pins.items())), added_event_count=len(added),
        added_weight_count=sum(row['weights_evaluated'] for row in added),
        added_events_zero=not new_bounds,
        added_events_union_log2_upper_rational=str(index.union_ceiling(new_bounds, 384)) if new_bounds else None,
        latest_audited_multiweight_events=added,
        prior_entry_records_and_order_preserved=True,
        prior_semantic_audits_inherited=True, new_count_audits_inherited=True,
        new_integer_partitions_rechecked=True, new_cell_unions_reaggregated=True,
        coefficient_replay=False, normalized_regional_replay=False,
        input_file_pins_verified=len(loader.files),
        non_json_input_file_pins_verified=len(source_inputs.non_json_pins),
        checked_non_json_file_pins=list(source_inputs.non_json_pins.values()),
        new_audit_schemas_allowed=list(AUDIT_SCHEMAS), helper_source_files_verified=len(pins),
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
