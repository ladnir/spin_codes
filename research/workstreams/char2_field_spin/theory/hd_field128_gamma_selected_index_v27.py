"""Append explicit typed-audited bands to the byte-pinned v26 checkpoint.

Prior semantic audits are inherited. The frozen band validator rechecks new
integer partitions and saved cell unions. This module returns an object only;
it never writes an index, generates coefficients, or claims a whole-code proof.
"""
from copy import deepcopy
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path

import hd_field128_gamma_selected_index_v23 as bands

index, frozen, band, midpoint = bands.index, bands.frozen, bands.band, bands.midpoint
require, same, fraction = bands.require, bands.same, bands.fraction
DIRECTORY = Path(__file__).resolve().parent
PARENT_SHA256 = '45d2854157e47e23145f76df445bc8cae4f318ce828003e8322eb36307ea1f20'
PARENT_EVENT_COUNT, PARENT_WEIGHT_COUNT = 253, 7928
BANDS_SHA256 = '5543d4984d341204cc8de168a03a18cc1d7dbd700e3890debe097de7305808c8'
FAILURE_BITS = 54


def source_pins():
    require(sha256(Path(bands.__file__).read_bytes()).hexdigest() == BANDS_SHA256,
            'Frozen typed band validator changed')
    pins = index.source_closure(Path(__file__).name)
    frozen._merge(pins, bands.source_pins())
    return dict(sorted(pins.items()))


def _local(path):
    path = Path(path)
    return path if path.is_absolute() else DIRECTORY/path


def extend_v26(parent_path, audit_paths):
    """Return a checked v27 object for 1..64 new bands; never publish it."""
    require(isinstance(audit_paths, (list, tuple)) and 1 <= len(audit_paths) <= 64,
            'Expected one to 64 explicit typed band audits')
    loader, pins = index.Validator(F(2, 5)), source_pins()
    source_inputs = bands.previous.audited.Inputs()
    parent, parent_pin = loader.load(parent_path, expected=PARENT_SHA256, role='frozen-v26-index')
    same(parent, dict(schema='field128-gamma-partial-cover-index-v26',
        construction=midpoint.previous.EXPECTED['construction'], field_order=1 << 128,
        source_dimension=1 << 20, output_length=1 << 21, target_distance='2/5',
        whole_code_certificate=False, directed_aggregation=True,
        receipt_count=PARENT_EVENT_COUNT, covered_weight_count=PARENT_WEIGHT_COUNT,
        pending_weight_count=8192-PARENT_WEIGHT_COUNT), 'frozen v26 checkpoint')
    require(len(parent['entries']) == PARENT_EVENT_COUNT, 'Unexpected parent event count')
    frozen._merge(pins, parent['aggregation_source_sha256'])
    # Keep and recheck inherited source links rather than replacing their list
    # with only the links introduced by this extension.
    for link in parent.get('checked_non_json_file_pins', []):
        checked = source_inputs.load_source(_local(link['path']), link['sha256'])
        same(checked, link, 'inherited non-JSON source pin')
    entries, covered = deepcopy(parent['entries']), set()
    for entry in entries:
        loader.load(_local(entry['path']), expected=entry['sha256'], role='prior-event-bytes')
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
        entry = bands._band_event(path, loader, pins, source_inputs)
        lo, hi = index.interval(entry['weight_range'])
        weights = set(range(lo, hi+1))
        require(not covered & weights, 'New band overlaps prior or newly added coverage')
        covered |= weights
        entries.append(entry)
        added.append(entry)
    require(entries[:PARENT_EVENT_COUNT] == parent['entries'], 'Parent entry prefix changed')
    union = index.union_ceiling(frozen._bounds(entries), 384)
    require(union <= -FAILURE_BITS, 'Covered-event union fails 54-bit budget')
    new_bounds = frozen._bounds(added)
    result = deepcopy(parent)
    result.pop('covered_events_union_log2_upper_enclosure', None)
    result.update(schema='field128-gamma-partial-cover-index-v27', parent_index=parent_pin,
        entries=entries, receipt_count=len(entries), whole_code_certificate=False,
        covered_weight_ranges=index.ranges(covered), covered_weight_count=len(covered),
        pending_weight_ranges=index.ranges(set(range(1, 8193))-covered), pending_weight_count=8192-len(covered),
        aggregation_precision_bits=384, covered_events_union_log2_upper_rational=str(union),
        covered_events_failure_bits=FAILURE_BITS, meets_covered_events_target=True,
        aggregation_source_sha256=dict(sorted(pins.items())), added_event_count=len(added),
        added_weight_count=sum(row['weights_evaluated'] for row in added), added_events_zero=not new_bounds,
        added_events_union_log2_upper_rational=str(index.union_ceiling(new_bounds, 384)) if new_bounds else None,
        latest_added_events=added, latest_explicit_selected_events=[], latest_audited_multiweight_events=added,
        prior_entry_records_and_order_preserved=True, prior_semantic_audits_inherited=True,
        prior_non_json_source_pins_rechecked=True, new_count_audits_inherited=True,
        new_selected_rows_and_prefix_reaggregated=False,
        new_integer_partitions_rechecked=True, new_cell_unions_reaggregated=True,
        coefficient_replay=False, normalized_regional_replay=False, transfer_recomputed=False,
        input_file_pins_verified=len(loader.files), non_json_input_file_pins_verified=len(source_inputs.non_json_pins),
        checked_non_json_file_pins=list(source_inputs.non_json_pins.values()),
        new_audit_schemas_allowed=list(bands.AUDIT_SCHEMAS), new_selected_audit_adapter_schemas_allowed=[],
        helper_source_files_verified=len(pins), all_receipt_bytes_match=True,
        all_helper_source_closures_match=True, source_pin_conflicts=False, no_weight_interpolation=True,
        zero_event_semantics='zero=true and null log ceiling; omitted from positive union only',
        aggregation_method=f'Directed384-bit union of {len(entries)} events on nonoverlapping weight ranges',
        scope=f'Only {len(covered)} listed weights; {8192-len(covered)} remain pending; no whole-code certificate')
    counts = result['producer_receipt_counts']
    name = 'hd_field128_gamma_midpoint_weight_bands'
    counts[name] = counts.get(name, 0)+len(added)
    band._check_links(tuple(loader.files.values()))
    source_inputs.check()
    frozen.producer.check_source_pins(pins)
    return result
