"""Extend frozen v25 with audited913 and optional typed-audited bands.

Counts are inherited from the explicit independent audits. Frozen validators
recheck the selected rows/prefix and band partitions/unions at384 bits. This
module returns an object only; it never writes an index or generates caps.
"""
from copy import deepcopy
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

import hd_field128_gamma_selected_index_v23 as bands
import hd_field128_gamma_selected_index_v14 as selected
import hd_field128_gamma30_d40_r913_portable_audit_v25 as portable_audit

index, midpoint, band = bands.index, selected.midpoint, bands.band
require, same, fraction = bands.require, bands.same, bands.fraction
DIRECTORY = Path(__file__).resolve().parent
PARENT_SHA256 = 'c92b67290ff9b1a729496598e071a7599a3c76ce2486e6441023ab3dee59aea0'
PARENT_EVENT_COUNT, PARENT_WEIGHT_COUNT = 249, 7744
BANDS_SHA256 = '5543d4984d341204cc8de168a03a18cc1d7dbd700e3890debe097de7305808c8'
SELECTED_SHA256 = 'd3006bd967827e1facefdac72d71088744a007e8b71d62c8398fff2345f68784'
PORTABLE_AUDIT_SHA256 = 'ae2c4364700d2d2a44bed54ec0eabed7840403bc2689594426b58e6fbf9ae96a'


def source_pins():
    require(sha256(Path(bands.__file__).read_bytes()).hexdigest() == BANDS_SHA256,
            'Frozen typed band validator changed')
    require(sha256(Path(selected.__file__).read_bytes()).hexdigest() == SELECTED_SHA256,
            'Frozen selected receipt validator changed')
    require(sha256(Path(portable_audit.__file__).read_bytes()).hexdigest() == PORTABLE_AUDIT_SHA256,
            'Reviewed portable selected audit adapter changed')
    return dict(sorted({**index.source_closure(Path(__file__).name),
                        **bands.source_pins(), **midpoint.source_pins(),
                        **portable_audit.source_pins()}.items()))


def _local(path):
    path = Path(path)
    return (path if path.is_absolute() else DIRECTORY/path).resolve()


def _midpoint_event(audit_path, loader, pins, source_inputs):
    """Only original913; portable links resolve at theory, never by rebasing JSON."""
    report, audit_pin = loader.load(audit_path, role='independent-midpoint913-audit')
    same(report, dict(schema='gamma-midpoint-selected-receipts-audit-v1', audit_result='PASS',
        construction=midpoint.previous.EXPECTED['construction'], state_dimension=30,
        target_distance='2/5', selected_weights=[913], precision_bits=384,
        source_sha256=midpoint.source_pins(), coefficient_replay=False,
        whole_code_certificate=False, external_checkpoint_context_reconstructed=True,
        external_rows_checked=457), 'independent913 midpoint audit')
    adapter_pins = portable_audit.source_pins()
    same(report, dict(audit_adapter_schema='gamma-r913-portable-midpoint-audit-adapter-v1',
        adapter=dict(path=portable_audit.NAME, sha256=PORTABLE_AUDIT_SHA256),
        adapter_source_sha256=adapter_pins, portable_proposal_links_resolved_from_theory=True,
        wrapper_limits_and_exact_proposal_checked=True, external_rows_equal_complete_inline_rows=True,
        source_and_input_bytes_checked_at_exit=True), 'portable audit adapter identity')
    source_inputs.load_source(Path(portable_audit.__file__), PORTABLE_AUDIT_SHA256)
    require(len(report['weight_audits']) == 1, 'Expected exactly one audited913 event')
    verified = report['weight_audits'][0]
    record, pin = loader.dependency(verified['receipt'], audit_pin, 'audited913-summary')
    bound = selected.complete_selected(record, 913, midpoint.producer.source_pins())
    same(record, dict(requested_failure_bits=54, meets_requested_selected_weight_target=True),
         'selected54-bit target')
    same(verified, dict(active_groups=913, even_rows=457, profile_pairs=208849,
        saved_log2_upper_rational=str(bound), all_producer_midpoint_and_final_counts_match=True,
        saved_ceiling_dominates_row_enclosures=True), 'selected count and row audit')
    # Counts are not recomputed here. The audited receipt pins the full exact
    # policy; all457 saved row enclosures and the original-r prefix were just
    # reaggregated by frozen complete_selected.
    checked = set()
    for reference in report['checked_file_pins']:
        unused, linked = loader.dependency(reference, audit_pin, 'independent-selected-audit-input')
        checked.add((linked['path'], linked['sha256']))
    require((pin['path'], pin['sha256']) in checked, 'Audit did not pin its selected receipt')
    external = {}
    for name, digest in report['external_input_sha256'].items():
        value, linked = loader.load(_local(name), expected=digest, role='external-midpoint-checkpoint')
        require(linked['path'] not in external, 'Duplicate resolved external checkpoint')
        external[linked['path']] = (value, linked)
    directory = Path(pin['path']).parent

    def checkpoint(name):
        path = str((directory/name).resolve())
        require(path in external, 'Missing external checkpoint: '+name)
        return external[path]

    unused, summary_pin = checkpoint('summary.json')
    require(summary_pin['sha256'] == pin['sha256'], 'External selected summary differs')
    context_record, unused = checkpoint('context.json')
    context = context_record['context']
    digest = sha256(json.dumps(context, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()
    require(digest == context_record['context_sha256'] == record['context_sha256'] == report['context_sha256'],
            'Selected context digest differs')
    same(record, {key:value for key,value in context.items() if key != 'scope'}, 'selected metadata context')
    for first, row in enumerate(record['rows']):
        saved, unused = checkpoint(f'row_{first:04d}.json')
        same(saved, row, 'external complete selected row')
        same(saved, dict(active_groups=913, context_sha256=digest,
                         directed=True, whole_code_certificate=False), 'external row context')
    proposal, proposal_pin = loader.load(_local(record['witness_source']['path']),
        expected=record['witness_source']['sha256'], role='portable-exact913-proposal')
    audited_proposal, verified_pin = loader.dependency(verified['proposal'], audit_pin, 'audited913-proposal')
    require(proposal_pin['path'] == verified_pin['path'] and proposal_pin['sha256'] == verified_pin['sha256'],
            'Audit proposal differs from original portable literal link')
    same(record['witness'], proposal['witness'], 'exact913 proposal literals')
    require((proposal_pin['path'], proposal_pin['sha256']) in checked, 'Audit did not pin exact proposal')
    bank_record, unused = checkpoint('bank.json')
    same(bank_record['witness'], record['witness'], 'bank exact913 witness')
    driver, unused = checkpoint('driver.json')
    same(driver, dict(active_groups=913, marker_count=83, expected_even_rows=457,
        expected_profile_pairs=208849, requested_failure_bits=54,
        fresh_bank_required=True, fresh_bank_rebuilt=False,
        old_numerical_bounds_reused=False, directed=False, whole_code_certificate=False), 'execution scope')
    same(bank_record, dict(fresh_bank_rebuilt=True), 'fresh execution bank')
    require(driver['driver'] == record['execution_driver'], 'Execution driver identity differs')
    require(report['execution_source_sha256'] == record['execution_source_sha256']
            == driver['source_sha256'] == bank_record['runner_source_sha256'], 'Execution source closure differs')
    execution = record['execution_driver']
    source_inputs.load_source(_local(execution['path']), execution['sha256'])
    name = Path(execution['path']).name
    require(record['execution_source_sha256'].get(name) == execution['sha256'], 'Driver absent from execution closure')
    for closure in (report['source_sha256'], adapter_pins, report['execution_source_sha256'],
                    record['source_sha256'], record['bank_source_sha256']):
        bands.frozen._merge(pins, closure)
    return dict(pin, weight_range=[913,913], weights_evaluated=1,
        log2_upper_rational=str(bound), receipt_kind='audited-midpoint-selected',
        audit=audit_pin, context_sha256=digest)


def extend_v25(parent_path, midpoint_audit, band_audits=()):
    """One mandatory913 point plus0..63 explicit bands; return, never publish."""
    require(isinstance(band_audits,(list,tuple)) and len(band_audits) <= 63,
            'Expected zero to63 explicit typed bands plus one selected point')
    loader, pins = index.Validator(F(2,5)), source_pins()
    source_inputs = bands.previous.audited.Inputs()
    parent, parent_pin = loader.load(parent_path, expected=PARENT_SHA256, role='frozen-v25-index')
    same(parent, dict(schema='field128-gamma-partial-cover-index-v23',
        construction=midpoint.previous.EXPECTED['construction'], field_order=1<<128,
        source_dimension=1<<20, output_length=1<<21, target_distance='2/5',
        whole_code_certificate=False, directed_aggregation=True,
        receipt_count=PARENT_EVENT_COUNT, covered_weight_count=PARENT_WEIGHT_COUNT,
        pending_weight_count=8192-PARENT_WEIGHT_COUNT), 'frozen v25 checkpoint')
    require(len(parent['entries']) == PARENT_EVENT_COUNT, 'Unexpected parent event count')
    bands.frozen._merge(pins,parent['aggregation_source_sha256'])
    entries, covered = deepcopy(parent['entries']), set()
    for entry in entries:
        loader.load(_local(entry['path']), expected=entry['sha256'], role='prior-event-bytes')
        if 'audit' in entry:
            loader.dependency(entry['audit'],parent_pin,'inherited-prior-audit')
        lo,hi=index.interval(entry['weight_range']); weights=set(range(lo,hi+1))
        require(not covered & weights,'Prior event overlap');covered |= weights
    same(parent,dict(covered_weight_ranges=index.ranges(covered),covered_weight_count=len(covered),
        pending_weight_ranges=index.ranges(set(range(1,8193))-covered),
        pending_weight_count=8192-len(covered),receipt_count=len(entries)),'prior coverage')
    require(index.union_ceiling(bands.frozen._bounds(entries),384)
            <= fraction(parent['covered_events_union_log2_upper_rational']),'Prior union understated')
    point = _midpoint_event(midpoint_audit,loader,pins,source_inputs)
    added=[point]
    for path in band_audits:
        added.append(bands._band_event(path,loader,pins,source_inputs))
    for entry in added:
        lo,hi=index.interval(entry['weight_range']); weights=set(range(lo,hi+1))
        require(not covered & weights,'New event overlaps prior or newly added coverage')
        covered |= weights;entries.append(entry)
    require(entries[:PARENT_EVENT_COUNT] == parent['entries'],'Parent entry prefix changed')
    union=index.union_ceiling(bands.frozen._bounds(entries),384)
    require(union <= -54,'Covered-event union misses54-bit target')
    result=deepcopy(parent)
    result.pop('covered_events_union_log2_upper_enclosure',None)
    result.update(schema='field128-gamma-partial-cover-index-v26',parent_index=parent_pin,
        entries=entries,receipt_count=len(entries),whole_code_certificate=False,
        covered_weight_ranges=index.ranges(covered),covered_weight_count=len(covered),
        pending_weight_ranges=index.ranges(set(range(1,8193))-covered),pending_weight_count=8192-len(covered),
        aggregation_precision_bits=384,covered_events_union_log2_upper_rational=str(union),
        covered_events_failure_bits=54,meets_covered_events_target=True,
        aggregation_source_sha256=dict(sorted(pins.items())),added_event_count=len(added),
        added_weight_count=sum(row['weights_evaluated'] for row in added),added_events_zero=False,
        added_events_union_log2_upper_rational=str(index.union_ceiling(bands.frozen._bounds(added),384)),
        latest_added_events=added,latest_explicit_selected_events=[point],latest_audited_multiweight_events=added[1:],
        prior_entry_records_and_order_preserved=True,prior_semantic_audits_inherited=True,
        new_count_audits_inherited=True,new_selected_rows_and_prefix_reaggregated=True,
        new_integer_partitions_rechecked=True,new_cell_unions_reaggregated=True,
        coefficient_replay=False,normalized_regional_replay=False,transfer_recomputed=False,
        input_file_pins_verified=len(loader.files),non_json_input_file_pins_verified=len(source_inputs.non_json_pins),
        checked_non_json_file_pins=list(source_inputs.non_json_pins.values()),
        new_audit_schemas_allowed=[*bands.AUDIT_SCHEMAS,'gamma-midpoint-selected-receipts-audit-v1'],
        new_selected_audit_adapter_schemas_allowed=['gamma-r913-portable-midpoint-audit-adapter-v1'],
        helper_source_files_verified=len(pins),all_receipt_bytes_match=True,
        all_helper_source_closures_match=True,source_pin_conflicts=False,no_weight_interpolation=True,
        aggregation_method=f'Directed384-bit union of {len(entries)} events on nonoverlapping weight ranges',
        scope=f'Only {len(covered)} listed weights; {8192-len(covered)} remain pending; no whole-code certificate')
    counts=result['producer_receipt_counts']
    counts['hd_field128_gamma_midpoint_certificate']=counts.get('hd_field128_gamma_midpoint_certificate',0)+1
    name='hd_field128_gamma_midpoint_weight_bands';counts[name]=counts.get(name,0)+len(band_audits)
    band._check_links(tuple(loader.files.values()))
    source_inputs.check()
    bands.frozen.producer.check_source_pins(pins)
    return result
